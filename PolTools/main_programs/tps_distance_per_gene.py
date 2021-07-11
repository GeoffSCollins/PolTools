import multiprocessing
import sys
import argparse

from PolTools.utils.get_region_length import determine_region_length
from PolTools.utils.make_transcripts_dict import build_transcripts_dict
from PolTools.utils.print_tab_delimited import print_tab_delimited
from PolTools.utils.verify_region_length_is_even import verify_region_length_is_even


def get_pausing_distances_helper(region_filename, transcripts_dict, region_length):
    # Go through each gene and get the distances from each one
    ret_dict = {}
    with open(region_filename) as file:
        for line in file:
            chromosome, left, right, gene_name, fold_change, strand = line.rstrip().split()

            five_prime_position = int(left) + int(region_length / 2)

            # Get all of the transcript lengths at this position
            if five_prime_position in transcripts_dict[chromosome][strand]:
                curr_dict = transcripts_dict[chromosome][strand][five_prime_position]
                most_common_pausing_distance = abs(max(curr_dict, key=curr_dict.get) - five_prime_position)

                # If the strand is positive, we add 1 because the transcripts dict is inclusive
                if strand == "+":
                    most_common_pausing_distance += 1

            else:
                most_common_pausing_distance = "N/A"

            ret_dict[gene_name] = most_common_pausing_distance

    return ret_dict

def get_pausing_distances(sequencing_filename, region_filename, region_length):
    transcripts_dict = build_transcripts_dict(sequencing_filename)
    ret_dict = get_pausing_distances_helper(region_filename, transcripts_dict, region_length)

    return ret_dict, sequencing_filename


def parse_args(args):
    def positive_int(num):
        try:
            val = int(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")

        return val

    parser = argparse.ArgumentParser(prog='PolTools tps_distance_per_gene',
                                     description='Determine the most common pausing distance from the max TSS for each gene.\n' +
                                                 "More information can be found at " +
                                                 "https://geoffscollins.github.io/PolTools/tps_distance_per_gene.html")

    parser.add_argument('regions_filename', metavar='regions_filename', type=str,
                        help='Bed formatted regions file with an even region length or a region length of one.')

    parser.add_argument('seq_files', metavar='sequencing_files', nargs='+', type=str,
                        help='Bed formatted files from the sequencing experiment')

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    args = parser.parse_args(args)

    regions_filename = args.regions_filename
    sequencing_files_list = args.seq_files
    max_threads = args.threads

    if len(sequencing_files_list) == 0:
        sys.stderr.write("")

    region_length = determine_region_length(regions_filename)
    if region_length != 1:
        verify_region_length_is_even(region_length)

    return regions_filename, sequencing_files_list, region_length, max_threads


def output_data(pausing_distances):
    output_dict = {}
    for tup in pausing_distances:
        ret_dict, seq_filename = tup

        for gene_name in ret_dict:
            if gene_name not in output_dict:
                output_dict[gene_name] = {}

            output_dict[gene_name][seq_filename] = ret_dict[gene_name]


    for i, gene_name in enumerate(output_dict):
        if i == 0:
            # We print the headers
            print_tab_delimited(["Gene"] + [x.split("/")[-1] for x in output_dict[gene_name].keys()])

        print_tab_delimited([gene_name] + [output_dict[gene_name][sequencing_filename] for sequencing_filename in output_dict[gene_name]])


def run_tps_distance_per_gene(regions_filename, sequencing_files_list, region_length, max_threads):
    with multiprocessing.Pool(processes=max_threads) as pool:
        pausing_distances = pool.starmap(get_pausing_distances, [(sequencing_filename, regions_filename, region_length) for sequencing_filename in sequencing_files_list])

    output_data(pausing_distances)


def main(args):
    """
    PolTools tps_distance_per_gene <Regions Filename> <Sequencing Files>
    More information can be found at https://geoffscollins.github.io/PolTools/tps_distance_per_gene.html

    :param args:
    :return:
    """
    regions_filename, sequencing_files_list, region_length, max_threads = parse_args(args)
    run_tps_distance_per_gene(regions_filename, sequencing_files_list, region_length, max_threads)


if __name__ == '__main__':
    main(sys.argv[1:])
