import sys
import multiprocessing

from utils.make_transcripts_dict import build_transcripts_dict
from utils.get_region_length import determine_region_length
from utils.print_tab_delimited import print_tab_delimited
from utils.verify_region_length_is_even import verify_region_length_is_even


def most_common(lst):
    return max(set(lst), key=lst.count)


def get_pausing_distances_helper(region_filename, transcripts_dict):
    # Go through each gene and get the distances from each one
    ret_dict = {}
    with open(region_filename) as file:
        for line in file:
            chromosome, left, right, gene_name, fold_change, strand = line.rstrip().split()

            five_prime_position = int(left) + region_length / 2

            # Get all of the transcript lengths at this position
            if five_prime_position in transcripts_dict[chromosome][strand]:
                most_common_pausing_position = most_common(transcripts_dict[chromosome][strand][five_prime_position])
                most_common_pausing_distance = abs(five_prime_position - most_common_pausing_position)

            else:
                most_common_pausing_distance = "N/A"

            ret_dict[gene_name] = most_common_pausing_distance

    return ret_dict

def get_pausing_distances(sequencing_filename, region_filename):
    transcripts_dict = build_transcripts_dict(sequencing_filename)
    ret_dict = get_pausing_distances_helper(region_filename, transcripts_dict)

    return ret_dict, sequencing_filename


def print_usage():
    print("Usage: ")
    print("python3 tps_distance_per_gene.py <Regions Filename> <Sequencing Files>")
    print("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/tps_distance_per_gene.rst")



def parse_args(args):
    if len(args) < 2:
        print_usage()
        sys.exit(1)

    regions_filename = args[0]

    global region_length
    region_length = determine_region_length(regions_filename)
    verify_region_length_is_even(region_length)

    sequencing_files_list = args[1:]

    return regions_filename, sequencing_files_list


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


def run_tps_distance_per_gene(regions_filename, sequencing_files_list):
    pool = multiprocessing.Pool(processes=len(sequencing_files_list))
    pausing_distances = pool.starmap(get_pausing_distances, [(sequencing_filename, regions_filename) for sequencing_filename in sequencing_files_list])
    output_data(pausing_distances)


def main(args):
    """
    python3 tps_distance_per_gene.py <Regions Filename> <Sequencing Files>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/tps_distance_per_gene.rst

    :param args:
    :return:
    """
    regions_filename, sequencing_files_list = parse_args(args)
    run_tps_distance_per_gene(regions_filename, sequencing_files_list)


if __name__ == '__main__':
    main(sys.argv[1:])
