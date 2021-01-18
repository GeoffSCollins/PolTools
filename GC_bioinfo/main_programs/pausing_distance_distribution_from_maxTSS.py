import multiprocessing
import sys
import argparse

from GC_bioinfo.utils.get_region_length import determine_region_length
from GC_bioinfo.utils.make_transcripts_dict import build_transcripts_dict
from GC_bioinfo.utils.print_tab_delimited import print_tab_delimited
from GC_bioinfo.utils.verify_region_length_is_even import verify_region_length_is_even


def get_pausing_distances_helper(transcripts_dict, regions_filename):
    all_pause_distances = [0] * 101

    # Go through each gene and get the distances from each one
    with open(regions_filename) as file:
        for line in file:
            chromosome, left, right, gene_name, fold_change, strand = line.rstrip().split()

            region_length = int(right) - int(left)
            five_prime_position = int(left) + int(region_length / 2)

            if strand == "-" and region_length != 1:
                # Subtract one if the strand is negative because the +1 is on the left side
                five_prime_position -= 1

            # Get all of the transcript lengths at this position
            if five_prime_position in transcripts_dict[chromosome][strand]:
                # Go through all the transcripts that have this 5' end
                for three_prime_end in transcripts_dict[chromosome][strand][five_prime_position]:
                    # The three prime ends are inclusive, so we need to add 1 to the transcript length
                    transcript_length = abs(five_prime_position - three_prime_end) + 1

                    if transcript_length <= 100:
                        all_pause_distances[transcript_length] += 1

    return all_pause_distances


def get_pausing_distances(regions_filename, sequencing_filename):
    transcripts_dict = build_transcripts_dict(sequencing_filename)
    all_pause_distances = get_pausing_distances_helper(transcripts_dict, regions_filename)

    return all_pause_distances


def parse_input(args):
    def positive_int(num):
        try:
            val = int(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")

        return val

    parser = argparse.ArgumentParser(prog='GC_bioinfo pausing_distance_distribution_from_maxTSS',
                description='Quantify the number of transcripts at each length originating from the max TSS\n' +
                            "More information can be found at " +
                            "https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/pausing_distance_distribution_from_maxTSS.rst"
    )

    parser.add_argument('regions_filename', metavar='regions_filename', type=str,
                        help='Bed formatted regions file with an even region length or a region length of one.')

    parser.add_argument('seq_files', metavar='sequencing_files', nargs='+', type=str,
                        help='Bed formatted files from the sequencing experiment')

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    args = parser.parse_args(args)
    regions_filename = args.regions_filename
    sequencing_files = args.seq_files
    max_threads = args.threads

    region_length = determine_region_length(regions_filename)

    if region_length != 1:
        verify_region_length_is_even(region_length)

    return regions_filename, sequencing_files, max_threads


def output_pausing_distances(pausing_distances, sequencing_files):
    # Print the headers first
    print_tab_delimited(["Transcript Length"] + [seq_filename.split("/")[-1] for seq_filename in sequencing_files])

    for i in range(len(pausing_distances[0])):
        print_tab_delimited([i] + [x[i] for x in pausing_distances])


def pausing_distance_distribution_from_maxTSS(regions_filename, sequencing_files, max_threads):
    pool = multiprocessing.Pool(processes=max_threads)
    pausing_distances = pool.starmap(get_pausing_distances, [(regions_filename, seq_filename) for seq_filename in sequencing_files])
    pool.close()
    output_pausing_distances(pausing_distances, sequencing_files)


def main(args):
    """
    pausing_distance_distribution_from_maxTSS <regions file> <sequencing files>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/pausing_distance_distribution_from_maxTSS.rst

    :param args: list of the sequencing files
    :type args: list
    :return:
    """
    regions_filename, sequencing_files, max_threads = parse_input(args)
    pausing_distance_distribution_from_maxTSS(regions_filename, sequencing_files, max_threads)

if __name__ == '__main__':
    main(sys.argv[1:])
