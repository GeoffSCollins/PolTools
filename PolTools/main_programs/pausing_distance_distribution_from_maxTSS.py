import multiprocessing
import sys
import argparse

from PolTools.utils.get_region_length import determine_region_length
from PolTools.utils.make_transcripts_dict import build_transcripts_dict
from PolTools.utils.print_tab_delimited import print_tab_delimited
from PolTools.utils.verify_region_length_is_even import verify_region_length_is_even


def get_pausing_distances_helper(transcripts_dict, regions_filename, max_transcript_length):
    all_pause_distances = [0] * (max_transcript_length + 1)

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
                for three_prime_end, amount in transcripts_dict[chromosome][strand][five_prime_position].items():
                    # The three prime ends are inclusive, so we need to add 1 to the transcript length
                    transcript_length = abs(five_prime_position - three_prime_end) + 1

                    if transcript_length <= max_transcript_length:
                        all_pause_distances[transcript_length] += amount

    return all_pause_distances


def get_pausing_distances(regions_filename, sequencing_filename, max_transcript_length):
    transcripts_dict = build_transcripts_dict(sequencing_filename)
    all_pause_distances = get_pausing_distances_helper(transcripts_dict, regions_filename, max_transcript_length)

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

    parser = argparse.ArgumentParser(prog='PolTools pausing_distance_distribution_from_maxTSS',
                description='Quantify the number of transcripts at each length originating from the max TSS\n' +
                            "More information can be found at " +
                            "https://geoffscollins.github.io/PolTools/pausing_distance_distribution_from_maxTSS.html"
    )

    parser.add_argument('regions_filename', metavar='regions_filename', type=str,
                        help='Bed formatted regions file with an even region length or a region length of one.')

    parser.add_argument('seq_files', metavar='sequencing_files', nargs='+', type=str,
                        help='Bed formatted files from the sequencing experiment')

    parser.add_argument('-m', '--max_transcript_length', dest='max_transcript_length', metavar='max_transcript_length',
                        type=positive_int, nargs='?',
                        help='Set the maximum transcript length to be outputted. Default is 100', default=100)

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    args = parser.parse_args(args)
    regions_filename = args.regions_filename
    sequencing_files = args.seq_files
    max_transcript_length = args.max_transcript_length
    max_threads = args.threads

    region_length = determine_region_length(regions_filename)

    if region_length != 1:
        verify_region_length_is_even(region_length)

    return regions_filename, sequencing_files, max_transcript_length, max_threads


def output_pausing_distances(pausing_distances, sequencing_files):
    # Print the headers first
    print_tab_delimited(["Transcript Length"] + [seq_filename.split("/")[-1] for seq_filename in sequencing_files])

    for i in range(len(pausing_distances[0])):
        print_tab_delimited([i] + [x[i] for x in pausing_distances])


def pausing_distance_distribution_from_maxTSS(regions_filename, sequencing_files, max_transcript_length, max_threads):
    with multiprocessing.Pool(processes=max_threads) as pool:
        pausing_distances = pool.starmap(
            get_pausing_distances,
            [(regions_filename, seq_filename, max_transcript_length) for seq_filename in sequencing_files]
        )

    output_pausing_distances(pausing_distances, sequencing_files)


def main(args):
    """
    pausing_distance_distribution_from_maxTSS <regions file> <sequencing files>
    More information can be found at https://geoffscollins.github.io/PolTools/pausing_distance_distribution_from_maxTSS.html

    :param args: list of the sequencing files
    :type args: list
    :return:
    """
    regions_filename, sequencing_files, max_transcript_length, max_threads = parse_input(args)
    pausing_distance_distribution_from_maxTSS(regions_filename, sequencing_files, max_transcript_length, max_threads)

if __name__ == '__main__':
    main(sys.argv[1:])
