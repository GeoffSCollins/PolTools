import sys
import multiprocessing

from utils.make_transcripts_dict import build_transcripts_dict
from utils.get_region_length import determine_region_length
from utils.print_tab_delimited import print_tab_delimited
from utils.verify_region_length_is_even import verify_region_length_is_even


def get_pausing_distances_helper(transcripts_dict, regions_filename):
    all_pause_distances = [0] * 101

    # Go through each gene and get the distances from each one
    with open(regions_filename) as file:
        for line in file:
            chromosome, left, right, gene_name, fold_change, strand = line.rstrip().split()

            five_prime_position = int(left) + int(region_length / 2)

            if strand == "-":
                # Subtract one if the strand is negative because the +1 is on the left side
                five_prime_position -= 1

            # Get all of the transcript lengths at this position
            if five_prime_position in transcripts_dict[chromosome][strand]:
                # Go through all the transcripts that have this 5' end
                for three_prime_end in transcripts_dict[chromosome][strand][five_prime_position]:
                    transcript_length = abs(five_prime_position - three_prime_end)

                    # if transcript_length <= 100 and transcript_length > 17:
                    if transcript_length <= 100:
                        all_pause_distances[transcript_length] += 1

    return all_pause_distances


def get_pausing_distances(regions_filename, sequencing_filename):
    transcripts_dict = build_transcripts_dict(sequencing_filename)

    all_pause_distances = get_pausing_distances_helper(transcripts_dict, regions_filename)
    return all_pause_distances


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("python3 pausing_distance_distribution_from_maxTSS.py <regions file> <sequencing files>\n")
    sys.stderr.write("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/pausing_distance_distribution_from_maxTSS.rst\n")


def parse_input(args):
    if len(args) < 2:
        print_usage()
        sys.exit(1)

    regions_filename = args[0]
    global region_length

    region_length = determine_region_length(regions_filename)
    verify_region_length_is_even(region_length)

    sequencing_files = args[1:]

    return regions_filename, region_length, sequencing_files


def output_pausing_distances(pausing_distances, sequencing_files):
    # Print the headers first
    print_tab_delimited(["Transcript Length"] + [seq_filename.split("/")[-1] for seq_filename in sequencing_files])

    for i in range(len(pausing_distances[0])):
        print_tab_delimited([i] + [x[i] for x in pausing_distances])


def pausing_distance_distribution_from_maxTSS(regions_filename, sequencing_files):
    pool = multiprocessing.Pool(processes=len(sequencing_files))
    pausing_distances = pool.starmap(get_pausing_distances, [(regions_filename, seq_filename) for seq_filename in sequencing_files])
    output_pausing_distances(pausing_distances, sequencing_files)


def main(args):
    """
    pausing_distance_distribution_from_maxTSS.py <regions file> <sequencing files>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/pausing_distance_distribution_from_maxTSS.rst

    :param args: list of the sequencing files
    :type args: list
    :return:
    """
    regions_filename, region_length, sequencing_files = parse_input(args)
    pausing_distance_distribution_from_maxTSS(regions_filename, sequencing_files)

if __name__ == '__main__':
    main(sys.argv[1:])
