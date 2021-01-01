import math
import multiprocessing
import sys

from GC_bioinfo.utils.get_region_length import determine_region_length
from GC_bioinfo.utils.make_five_and_three_dict import build_counts_dict
from GC_bioinfo.utils.print_tab_delimited import print_tab_delimited
from GC_bioinfo.utils.verify_bed_file import verify_bed_files
from GC_bioinfo.utils.verify_region_length_is_even import verify_region_length_is_even


def get_counts_helper(five_prime_counts_dict, regions_filename):
    # Counts_dict has keys of gene_names and values of the number of inr_counts
    counts_dict = {}

    # Loop through each region
    with open(regions_filename) as file:
        for line in file:
            chromosome, left, right, gene_name, _, strand = line.split()

            region_length = int(right) - int(left)
            inr_position = int(left) + int(region_length / 2)

            if strand == "-" and region_length != 1:
                # Subtract one if the strand is negative because the +1 is on the left side
                inr_position -= 1

            # Check if the there is a read at that position
            if chromosome in five_prime_counts_dict:
                inr_counts = five_prime_counts_dict[chromosome][strand][inr_position]
            else:
                inr_counts = 0

            counts_dict[gene_name] = inr_counts

    return counts_dict


def get_counts(regions_filename, sequencing_filename):
    five_prime_counts_dict, _ = build_counts_dict(sequencing_filename)
    counts_dict = get_counts_helper(five_prime_counts_dict, regions_filename)

    return counts_dict


def output_data(output_list, sequencing_filenames):

    # First print out the header
    print_tab_delimited(["Gene"] + sequencing_filenames)

    # We will loop through every gene using the first dictionary in the list
    for gene_name in output_list[0]:
        # Print the counts for each of the datasets
        print_tab_delimited([gene_name] + [output_dict[gene_name] for output_dict in output_list])


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo inr_reads.py <regions file> <sequencing files>\n")
    sys.stderr.write("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/inr_reads.rst\n")


def parse_args(args):
    if len(args) < 2:
        print_usage()
        sys.exit(1)

    regions_filename = args[0]
    sequencing_files_list = args[1:]

    verify_bed_files(regions_filename, sequencing_files_list)

    region_length = determine_region_length(regions_filename)

    if region_length != 1:
        verify_region_length_is_even(region_length)

    return regions_filename, sequencing_files_list


def run_inr_reads(regions_filename, sequencing_files_list):
    pool = multiprocessing.Pool(processes=len(sequencing_files_list))
    output_list = pool.starmap(get_counts, [(regions_filename, seq_file) for seq_file in sequencing_files_list])

    output_data(output_list, [seq_file.split("/")[-1] for seq_file in sequencing_files_list])


def main(args):
    """
    GC_bioinfo inr_reads.py <regions file> <sequencing files>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/inr_reads.rst

    :param args:
    :return:
    """
    regions_filename, sequencing_files_list = parse_args(args)
    run_inr_reads(regions_filename, sequencing_files_list)


if __name__ == '__main__':
    main(sys.argv[1:])
