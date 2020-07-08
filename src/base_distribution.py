import sys
import csv

from utils.fasta_reader import read_fasta
from utils.verify_bed_file import verify_bed_files
from utils.check_dependencies import check_dependencies
from utils.run_bedtools_getfasta import run_getfasta
from utils.remove_files import remove_files
from utils.get_region_length import determine_region_length
from utils.print_tab_delimited import print_tab_delimited


def calculate_averages(sequences):
    sums_dict = {
        "A": [0] * len(sequences[0]),
        "T": [0] * len(sequences[0]),
        "G": [0] * len(sequences[0]),
        "C": [0] * len(sequences[0])
    }

    for sequence in sequences:
        for i, char in enumerate(sequence):
            if char.upper() in sums_dict:
                sums_dict[char.upper()][i] += 1

    averages_dict = {
        "A": [0] * len(sequences[0]),
        "T": [0] * len(sequences[0]),
        "G": [0] * len(sequences[0]),
        "C": [0] * len(sequences[0])
    }

    # We loop through each position in the list
    # Loop through each base in the dictionary
    # Do the division

    for i in range(len(sequences[0])):
        # Looping through each position in the region
        for nt in averages_dict:
            # Looping through nucleotides (A, T, G, C)
            averages_dict[nt][i] = sums_dict[nt][i] / len(sequences)

    return averages_dict


def output_data(avgs_dict):
    # Write the header
    print_tab_delimited(["Position"] + list(avgs_dict.keys()))

    position = region_length / 2 * -1

    # We go through each position and output the averages
    for i in range(region_length):
        if position == 0:
            position += 1

        print_tab_delimited([position] + [avgs_dict[nt][i] for nt in avgs_dict.keys()])

        position += 1


def print_usage():
    print("Usage: ")
    print("python3 base_distribution.py <Regions Filename>")
    print("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/base_distribution.rst")


def parse_input(args):
    if len(args) != 1:
        print_usage()
        sys.exit(1)

    regions_file = args[0]
    verify_bed_files(regions_file)

    global region_length
    region_length = determine_region_length(regions_file)

    if region_length % 2 != 0:
        print("The region length is not even, so the +1 nt position cannot be determined. Exiting ...")
        sys.exit(1)

    return regions_file


def run_base_distribution(regions_file):
    # 1. Get the sequences of the region
    fasta_file = run_getfasta(regions_file)
    sequences = read_fasta(fasta_file)
    remove_files(fasta_file)

    # 2. Get the percentages at each position
    avgs_dict = calculate_averages(sequences)

    # 3. Output into a file
    output_data(avgs_dict)


def main(args):
    """
    base_distribution.py <Regions Filename>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/base_distribution.rst

    :param args: list of the sequencing files
    :type args: list
    :return:
    """
    check_dependencies("bedtools")
    regions_file = parse_input(args)
    run_base_distribution(regions_file)


if __name__ == '__main__':
    main(sys.argv[1:])
