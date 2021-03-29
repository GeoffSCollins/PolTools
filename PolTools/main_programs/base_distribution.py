import sys
import argparse

from PolTools.utils.check_dependencies import check_dependencies
from PolTools.utils.get_region_length import determine_region_length
from PolTools.utils.print_tab_delimited import print_tab_delimited
from PolTools.utils.remove_files import remove_files
from PolTools.utils.bedtools_utils.run_bedtools_getfasta import run_getfasta
from PolTools.utils.verify_bed_file import verify_bed_files
from PolTools.utils.verify_region_length_is_even import verify_region_length_is_even


def read_fasta(filename):
    """
    Returns the sequences from a fasta file.

    :param filename: filename of the fasta file
    :type filename: str
    :return: list
    """
    sequences = []
    with open(filename) as file:
        for i, line in enumerate(file):
            if i % 2 == 0:
                pass
            else:
                sequences.append(line.rstrip())

    return sequences

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


def output_data(avgs_dict, region_length):
    # Write the header
    print_tab_delimited(["Position"] + list(avgs_dict.keys()))

    position = region_length / 2 * -1

    # We go through each position and output the averages
    for i in range(region_length):
        if position == 0:
            position += 1

        print_tab_delimited([position] + [avgs_dict[nt][i] for nt in avgs_dict.keys()])

        position += 1


def run_base_distribution(regions_file, region_length):
    # 1. Get the sequences of the region
    fasta_file = run_getfasta(regions_file)
    sequences = read_fasta(fasta_file)
    remove_files(fasta_file)

    # 2. Get the percentages at each position
    avgs_dict = calculate_averages(sequences)

    # 3. Output into a file
    output_data(avgs_dict, region_length)


def parse_args(args):
    parser = argparse.ArgumentParser(prog='PolTools base_distribution',
        description='Compute the average base composition at each position of the given region.\n' +
                    "More information can be found at " +
                    "https://geoffscollins.github.io/PolTools/base_distribution.html")

    parser.add_argument('regions_file', metavar='regions_file', type=str,
                        help='Bed formatted file containing all the regions you want to average the sequences')

    args = parser.parse_args(args)

    regions_file = args.regions_file

    verify_bed_files(regions_file)
    region_length = determine_region_length(regions_file)
    verify_region_length_is_even(region_length)
    return args.regions_file, region_length


def main(args):
    """
    base_distribution <Regions Filename>
    More information can be found at https://geoffscollins.github.io/PolTools/base_distribution.html

    :param args: list of the sequencing files
    :type args: list
    :return:
    """
    check_dependencies("bedtools")
    regions_file, region_length = parse_args(args)
    run_base_distribution(regions_file, region_length)


if __name__ == '__main__':
    main(sys.argv[1:])
