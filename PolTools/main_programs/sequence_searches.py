import sys
import argparse

from itertools import product
import Bio.Data.IUPACData as iupac

from PolTools.utils.get_region_length import determine_region_length
from PolTools.utils.print_tab_delimited import print_tab_delimited
from PolTools.utils.remove_files import remove_files
from PolTools.utils.bedtools_utils.run_bedtools_getfasta import run_getfasta
from PolTools.utils.verify_region_length_is_even import verify_region_length_is_even


class InvalidSearchException(Exception):
    def __init__(self, search):
        super().__init__(search)
        self.search = search

    def __repr__(self):
        return "InvalidSearchException: The following search was not valid: " + self.search


def extend_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    d = iupac.ambiguous_dna_values
    return list(map("".join, product(*map(d.get, seq))))


def find_sequences(gene_dict, search, region_length):
    sequence, left, right = search

    left = int(region_length / 2) + int(left)
    right = int(region_length / 2) + int(right)

    expanded_sequence_list = extend_ambiguous_dna(sequence)

    for seq in expanded_sequence_list:
        if seq in gene_dict["Sequence"][left:right + 1]:
            gene_dict["Motifs"][sequence] = True
            return


def clean_sequence_searches(sequence_searches):
    cleaned_sequence_searches = []

    for search in sequence_searches:
        sequence, region = search.split(",")
        sequence = sequence.upper()
        left, right = region.split("_")

        cleaned_sequence_searches.append([sequence, left, right])

    return cleaned_sequence_searches


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


def parse_input(args):
    parser = argparse.ArgumentParser(prog='PolTools sequence_searches',
                                     description='Determine if a sequence in present in a region around the max TSS\n' +
                                     'More information can be found at ' +
                                     'https://geoffscollins.github.io/PolTools/sequence_searches.html')

    parser.add_argument('regions_filename', metavar='regions_filename', type=str,
                        help='Bed formatted regions file with an even region length or a region length of one.')

    parser.add_argument('search', metavar='search', type=str, nargs='+',
                        help='Search region formatted as follows: (Sequence),(-/+)left:(-/+)right. Ex: TATA,-30:-20')

    args = parser.parse_args(args)

    regions_file = args.regions_filename
    searching_sequences = args.search

    region_length = determine_region_length(regions_file)

    verify_region_length_is_even(region_length)

    # Max position can be either positive or negative
    max_position = int(region_length / 2)

    cleaned_search_sequences = []
    for search in searching_sequences:
        try:
            sequence, search_region = search.split(",")
            left, right = search_region.split(":")

            left = int(left)
            right = int(right)

        except Exception as _:
            raise InvalidSearchException(search)

        if right < left:
            raise InvalidSearchException(search)

        # Check in the left and right positions are in the area
        if left < (max_position * -1):
            raise InvalidSearchException(search)
        if  right > max_position:
            raise InvalidSearchException(search)

        cleaned_search_sequences.append([sequence, left, right])

    return regions_file, cleaned_search_sequences, region_length


def run_sequence_searches(regions_file, searching_sequences, region_length):
    # Has keys of gene names and values of dictionary which has keys of "Sequence" and "Region" and "Motifs".
    master_dict = {}

    # 1. Read in the contents of the bed file
    with open(regions_file) as file:
        bed_lines = file.readlines()

    fasta_file = run_getfasta(regions_file)
    fasta_sequences = read_fasta(fasta_file)
    remove_files(fasta_file)

    # Fill the dictionary
    for i, line in enumerate(bed_lines):
        chromosome, left, right, gene_name, _, strand = line.split()

        master_dict[gene_name] = {"Sequence": "", "Region": line.split(), "Motifs": {}}
        master_dict[gene_name]["Sequence"] = fasta_sequences[i]

        for search in searching_sequences:
            sequence, _, _ = search

            master_dict[gene_name]["Motifs"][sequence] = False

    for gene_name in master_dict:
        for search in searching_sequences:
            find_sequences(master_dict[gene_name], search, region_length)

    # Output the results
    print_tab_delimited(["Chromosome", "Left", "Right", "Gene", "Score", "Strand"] +
                        [search[0] for search in searching_sequences])

    for gene_name in master_dict:
        print_tab_delimited(master_dict[gene_name]["Region"] +
                            [has_motif for _, has_motif in master_dict[gene_name]["Motifs"].items()])


def main(args):
    """
    PolTools sequence_searches <regions file> <Sequence,startPosition:endPosition>
    Ex. PolTools sequence_searches genes.bed TATA,-34:-28
    More information can be found at https://geoffscollins.github.io/PolTools/sequence_searches.html

    :param args: list of the sequencing files
    :type args: list
    :return:
    """
    regions_file, searching_sequences, region_length = parse_input(args)
    run_sequence_searches(regions_file, searching_sequences, region_length)


if __name__ == '__main__':
    main(sys.argv[1:])
