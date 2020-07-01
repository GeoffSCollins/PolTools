'''
Takes in a regions file in bed format and a fasta file and outputs the bed file and plus sequence columns True / False on the right

We are searching for TATA (-25 to -35), BRE (CGCC -30 to -40), DPE (RGWYV +25 to +35), MTE (CSARSSAAC +15 to +30)
'''

# Todo: Test this program

import sys
import csv
from Bio import Seq
from itertools import product

from utils import run_bedtools_getfasta, fasta_reader, get_region_length


def extend_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    d = Seq.IUPAC.IUPACData.ambiguous_dna_values
    return [list(map("".join, product(*map(d.get, seq))))]


def find_sequences(gene_dict, search, region_length):
    sequence, left, right = search

    left = (region_length / 2) + int(left)
    right = (region_length / 2) + int(right)

    expanded_sequence_list = extend_ambiguous_dna(sequence)

    for seq in expanded_sequence_list:
        if seq in gene_dict["Sequence"][left:right+1]:
            gene_dict["Motifs"][sequence] = True
            return

    master_dict[gene_name]["Motifs"][sequence] = False


def clean_sequence_searches(sequence_searches):
    cleaned_sequence_searches = []

    for search in sequence_searches:
        sequence, region = search.split(",")
        sequence = sequence.upper()
        left, right = region.split("_")

        cleaned_sequence_searches.append([sequence, left, right])

    return cleaned_sequence_searches


if __name__ == '__main__':
    regions_file = sys.argv[1]
    searching_sequences = list(sys.argv[2:])

    # How to call the program
    # python3 sequence_searches.py <regions file> <Sequence,startPosition_endPosition>
    # Ex. python3 sequence_searches.py genes.bed TATA,-34_-28

    cleaned_sequence_searches = clean_sequence_searches(searching_sequences)
    region_length = get_region_length.determine_region_length(regions_file)

    output_filename = regions_file.replace(".bed", "sequence_search.tsv")

    # Has keys of gene names and values of dictionary which has keys of "Sequence" and "Region" and "Motifs".
    master_dict = {}

    # 1. Read in the contents of the bed file
    with open(regions_file) as file:
        bed_lines = file.readlines()

    fasta_file = run_bedtools_getfasta.run_getfasta(regions_file)
    fasta_sequences = fasta_reader.read_fasta(fasta_file)

    # Fill the dictionary
    for i, line in enumerate(bed_lines):
        chromosome, left, right, gene_name, _, strand = line.split()

        master_dict[gene_name] = {"Sequence": "", "Region": line.split(), "Motifs": {}}

        master_dict[gene_name]["Sequence"] = fasta_sequences

        for search in cleaned_sequence_searches:
            sequence, _, _ = search

            master_dict[gene_name]["Motifs"][sequence] = False

    for gene_name in master_dict:
        for search in cleaned_sequence_searches:
            find_sequences(master_dict[gene_name], search, region_length)

    # Output the results
    with open(output_filename, 'w') as output_file:
        output_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')
        output_writer.writerow(
            ["Chromosome", "Left", "Right", "Gene", "Fold Change", "Strand", "TATA", "TATAWAAR", "BRE", "DPE", "MTE",
             "INR"])
        for gene_name in master_dict:
            to_output = master_dict[gene_name]["Region"] + [master_dict[gene_name]["Motifs"]["TATA"]] + [
                master_dict[gene_name]["Motifs"]["TATAWAAR"]] + \
                        [master_dict[gene_name]["Motifs"]["BRE"]] + [master_dict[gene_name]["Motifs"]["DPE"]] + \
                        [master_dict[gene_name]["Motifs"]["MTE"]] + [master_dict[gene_name]["Motifs"]["INR"]]

            output_writer.writerow(to_output)
