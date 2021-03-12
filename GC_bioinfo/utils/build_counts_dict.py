import os

from collections import defaultdict

from GC_bioinfo.utils.constants import hg38_chrom_sizes_random_file
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files

def build_counts_dict(sequencing_filename, read_type):
    """
    Builds two dictionaries for the 5' ends and 3' ends of transcripts. Dictionaries can be accessed like so:
    five_prime_counts_dict[chromosome][strand][five_prime_position], where five_prime_position is an integer

    :param sequencing_filename: filename of the sequencing data collected
    :type sequencing_filename: str
    :return: two dictionaries, five_prime_counts_dict and three_prime_counts_dict
    :rtype: (dict, dict)
    """

    counts_dict = {}

    if read_type not in ["five", "three", "whole"]:
        raise ValueError("Read type must be either five, three, or whole")

    if read_type == 'whole':
        # Use bedtools genome coverage because it is much faster
        fw_bedgraph = generate_random_filename(".bedGraph")
        rv_bedgraph = generate_random_filename(".bedGraph")

        os.system("bedtools genomecov -i " + sequencing_filename + " -g " + hg38_chrom_sizes_random_file + \
                  " -bg -strand + > " + fw_bedgraph)

        os.system("bedtools genomecov -i " + sequencing_filename + " -g " + hg38_chrom_sizes_random_file + \
                  " -bg -strand - > " + rv_bedgraph)

        with open(fw_bedgraph) as file:
            for line in file:
                chromosome, left, right, counts = line.split()

                left, right = int(left), int(right)

                if chromosome not in counts_dict:
                    counts_dict[chromosome] = {
                        "+": defaultdict(int),
                        "-": defaultdict(int)
                    }

                for position in range(left, right):
                    counts_dict[chromosome]["+"][position] += 1

        with open(rv_bedgraph) as file:
            for line in file:
                chromosome, left, right, counts = line.split()

                left, right = int(left), int(right)

                for position in range(left, right):
                    counts_dict[chromosome]["-"][position] += 1

        remove_files(fw_bedgraph, rv_bedgraph)

        return counts_dict


    with open(sequencing_filename) as file:
        for i, line in enumerate(file):
            chromosome, left, right, _, _, strand = line.rstrip().split()

            left = int(left)
            right = int(right)

            if chromosome not in counts_dict:
                counts_dict[chromosome] = {
                    "+": defaultdict(int),
                    "-": defaultdict(int)
                }

            if read_type == "five":
                if strand == "+":
                    position = left
                else:
                    position = right - 1

            elif read_type == "three":
                if strand == "+":
                    position = right - 1
                else:
                    position = left

            counts_dict[chromosome][strand][position] += 1

    return counts_dict