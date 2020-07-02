import sys
import os
import csv

from pathlib import Path

from utils.fasta_reader import read_fasta
from utils.make_random_filename import generate_random_filename
from utils.verify_bed_file import verify_bed_files


def calculate_averages(sequences):
    a_sums = [0] * len(sequences[0])
    t_sums = [0] * len(sequences[0])
    g_sums = [0] * len(sequences[0])
    c_sums = [0] * len(sequences[0])

    for sequence in sequences:
        for i, char in enumerate(sequence):
            if char.upper() == "A":
                a_sums[i] += 1
            elif char.upper() == "T":
                t_sums[i] += 1
            elif char.upper() == "G":
                g_sums[i] += 1
            elif char.upper() == "C":
                c_sums[i] += 1

    a_avgs = [0] * len(sequences[0])
    g_avgs = [0] * len(sequences[0])
    c_avgs = [0] * len(sequences[0])
    t_avgs = [0] * len(sequences[0])

    # Loop through all of the values in the sums and compute the average
    for i in range(len(a_sums)):
        # The averages are equal to the sum / the number of sequences
        a_avgs[i] = a_sums[i] / len(sequences)
        g_avgs[i] = g_sums[i] / len(sequences)
        c_avgs[i] = c_sums[i] / len(sequences)
        t_avgs[i] = t_sums[i] / len(sequences)

    return a_avgs, t_avgs, g_avgs, c_avgs


def output_to_file(filename, a_avgs, t_avgs, g_avgs, c_avgs):
    with open(filename, 'w') as output_file:
        tsv_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')
        # Write a header
        tsv_writer.writerow(["Position", "A", "T", "G", "C"])

        position = len(a_avgs) / 2 * -1

        # We go through each position and output the averages
        for i in range(len(a_avgs)):
            if position == 0:
                position += 1

            tsv_writer.writerow([position, a_avgs[i], t_avgs[i], g_avgs[i], c_avgs[i]])

            position += 1


def print_usage():
    print("Usage: ")
    print("python3 base_distribution.py <Regions Filename>")
    print("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/base_distribution.rst")

def parse_input():
    if len(sys.argv) == 1 or len(sys.argv) > 2:
        print_usage()
        sys.exit(1)

    regions_file = sys.argv[1]
    verify_bed_files(regions_file)

    return regions_file


if __name__ == '__main__':
    regions_file = parse_input()

    if regions_file == "":
        print("No region file was given. Exiting...")
        exit(1)

    _file_path = str(Path(__file__).parent.absolute())
    hg38_fasta_file = _file_path + "/static/hg38.fa"

    # 1. Make a fasta file
    random_filename = generate_random_filename()
    os.system("bedtools getfasta -s -fi " + hg38_fasta_file + " -bed " + regions_file + " > " + random_filename)
    sequences = read_fasta(random_filename)
    os.system("rm " + random_filename)

    # 2. Get the percentages at each position
    a_avgs, t_avgs, g_avgs, c_avgs = calculate_averages(sequences)

    output_filename = regions_file.replace(".bed", "-base_distribution_plot.tsv")

    # 3. Output into a file
    output_to_file(output_filename, a_avgs, t_avgs, g_avgs, c_avgs)
