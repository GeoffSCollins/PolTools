"""
This program is the interface and driver for tsrFinder
"""

import os
import sys
import argparse
import multiprocessing

from collections import defaultdict
from multiprocessing import Pool

from PolTools.utils.constants import tsr_finder_location
from PolTools.utils.tsr_finder_step_four_from_rocky import run_step_four
from PolTools.utils.remove_files import remove_files

def positive_int(num):
    try:
        val = int(num)
        if val <= 0:
            raise Exception("Go to the except")
    except:
        raise argparse.ArgumentTypeError(num + " must be positive")

    return val

parser = argparse.ArgumentParser(prog='PolTools tsrFinder',
                                 description='Find transcription start regions\n' +
                                 "More information can be found at " +
                                 "https://geoffscollins.github.io/PolTools/tsrFinder.html")

parser.add_argument('seq_file', metavar='seq_file', type=str, help='Bed formatted sequencing file to find the TSRs')
parser.add_argument('window_size', metavar='window_size', type=positive_int, help='Base pair size of the sliding window')
parser.add_argument('min_seq_depth', metavar='min_seq_depth', type=positive_int, help="Minimum number of 5' ends to be considered a TSR")
parser.add_argument('min_avg_transcript_length', metavar='min_avg_transcript_length', type=positive_int, help="Minimum average transcript length to be considered a TSR")
parser.add_argument('max_fragment_size', metavar='max_fragment_size', type=positive_int, help="Maximum fragment size for a read to be counted in tsrFinder")
parser.add_argument('chrom_size_file', metavar='chrom_size_file', type=str, help="Chromosome sizes file")
parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?', default=multiprocessing.cpu_count())

args = parser.parse_args(sys.argv[1:])
bed_file = args.seq_file
window_size = args.window_size
min_seq_depth = args.min_seq_depth
min_avg_transcript_length = args.min_avg_transcript_length
max_fragment_size = args.max_fragment_size
chrom_size_file = args.chrom_size_file
max_threads = args.threads


# Make sure bed_file and chrom_size_file exist
if not os.path.isfile(bed_file):
    sys.stderr.write(bed_file + " was not found. Exiting ...\n")
    sys.exit(1)

if not os.path.isfile(chrom_size_file):
    sys.stderr.write(chrom_size_file + " was not found. Exiting ...\n")
    sys.exit(1)

if not bed_file.endswith(".bed"):
    sys.stderr.write("The sequencing file must end in .bed. Exiting ...\n")
    sys.exit(1)

chromosome_sizes = defaultdict(int)

with open(chrom_size_file) as file:
    for line in file:
        chromosome, size = line.split()
        chromosome_sizes[chromosome] = int(size)

# Step 1. Split the bed file into files by chromosome and strands
fw_filename = bed_file.replace(".bed", "-FW.bed")
rv_filename = bed_file.replace(".bed", "-RV.bed")

parameters_string = "_".join([str(window_size), str(min_seq_depth), str(min_avg_transcript_length), str(max_fragment_size)])

output_filename = bed_file.replace(".bed", "_" + parameters_string + "-TSR.tab")

chromosome_file_writers = defaultdict(lambda : {"+": None, "-": None})

chromosome_files = []
tsr_finder_step_files = []
output_files = []

with open(bed_file) as file:
    for line in file:
        chromosome, left, right, name, score, strand = line.split()

        if chromosome in chromosome_sizes:
            if chromosome not in chromosome_file_writers:
                fw_filename = bed_file.replace(".bed", "-" + chromosome + "-FW.bed")
                rv_filename = bed_file.replace(".bed", "-" + chromosome + "-RV.bed")

                chromosome_file_writers[chromosome]["+"] = open(fw_filename, 'w')
                chromosome_file_writers[chromosome]["-"] = open(rv_filename, 'w')

                chromosome_files.extend([fw_filename, rv_filename])

                for i in range(2, 5):
                    tsr_finder_step_files.append(fw_filename.replace(".bed", "-" + str(i) + "-output.txt"))
                    tsr_finder_step_files.append(rv_filename.replace(".bed", "-" + str(i) + "-output.txt"))

                output_files.append(fw_filename.replace(".bed", "-4-output.txt"))
                output_files.append(rv_filename.replace(".bed", "-4-output.txt"))

            chromosome_file_writers[chromosome][strand].write(line)

# Need to close all the writers
for chromosome in chromosome_file_writers:
    chromosome_file_writers[chromosome]["+"].close()
    chromosome_file_writers[chromosome]["-"].close()


# Step 2: Run tsrFinder on both files concurrently
def run_tsrFinderGC(filename):
    os.system(tsr_finder_location + " " + filename + " " +
              " ".join([str(window_size), str(min_seq_depth), str(min_avg_transcript_length),
                        str(max_fragment_size), chrom_size_file]))

    step_three_filename = filename.replace(".bed", "-3-output.txt")
    step_four_filename = filename.replace(".bed", "-4-output.txt")

    run_step_four(step_three_filename, window_size, chromosome_sizes, step_four_filename)

with Pool(max_threads) as pool:
    pool.map(run_tsrFinderGC, (filename for filename in chromosome_files))

# # Step 3: Combine the output files and delete intermediate files
os.system("cat " + " ".join(output_files) + " > " + output_filename)

remove_files(tsr_finder_step_files)
remove_files(chromosome_files)
