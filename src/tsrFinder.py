"""
This program is the interface and driver for tsrFinder
"""

import os
import sys
from multiprocessing import Process
from collections import defaultdict

from utils.tsr_finder_step_four_from_rocky import run_step_four


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo tsrFinder <Sequencing File> <TSR Window Size> <TSR Read Depth> <Minimum Average Read Length> ")
    sys.stderr.write("<Maximum Fragment Length> <Chromosome Sizes File>\n")
    sys.stderr.write("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/tsrFinder.rst\n")


# Step 0: Setup. Get all the command line arguments. Build the chromosome sizes dictionary.
try:
    args = sys.argv[1:]

    if len(args) != 6:
        print_usage()
        sys.exit(1)

    bed_file = args[0]
    window_size = int(args[1])
    min_seq_depth = int(args[2])
    min_avg_transcript_length = int(args[3])
    max_fragment_size = int(args[4])
    chrom_size_file = args[5]

    # Make sure bed_file and chrom_size_file exist
    if not os.path.isfile(bed_file) or not os.path.isfile(chrom_size_file):
        print_usage()
        sys.exit(1)


except:
    # If the user did not input an argument correctly
    print_usage()
    sys.exit(1)


chromosome_sizes = defaultdict(int)

with open(chrom_size_file) as file:
    for line in file:
        chromosome, size = line.split()

        chromosome_sizes[chromosome] = int(size)



# Step 1. Split the bed file into files by chromosome and strands

fw_filename = bed_file.replace(".bed", "-FW.bed")
rv_filename = bed_file.replace(".bed", "-RV.bed")

output_filename = bed_file.replace(".bed", "-TSR.txt")

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
    os.system("utils/tsrFinder " + filename + " " + " ".join(args[1:]))

    step_three_filename = filename.replace(".bed", "-3-output.txt")
    step_four_filename = filename.replace(".bed", "-4-output.txt")

    run_step_four(step_three_filename, window_size, chromosome_sizes, step_four_filename)



jobs = []
for filename in chromosome_files:
    j = Process(target=run_tsrFinderGC, args=(filename, ))
    j.start()
    jobs.append(j)

for job in jobs:
    job.join()



# # Step 3: Combine the output files and delete intermediate files
os.system("cat " + " ".join(output_files) + " > " + output_filename)
os.system("rm " + " ".join(tsr_finder_step_files + chromosome_files))