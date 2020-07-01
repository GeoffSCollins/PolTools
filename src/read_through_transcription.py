'''
We are going to look at the TES and go to 10 kb past that and get all of the 3' ends
We will output for every 50 bp
'''

import csv
import sys

from utils.run_bedtools_coverage import run_coverage
from utils.run_bedtools_subtract import run_subtract
from utils.make_three_prime_bed_file import make_three_bed_file
from utils.make_random_filename import generate_random_filename
from utils.remove_files import remove_files

# The user must give us a bed file with the regions and a list of the sequencing files
regions_file, sequencing_files, tsr_file, output_filename = sys.argv[1:]

sequencing_files = sequencing_files.split(",")

# 1. Make the region from -1 kb to +100 kb
regions = []

upstream_distance = 1_000
downstream_distance = 100_000


with open(regions_file) as file:
    for line in file:
        chromosome, left, right, gene_name, score, strand = line.split()

        if strand == "+":
            # We use the right position because that is the TES
            region_left = int(right) - upstream_distance
            region_right = int(right) + downstream_distance
        else:
            # We use the left position because that is the TES
            region_left = int(left) - downstream_distance
            region_right = int(left) + upstream_distance

        # Add the region to regions
        regions.append([chromosome, region_left, region_right, gene_name, score, strand])

# 2. Split up the regions into 50 bp increments and order them from in the gene to past the TES
incremented_regions = []
for region in regions:
    chromosome, left, right, gene_name, score, strand = region

    if strand == "+":
        # We work from left to right
        for i in range(left + 50, right + 1, 50):
            # Looping through each 50 bp region
            incremented_regions.append([chromosome, i - 50, i, gene_name, score, strand])

    else:
        # We work from right to left
        if left + 49 > 0:
            for i in range(right, left + 49, -50):
                # Looping through each 50 bp region
                incremented_regions.append([chromosome, i - 50, i, gene_name, score, strand])

# With all the regions defined, we can get all of the counts for this region and output to a temp file
region_intervals_file = generate_random_filename()

with open(region_intervals_file, 'w') as tmp_region_file:
    output_writer = csv.writer(tmp_region_file, delimiter='\t', lineterminator='\n')
    for region in incremented_regions:
        output_writer.writerow(region)

# Go through all of the sequencing files and blacklist TSRs
blacklisted_filenames = []
for seq_filename in sequencing_files:
    curr_filename = generate_random_filename()
    blacklisted_filenames.append(curr_filename)
    run_subtract(seq_filename, tsr_file, output_filename=curr_filename)

# Go through all of the blacklisted sequencing files and make 3' end files
three_end_files = []
for filename in blacklisted_filenames:
    three_end_files.append(make_three_bed_file(filename))

# Now run bedtools coverage on each 3' end file
coverage_files = []
for i, seq_filename in enumerate(three_end_files):
    curr_filename = generate_random_filename()
    coverage_files.append(curr_filename)
    run_coverage(region_intervals_file, three_end_files[i], output_filename=curr_filename)


# Has keys of the region number and values of a list containing the counts for each file
combined_dict = {}

for i, filename in enumerate(coverage_files):
    with open(filename) as file:
        # We have open the bedtools coverage results file
        region_number = 0
        previous_gene_name = ""

        for j, line in enumerate(file):
            chromosome, left, right, gene_name, score, strand, counts, *_ = line.split()

            if previous_gene_name != gene_name:
                region_number = 0

            if region_number not in combined_dict:
                combined_dict[region_number] = [0] * len(sequencing_files)

            combined_dict[region_number][i] += int(counts)

            region_number += 1
            previous_gene_name = gene_name


# Now the data is in the combined_dict, we need to reduce it back down to positions again
with open(output_filename, 'w') as output_file:
    output_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

    # First write the headers
    output_writer.writerow(["Position"] + [seq_file.split("/")[-1] for seq_file in sequencing_files])

    real_position = upstream_distance * -1
    # Now output all of the data
    for position in combined_dict:
        output_writer.writerow([real_position] + combined_dict[position])
        real_position += 50

# Remove all of the temporary files
remove_files(region_intervals_file, blacklisted_filenames, three_end_files, coverage_files)
