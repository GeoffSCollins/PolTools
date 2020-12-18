"""
Do ±20 bp from the TES and plot every 200 bp the number of 3' ends
"""

import csv
import glob
import os
import sys
from collections import defaultdict

from GC_bioinfo.utils.average_matrix import average_matrix
from GC_bioinfo.utils.constants import rna_blacklist_file
from GC_bioinfo.utils.generate_blacklist_regions_for_gene_body_heatmap import blacklist_extended_gene_bodies
from GC_bioinfo.utils.generate_heatmap import generate_heatmap, Ticks
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.make_three_prime_bed_file import make_three_bed_file
from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.run_bedtools_coverage import run_coverage
from GC_bioinfo.utils.run_bedtools_subtract import run_subtract
from GC_bioinfo.utils.scale_matrix import scale_matrix


def make_incremented_regions(truQuant_output_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, interval_size):
    # Using the regions provided, make incremented regions centered on the TES
    with open(truQuant_output_file) as file:
        lines = []
        for i, line in enumerate(file):
            if i != 0:
                lines.append(line)

    # Sort the lines by gene length
    lines.sort(key=lambda x: int(x.split()[12]))

    regions = []
    for line in lines:
        gene_name, chromosome, pause_left, pause_right, strand, total_reads, max_tss, max_tss_five_prime_reads, avg_tss, \
        std_tss, gene_body_left, gene_body_right, gene_body_length, *_ = line.split()

        if strand == "+":
            region_left = int(max_tss) - upstream_distance
            region_right = int(gene_body_right) + downstream_distance

            # Cut off the region if it goes before the start of the gene body - distance_upstream_of_gene_body_start
            if region_left < (int(max_tss) - distance_upstream_of_gene_body_start):
                region_left = int(max_tss) - distance_upstream_of_gene_body_start
        else:
            region_left = int(gene_body_left) - downstream_distance
            region_right = int(max_tss) + upstream_distance

            # Cut off the region if it goes before the start of the gene body + distance_upstream_of_gene_body_start
            if region_right > int(max_tss) + distance_upstream_of_gene_body_start:
                region_right = int(max_tss) + distance_upstream_of_gene_body_start


        if region_left < 0:
            # If the region is negative, this is not possible so just don't add it the regions
            continue

        # Add the region to regions
        regions.append([chromosome, region_left, region_right, gene_name, max_tss_five_prime_reads, strand])

    # Go through all the regions and make the incremented ones
    incremented_regions = []
    for region in regions:
        chromosome, left, right, gene_name, score, strand = region

        if strand == "+":
            # We work from left to right
            for i in range(left + interval_size, right + 1, interval_size):
                # Looping through each interval region
                incremented_regions.append([chromosome, i - interval_size, i, gene_name, score, strand])

        else:
            # We work from right to left
            for i in range(right, left + (interval_size - 1), (-1 * interval_size)):
                # Looping through each interval region
                incremented_regions.append([chromosome, i - interval_size, i, gene_name, score, strand])

    region_intervals_filename = generate_random_filename()

    with open(region_intervals_filename, 'w') as tmp_region_file:
        output_writer = csv.writer(tmp_region_file, delimiter='\t', lineterminator='\n')
        for region in incremented_regions:
            output_writer.writerow(region)

    return region_intervals_filename


def quantify_intervals(sequencing_filename, blacklist_filename, intervals_file):
    # First make the sequencing file 3' ends only
    three_prime_filename = make_three_bed_file(sequencing_filename)

    # Now blacklist the sequencing file
    blacklisted_sequencing_filename = run_subtract(three_prime_filename, blacklist_filename, rna_blacklist_file)

    # Use bedtools coverage to get the number of 3' reads for each interval
    coverage_file = run_coverage(intervals_file, blacklisted_sequencing_filename)

    # We can now remove the 3' reads file
    remove_files(blacklisted_sequencing_filename, three_prime_filename)

    return coverage_file


def read_coverage_file(coverage_file, width):
    # Goal is to make a table with the gene name and then all of the values
    data = defaultdict(list)

    with open(coverage_file) as file:
        for line in file:
            chrom, left, right, gene_name, score, strand, counts, _, _, _ = line.split()
            data[gene_name].append(counts)


    # Write the dictionary to a file by sorting by gene length
    lines = ["\t".join(data[gene_name]) for gene_name in data]
    num_lines = len(lines)

    sorted_matrix_filename = generate_random_filename()

    with open(sorted_matrix_filename, 'w') as file:
        for line in lines:
            # Need to add 0's to get to the same length
            curr_length = len(line.split())

            if curr_length >= width:
                line = "\t".join(line.split()[(-1*width):])
                append_string = ""
            else:
                append_string = "\t".join(["0"] * (width - curr_length))

            file.write(append_string + "\t" + line + "\n")

    return sorted_matrix_filename, num_lines


def get_matrix(tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start,
               truQuant_output_file, interval_size, sequencing_filename, width, height, spike_in):

    blacklist_regions_file = blacklist_extended_gene_bodies(tsr_file, downstream_distance)

    # Make the intervals file to quantify
    intervals_filename = make_incremented_regions(truQuant_output_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, interval_size)

    # Quantify it
    quantified_regions_filename = quantify_intervals(sequencing_filename, blacklist_regions_file, intervals_filename)

    # Make a matrix for the data
    sorted_matrix_filename, num_lines = read_coverage_file(quantified_regions_filename, width)
    num_lines_to_average = int(num_lines / height)

    # Average the matrix
    averaged_matrix_filename = average_matrix(sorted_matrix_filename, num_lines_to_average)

    # Spike in normalize
    spike_in_matrix = scale_matrix(averaged_matrix_filename, spike_in)

    remove_files(blacklist_regions_file, averaged_matrix_filename, sorted_matrix_filename, quantified_regions_filename, intervals_filename)

    return spike_in_matrix


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo TES_heatmap <truQuant output file> <Width> <Height> <Downstream Distance> <Upstream Distance> <Distance Upstream of Gene Body Start> "
                     "<Gamma> <Max black value> <Spike in Correction> <Sequencing Filename> <Output Filename> \n")
    sys.stderr.write("\nMore information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/gene_body_heatmap.rst\n")


def get_args(args):
    if len(args) != 11:
        print_usage()
        sys.exit(1)

    truQuant_output_file, width, height, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, \
        gamma, max_black_value, spike_in, sequencing_filename, output_prefix = args

    gene_body_file = glob.glob(truQuant_output_file.replace("-truQuant_output.txt", "") + "*gene_body_regions.bed")

    if not gene_body_file:
        sys.stderr.write("No gene body file was found for that run of truQuant. Exiting ...\n")
        sys.exit(1)

    if len(gene_body_file) != 1:
        sys.stderr.write("More than one gene body file was found. Exiting ...\n")
        sys.exit(1)

    gene_body_file = gene_body_file[0]

    # Find all regions to blacklist
    tsr_file = glob.glob(gene_body_file.replace("gene_body_regions.bed", "") + "*TSR.tab")

    if not tsr_file:
        sys.stderr.write("No tsrFinder file was found. Exiting ...\n")
        sys.exit(1)

    if len(tsr_file) != 1:
        sys.stderr.write("More than one tsrFinder file was found for this run of truQuant. Exiting ...\n")
        sys.exit(1)

    tsr_file = tsr_file[0]

    if not os.path.isfile(sequencing_filename):
        sys.stderr.write("File " + sequencing_filename + " was not found.\n")
        sys.exit(1)

    def try_to_convert_to_int(var, var_name):
        try:
            int_var = int(var)
            return int_var
        except ValueError:
            sys.stderr.write("The " + var_name + " could not be converted to an integer")
            sys.exit(1)


    width = try_to_convert_to_int(width, "width")
    height = try_to_convert_to_int(height, "height")
    downstream_distance = try_to_convert_to_int(downstream_distance, "downstream distance")
    upstream_distance = try_to_convert_to_int(upstream_distance, "upstream distance")
    distance_upstream_of_gene_body_start = try_to_convert_to_int(distance_upstream_of_gene_body_start, "distance upstream of gene body start")
    interval_size = int((downstream_distance + upstream_distance) / width)

    # If the interval size is not an integer, then we can't use it
    if (downstream_distance + upstream_distance) % width:
        sys.stderr.write("The heatmap width in px must be a factor of the base pair width (bp width / px width must be an integer)")
        sys.exit(1)

    try:
        gamma = float(gamma)
    except ValueError:
        sys.stderr.write("The gamma could not be converted to a float")

    try:
        max_black_value = float(max_black_value)
    except ValueError:
        sys.stderr.write("The gamma could not be converted to a float")

    try:
        spike_in = float(spike_in)
    except ValueError:
        sys.stderr.write("The spike in correction factor could not be converted to a float")

    return truQuant_output_file, tsr_file, width, height, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, \
           gamma, max_black_value, interval_size, spike_in, sequencing_filename, output_prefix


def main(args):
    truQuant_output_file, tsr_file, width, height, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, \
        gamma, max_black_value, interval_size, spike_in, sequencing_filename, output_prefix = get_args(args)

    spike_in_matrix = get_matrix(tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, truQuant_output_file, interval_size, sequencing_filename, width, height, spike_in)

    # Make the heatmap
    output_filename = output_prefix + "_max_" + str(max_black_value) + "_gamma_" + str(gamma) + "_width_" + str(downstream_distance + upstream_distance) + "bp_TES_heatmap.tiff"

    t = Ticks(10_000 / interval_size, 50_000 / interval_size)
    generate_heatmap(spike_in_matrix, "gray", output_filename, gamma, 0, max_black_value, t)

    remove_files(spike_in_matrix)


if __name__ == '__main__':
    main(sys.argv[1:])