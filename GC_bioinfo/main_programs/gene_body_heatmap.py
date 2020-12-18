"""
The goal of this is to make the plot that David described
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


def make_incremented_regions(regions_filename, downstream_distance, interval_size, upstream_distance):
    # Using the regions provided, make incremented regions
    with open(regions_filename) as file:
        regions = []
        for i, line in enumerate(file):
            if i != 0:
                gene_name, chromosome, pause_left, pause_right, strand, total_reads, max_tss, max_tss_five_prime_reads, avg_tss, \
                                    std_tss, gene_body_left, gene_body_right, *_ = line.split()

                if strand == "+":
                    region_left = int(max_tss) - upstream_distance
                    region_right = int(gene_body_right) + downstream_distance
                else:
                    region_left = int(gene_body_left) - downstream_distance
                    region_right = int(max_tss) + upstream_distance

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
            if left + (interval_size - 1) > 0:
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
    remove_files(three_prime_filename)

    # Use bedtools coverage to get the number of 3' reads for each interval
    coverage_file = run_coverage(intervals_file, blacklisted_sequencing_filename)

    # We can now remove the 3' reads file
    remove_files(blacklisted_sequencing_filename)

    return coverage_file


def read_coverage_file(coverage_file, width):
    # Goal is to make a table with the gene name and then all of the values
    data = defaultdict(list)

    with open(coverage_file) as file:
        for line in file:
            chrom, left, right, gene_name, score, strand, counts, _, _, _ = line.split()
            data[gene_name].append(counts)


    # Write the dictionary to a file
    lines = ["\t".join(data[gene_name]) for gene_name in data]
    lines.sort(key=lambda x: len(x.split()))
    num_lines = len(lines)

    sorted_matrix_filename = generate_random_filename()

    with open(sorted_matrix_filename, 'w') as file:
        for line in lines:
            # Need to add 0's to get to the same length
            curr_length = len(line.split())

            if curr_length >= width:
                line = "\t".join(line.split()[:width])
                append_string = ""
            else:
                append_string = "\t".join(["0"]*(width - curr_length))

            file.write(line + "\t" + append_string + "\n")

    return sorted_matrix_filename, num_lines


def build_matrix(seq_file_data, matrix_params, filenames):
    sequencing_filename, spike_in = seq_file_data
    truQuant_output_file, tsr_file, output_filename_prefix = filenames
    upstream_distance, distance_past_tes, width, height, interval_size = matrix_params

    blacklist_regions_file = blacklist_extended_gene_bodies(tsr_file, distance_past_tes)

    # Step 1. Make regions to quantify
    intervals_filename = make_incremented_regions(truQuant_output_file, distance_past_tes, interval_size,
                                                  upstream_distance)

    # Step 2. Quantify them
    quantified_regions_filename = quantify_intervals(sequencing_filename, blacklist_regions_file, intervals_filename)

    # Step 3. Read the coverage data and add them to a 2d list
    sorted_matrix_filename, num_lines = read_coverage_file(quantified_regions_filename, width)

    remove_files(intervals_filename, quantified_regions_filename, blacklist_regions_file)

    # Make the dimensions correct
    num_lines_to_average = int(num_lines / height)
    averaged_matrix_filename = average_matrix(sorted_matrix_filename, num_lines_to_average)

    remove_files(sorted_matrix_filename)

    spike_in_normalized_matrix = scale_matrix(averaged_matrix_filename, spike_in)

    remove_files(averaged_matrix_filename)

    return spike_in_normalized_matrix


def make_heatmap(matrix, heatmap_params, output_filename_prefix):
    bp_width, width, height, gamma, max_black_value, interval_size = heatmap_params

    # Minor tick marks every 10 kb and major tick marks every 50 kb
    t = Ticks(minor_tick_mark_interval_size=(10_000 / interval_size),
              major_tick_mark_interval_size=(50_000 / interval_size))

    output_filename = output_filename_prefix + "_gamma_" + str(gamma) + "_max_" + str(
        max_black_value) + "_width_" + str(bp_width) + "bp_gene_body_heatmap.tiff"

    generate_heatmap(matrix, 'gray', output_filename, gamma, 0, max_black_value, ticks=t)


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo gene_body_heatmap <truQuant output file> <Upstream Distance>" +
                     " <Distance Past TES> <Width (bp)> <Width (px)> <Height> <Gamma> <Max black value> <Spike in Correction> <Sequencing Filename> <Output Filename> \n")
    sys.stderr.write("\nMore information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/gene_body_heatmap.rst\n")


def get_args(args):
    if len(args) != 11:
        print_usage()
        sys.exit(1)

    truQuant_output_file, upstream_distance, distance_past_tes, bp_width, width, height, gamma, max_black_value, spike_in,\
        sequencing_filename, output_filename_prefix = args

    tsr_file = glob.glob(truQuant_output_file.replace("-truQuant_output.txt", "") + "*TSR.tab")

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

    # Make sure the distance_past_tes, width, max_gene_length are all integers
    upstream_distance = try_to_convert_to_int(upstream_distance, "5' buffer distance")
    distance_past_tes = try_to_convert_to_int(distance_past_tes, "distance past the TES")
    width = try_to_convert_to_int(width, "width (px)")
    bp_width = try_to_convert_to_int(bp_width, "width (bp)")
    height = try_to_convert_to_int(height, "height")

    interval_size = bp_width / width
    interval_size = try_to_convert_to_int(interval_size, "interval size")

    if bp_width % width != 0:
        sys.stderr.write("The width (bp) must be evenly divisible by the width (px). Exiting ...")
        sys.exit(1)

    if bp_width < width:
        sys.stderr.write("The width (bp) must be greater than width (px). Exiting ...")
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

    seq_file_data = (sequencing_filename, spike_in)
    matrix_params = (upstream_distance, distance_past_tes, width, height, interval_size)
    heatmap_params = (bp_width, width, height, gamma, max_black_value, interval_size)
    filenames = (truQuant_output_file, tsr_file, output_filename_prefix)

    return seq_file_data, matrix_params, heatmap_params, filenames


def main(args):
    seq_file_data, matrix_params, heatmap_params, filenames = get_args(args)

    output_filename_prefix = filenames[-1]

    matrix = build_matrix(seq_file_data, matrix_params, filenames)

    make_heatmap(matrix, heatmap_params, output_filename_prefix)

    remove_files(matrix)

if __name__ == '__main__':
    main(sys.argv[1:])
