"""
The goal of this is to make the plot that David described
"""

import csv
import glob
import os
import sys
from collections import defaultdict

from GC_bioinfo.utils.constants import hg38_chrom_sizes_file
from GC_bioinfo.utils.constants import rna_blacklist_file
from GC_bioinfo.utils.generate_blacklist_regions_for_gene_body_heatmap import blacklist_extended_gene_bodies
from GC_bioinfo.utils.generate_heatmap import generate_heatmap, Ticks
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.make_three_prime_bed_file import make_three_bed_file
from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.run_bedtools_coverage import run_coverage
from GC_bioinfo.utils.run_bedtools_getfasta import run_getfasta
from GC_bioinfo.utils.run_bedtools_subtract import run_subtract
from GC_bioinfo.utils.scale_matrix import scale_matrix


def make_incremented_regions(regions_filename, downstream_distance, interval_size, upstream_distance, min_gene_length):
    # Get the chrom sizes
    chrom_sizes_dict = defaultdict(int)
    with open(hg38_chrom_sizes_file) as file:
        for line in file:
            chromosome, size = line.split()
            chrom_sizes_dict[chromosome] = int(size)

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
                if (int(gene_body_right) - int(gene_body_left)) > min_gene_length:
                    if int(region_left) > 0 and int(region_right) <= chrom_sizes_dict[chromosome]:
                        regions.append([chromosome, region_left, region_right, gene_name, max_tss_five_prime_reads, strand])

    # Go through all the regions and make the incremented ones
    incremented_regions = []
    for region in regions:
        chromosome, left, right, gene_name, score, strand = region

        if strand == "+":
            # We work from left to right
            for i in range(left + interval_size, right + 1, interval_size):
                # Looping through each interval region
                if i-interval_size > 0:
                    incremented_regions.append([chromosome, i - interval_size, i, gene_name, score, strand])

        else:
            # We work from right to left
            for i in range(right, left + (interval_size - 1), (-1 * interval_size)):
                # Looping through each interval region
                if i - interval_size > 0:
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
    lines = [gene_name + "\t" + "\t".join(data[gene_name]) for gene_name in data]
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
                append_string = "\t".join(["0"] * (width - curr_length))

            file.write(line + "\t" + append_string + "\n")

    return sorted_matrix_filename, num_lines


def get_nt_content_dict(nt, truQuant_output_file, content_distance):
    # Need to make a dictionary with the keys of the gene names and values of the nt content for the first distance nts

    # Get the chrom sizes
    chrom_sizes_dict = defaultdict(int)
    with open(hg38_chrom_sizes_file) as file:
        for line in file:
            chromosome, size = line.split()
            chrom_sizes_dict[chromosome] = int(size)

    # Read the genes and make the positions
    positions = defaultdict(list)

    with open(truQuant_output_file) as file:
        for i, line in enumerate(file):
            if i != 0:
                gene_name, chromosome, pause_left, pause_right, strand, total_reads, max_tss, max_tss_five_prime_reads, avg_tss, \
                std_tss, gene_body_left, gene_body_right, *_ = line.split()

                if strand == "+":
                    if int(max_tss) + content_distance <= chrom_sizes_dict[chromosome] and int(max_tss) >= 0:
                        positions[gene_name] = [chromosome, max_tss, str(int(max_tss) + content_distance), gene_name, "0", strand]
                else:
                    # If the strand is negative, we shift the right 1
                    if int(max_tss) - content_distance >= 0 and int(max_tss) + 1 <= chrom_sizes_dict[chromosome]:
                        positions[gene_name] = [chromosome, str(int(max_tss) - content_distance + 1), str(int(max_tss) + 1), gene_name, "0", strand]

    # Now we get the sequences
    sequences = defaultdict(str)

    regions_file = generate_random_filename()

    gene_names = []

    with open(regions_file, 'w') as file:
        for gene_name in positions:
            gene_names.append(gene_name)
            file.write("\t".join(positions[gene_name]) + "\n")

    fasta_file = run_getfasta(regions_file)

    with open(fasta_file) as file:
        for i, line in enumerate(file):
            if i % 2 == 1:
                sequence = line.rstrip().upper()

                curr_gene = gene_names[int(i / 2)]

                sequences[curr_gene] = sequence.count(nt)

    remove_files(regions_file, fasta_file)

    return sequences


def add_nt_content_column(sorted_matrix_filename, nt_content_dict):
    # The last column in this file will be the nt content and it will be sorted descendingly using this column
    with open(sorted_matrix_filename) as file:
        lines = file.readlines()

    output_file = generate_random_filename()

    output_lines = []

    for line in lines:
        gene_name = line.split()[0]

        # Now add the content at the end of the line
        if gene_name in nt_content_dict:
            output_lines.append(line.rstrip() + "\t" + str(nt_content_dict[gene_name]) + "\n")

    output_lines.sort(key=lambda x: float(x.split()[-1]))

    with open(output_file, 'w') as file:
        for line in output_lines:
            file.write(line)

    return output_file


def tmp_fun(mat):
    output_filename = generate_random_filename()
    with open(mat) as file:
        with open(output_filename, 'w') as output_file:
            for line in file:
                output_file.write("\t".join(line.split()[1:-1]) + "\n")

    return output_filename


def build_matrix(seq_file_data, matrix_params, filenames, nt, nt_content_distance, min_gene_length):
    sequencing_filename, spike_in = seq_file_data
    truQuant_output_file, tsr_file, output_filename_prefix = filenames
    upstream_distance, distance_past_tes, width, height, interval_size = matrix_params

    blacklist_regions_file = blacklist_extended_gene_bodies(tsr_file, distance_past_tes)

    # Step 1. Make regions to quantify
    intervals_filename = make_incremented_regions(truQuant_output_file, distance_past_tes, interval_size,
                                                  upstream_distance, min_gene_length)

    # Step 2. Quantify them
    quantified_regions_filename = quantify_intervals(sequencing_filename, blacklist_regions_file,
                                                     intervals_filename)

    # Step 3. Read the coverage data and add them to a 2d list
    sorted_matrix_filename, num_lines = read_coverage_file(quantified_regions_filename, width)

    # Add the final column with the nt data
    nt_content_dict = get_nt_content_dict(nt, truQuant_output_file, nt_content_distance)
    matrix_with_gene_and_nt_content = add_nt_content_column(sorted_matrix_filename, nt_content_dict)

    matrix = tmp_fun(matrix_with_gene_and_nt_content)

    remove_files(intervals_filename, quantified_regions_filename, blacklist_regions_file)

    spike_in_normalized_matrix = scale_matrix(matrix, spike_in)

    remove_files(sorted_matrix_filename, matrix, matrix_with_gene_and_nt_content)

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
    sys.stderr.write(
        "\nMore information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/gene_body_heatmap.rst\n")


def get_args(args):
    if len(args) != 14:
        print_usage()
        sys.exit(1)

    truQuant_output_file, upstream_distance, distance_past_tes, bp_width, width, height, gamma, max_black_value, spike_in, \
    sequencing_filename, output_filename_prefix, nt, nt_content_distance, min_gene_length = args

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

    return seq_file_data, matrix_params, heatmap_params, filenames, nt, nt_content_distance, min_gene_length


def main(args):
    seq_file_data, matrix_params, heatmap_params, filenames, nt, nt_content_distance, min_gene_length = get_args(args)

    output_filename_prefix = filenames[-1]

    matrix = build_matrix(seq_file_data, matrix_params, filenames, nt, nt_content_distance, min_gene_length)

    make_heatmap(matrix, heatmap_params, output_filename_prefix)

    remove_files(matrix)


if __name__ == '__main__':
    main(sys.argv[1:])

