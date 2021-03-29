import csv
import glob
import os
import sys
import argparse
import multiprocessing

from collections import defaultdict

from PolTools.utils.heatmap_utils.average_matrix import average_matrix
from PolTools.utils.constants import rna_blacklist_file
from PolTools.utils.generate_blacklist_regions_for_gene_body_heatmap import blacklist_extended_gene_bodies
from PolTools.utils.heatmap_utils.generate_heatmap import generate_heatmap, Ticks
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.remove_files import remove_files
from PolTools.utils.bedtools_utils.run_bedtools_coverage import run_coverage
from PolTools.utils.bedtools_utils.run_bedtools_subtract import run_subtract
from PolTools.utils.heatmap_utils.scale_matrix import scale_matrix
from PolTools.utils.make_read_end_file import make_read_end_file
from PolTools.utils.heatmap_utils.add_matrices import add_matrices


def make_incremented_regions(truQuant_output_file, downstream_distance, upstream_distance, bp_width, interval_size):
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

            # If the region will go past the image, cut it off at the end of the image using the 5' end
            if region_right - region_left > bp_width:
                region_left = region_right - bp_width
        else:
            region_left = int(gene_body_left) - downstream_distance
            region_right = int(max_tss) + upstream_distance

            # If the region will go past the image, cut it off at the end of the image
            if region_right - region_left > bp_width:
                region_right = region_left + bp_width

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
    three_prime_filename = make_read_end_file(sequencing_filename, 'three')

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

    sorted_matrix_filename = generate_random_filename(".sorted.matrix")

    with open(sorted_matrix_filename, 'w') as file:
        for line in lines:
            # Need to add 0's to get to the same length
            curr_length = len(line.split())
            append_string = "\t".join(["0"] * (width - curr_length))

            file.write(append_string + "\t" + line + "\n")

    return sorted_matrix_filename, num_lines


def get_individual_matrix(filenames, dimensions, spike_in):
    sequencing_filename, blacklist_regions_file, intervals_filename = filenames
    width, height = dimensions
    # Quantify it
    quantified_regions_filename = quantify_intervals(sequencing_filename, blacklist_regions_file, intervals_filename)

    # Make a matrix for the data
    sorted_matrix_filename, num_lines = read_coverage_file(quantified_regions_filename, width)
    num_lines_to_average = int(num_lines / height)

    # Average the matrix
    averaged_matrix_filename = average_matrix(sorted_matrix_filename, num_lines_to_average)

    # Spike in normalize
    spike_in_matrix = scale_matrix(averaged_matrix_filename, spike_in)

    remove_files(quantified_regions_filename, sorted_matrix_filename, averaged_matrix_filename)

    return spike_in_matrix


def get_matrix(seq_files_data, matrix_params, filenames, threads):
    upstream_distance, downstream_distance, bp_width, width, height, interval_size= matrix_params
    truQuant_output_file, tsr_file, output_filename_prefix = filenames

    blacklist_regions_file = blacklist_extended_gene_bodies(tsr_file, downstream_distance)

    # Make the intervals file to quantify
    intervals_filename = make_incremented_regions(truQuant_output_file, downstream_distance, upstream_distance, bp_width, interval_size)

    with multiprocessing.Pool(threads) as pool:
        args = []
        dimensions = width, height

        for dataset in seq_files_data:
            seq_filename, spike_in = dataset
            filenames = [seq_filename, blacklist_regions_file, intervals_filename]
            args.append(
                (filenames, dimensions, spike_in)
            )

        individual_matrices = pool.starmap(get_individual_matrix, args)

    combined_matrix = add_matrices(individual_matrices)
    remove_files(blacklist_regions_file, intervals_filename)

    return combined_matrix


def make_heatmap(matrix, heatmap_params, output_filename_prefix):
    bp_width, width, height, gamma, max_black_value, interval_size, minor_ticks_bp, major_ticks_bp = heatmap_params

    # Minor tick marks every 10 kb and major tick marks every 50 kb
    t = Ticks(minor_tick_mark_interval_size=(minor_ticks_bp / interval_size),
              major_tick_mark_interval_size=(major_ticks_bp / interval_size))

    output_filename = output_filename_prefix + "_max_" + str(max_black_value) + "_gamma_" + str(
        gamma) + "_width_" + str(bp_width) + "bp_TES_heatmap.tiff"

    generate_heatmap(matrix, 'gray', output_filename, gamma, 0, max_black_value, ticks=t)


def get_args(args):
    def positive_int(num):
        try:
            val = int(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")

        return val

    def positive_float(num):
        try:
            val = float(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")

        return val

    parser = argparse.ArgumentParser(prog='PolTools TES_heatmap',
                                     description="Generate a heatmap of 3' ends for each gene sorted by gene length " +
                                     "aligned by the transcription end site\n" +
                                     "More information can be found at " +
                                     "https://geoffscollins.github.io/PolTools/TES_heatmap.html")

    parser.add_argument('truQuant_output_file', metavar='truQuant_output_file', type=str,
                        help='truQuant output file which ends in -truQuant_output.txt')

    parser.add_argument('-s', '--seq_file', action='append', nargs=2, metavar=('seq_file', 'spike_in'),
                        required=True, help='Provide the sequencing file with its correction factor. You can supply '
                                            'more than one sequencing file by adding multiple -s arguments.')

    parser.add_argument('output_prefix', metavar='output_prefix', type=str, help='Prefix for the output filename')

    parser.add_argument('-w', '--width', metavar='width', dest='width',
                        type=positive_int, default=2_000, help='Width of the heatmap in pixels. Default is 2,000 px.')

    parser.add_argument('-e', '--height', metavar='height', dest='height',
                        type=positive_int, default=2_000, help='Height of the heatmap in pixels. Default is 2,000 px.')

    parser.add_argument('-d', '--downstream_distance', metavar='downstream_distance', dest='downstream_distance',
                        type=positive_int, default=50_000, help='Distance downstream from the transcription end site. Default is 50,000 bp.')

    parser.add_argument('-u', '--upstream_distance', metavar='upstream_distance', dest='upstream_distance',
                        type=positive_int, default=50_000, help='Distance upstream of the start of the gene body. Default is 50,000 bp.')

    parser.add_argument('-b', '--bp_width', metavar='bp_width', dest='bp_width', default=400_000, type=positive_int,
                        help='Total number of base pairs shown on the heatmap. This number must be greater than the ' +
                             'upstream distance + distance past TES. Default is 400,000 bp.')

    parser.add_argument('-g', '--gamma', metavar='gamma', dest='gamma',
                        type=positive_float, default=2.2, help='Gamma value of the heatmap. Default is 2.2, which is no gamma correction.')

    parser.add_argument('-m', '--max_black', metavar='max_black', dest='max_black',
                        type=positive_float, default=None, help='Max black value of the heatmap. Default is the max value found.')

    parser.add_argument('--minor_ticks', metavar='minor_ticks', dest='minor_ticks',
                        type=positive_int, default=10_000, help='Distance between minor ticks (bp). Default is 10,000 bp between minor tick marks.')

    parser.add_argument('--major_ticks', metavar='major_ticks', dest='major_ticks',
                        type=positive_int, default=50_000, help='Distance between major ticks (bp). Default is 50,000 bp between minor tick marks.')

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    args = parser.parse_args(args)

    truQuant_output_file = args.truQuant_output_file
    output_filename_prefix = args.output_prefix
    width = args.width
    height = args.height
    downstream_distance = args.downstream_distance
    upstream_distance = args.upstream_distance
    bp_width = args.bp_width
    gamma = args.gamma
    max_black_value = args.max_black
    minor_ticks = args.minor_ticks
    major_ticks = args.major_ticks

    # Find all regions to blacklist
    tsr_file = glob.glob(truQuant_output_file.replace("-truQuant_output.txt", "") + "*TSR.tab")

    if not tsr_file:
        sys.stderr.write("No tsrFinder file was found. Exiting ...\n")
        sys.exit(1)

    if len(tsr_file) != 1:
        sys.stderr.write("More than one tsrFinder file was found for this run of truQuant. Exiting ...\n")
        sys.exit(1)

    tsr_file = tsr_file[0]

    seq_files_data = []

    for dataset in args.seq_file:
        seq_file, corr_factor = dataset

        corr_factor = positive_float(corr_factor)
        seq_files_data.append((seq_file, corr_factor))

        if not os.path.isfile(seq_file):
            sys.stderr.write("File " + seq_file + " was not found.\n")
            sys.exit(1)

    # If the interval size is not an integer, then we can't use it
    if bp_width % width:
        sys.stderr.write("The heatmap width in px must be a factor of the base pair width (bp width / px width must be an integer)")
        sys.exit(1)

    interval_size = int(bp_width / width)

    matrix_params = (upstream_distance, downstream_distance, bp_width, width, height, interval_size)
    heatmap_params = (bp_width, width, height, gamma, max_black_value, interval_size, minor_ticks, major_ticks)
    filenames = (truQuant_output_file, tsr_file, output_filename_prefix)

    return seq_files_data, matrix_params, heatmap_params, filenames, args.threads


def main(args):
    seq_files_data, matrix_params, heatmap_params, filenames, threads = get_args(args)

    matrix = get_matrix(seq_files_data, matrix_params, filenames, threads)

    output_filename_prefix = filenames[-1]

    # Make the heatmap
    make_heatmap(matrix, heatmap_params, output_filename_prefix)

    remove_files(matrix)


if __name__ == '__main__':
    main(sys.argv[1:])