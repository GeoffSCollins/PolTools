import sys
import argparse
import os
import multiprocessing
import glob
import math
import csv
from collections import defaultdict
from GC_bioinfo.utils.nested_multiprocessing_pool import NestedPool
from GC_bioinfo.utils.make_random_filename import generate_random_filename

from GC_bioinfo.main_programs.gene_body_heatmap import blacklist_extended_gene_bodies, quantify_intervals

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

    parser = argparse.ArgumentParser(prog='GC_bioinfo gene_body_fold_change_heatmap',
                                     description="Generate a heatmap of 3' ends for each gene sorted by gene length " +
                                                 "aligned by the TSS\n" +
                                                 "More information can be found at " +
                                                 "https://geoffscollins.github.io/GC_bioinfo/gene_body_fold_change_heatmap.html")

    parser.add_argument('truQuant_output_file', metavar='truQuant_output_file', type=str,
                        help='truQuant output file which ends in -truQuant_output.txt')

    parser.add_argument('numerator_correction_factor_one', metavar='numerator_correction_factor_one', type=positive_float,
                        help='Correction factor for the first numerator dataset')

    parser.add_argument('numerator_seq_file_one', metavar='numerator_seq_file_one', type=str,
                        help='First numerator bed formatted sequencing file')

    parser.add_argument('numerator_correction_factor_two', metavar='numerator_correction_factor_two', type=positive_float,
                        help='Correction factor for the second numerator dataset')

    parser.add_argument('numerator_seq_file_two', metavar='numerator_seq_file_two', type=str,
                        help='Second numerator bed formatted sequencing file')

    parser.add_argument('denominator_correction_factor_one', metavar='denominator_correction_factor_one',
                        type=positive_float,
                        help='Correction factor for the first denominator dataset')

    parser.add_argument('denominator_seq_file_one', metavar='denominator_seq_file_one', type=str,
                        help='First denominator bed formatted sequencing file')

    parser.add_argument('denominator_correction_factor_two', metavar='denominator_correction_factor_two',
                        type=positive_float,
                        help='Correction factor for the second denominator dataset')

    parser.add_argument('denominator_seq_file_two', metavar='denominator_seq_file_two', type=str,
                        help='Second denominator bed formatted sequencing file')

    parser.add_argument('-u', '--upstream_distance', metavar='upstream_distance', dest='upstream_distance',
                        type=positive_int, default=50_000, help='Distance upstream of the max TSS')

    parser.add_argument('-d', '--distance_past_tes', metavar='distance_past_tes', dest='distance_past_tes',
                        type=positive_int, default=50_000, help='Distance downstream of the transcription end site')

    parser.add_argument('-b', '--bp_width', metavar='bp_width', dest='bp_width', default=400_000, type=positive_int,
                        help='Total number of base pairs shown on the heatmap. This number must be greater than the ' +
                             'upstream distance + distance past TES.')

    parser.add_argument('-w', '--width', metavar='width', dest='width',
                        type=positive_int, default=2_000, help='Width of the heatmap in pixels')

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    args = parser.parse_args(args)

    truQuant_output_file = args.truQuant_output_file
    upstream_distance = args.upstream_distance
    distance_past_tes = args.distance_past_tes
    bp_width = args.bp_width
    width = args.width

    numerator_spike_in_one = args.numerator_correction_factor_one
    numerator_sequencing_filename_one = args.numerator_seq_file_one
    numerator_spike_in_two = args.numerator_correction_factor_two
    numerator_sequencing_filename_two = args.numerator_seq_file_two

    denominator_spike_in_one = args.denominator_correction_factor_one
    denominator_sequencing_filename_one = args.denominator_seq_file_one
    denominator_spike_in_two = args.denominator_correction_factor_two
    denominator_sequencing_filename_two = args.denominator_seq_file_two

    max_threads = args.threads

    tsr_file = glob.glob(truQuant_output_file.replace("-truQuant_output.txt", "") + "*TSR.tab")

    if not tsr_file:
        sys.stderr.write("No tsrFinder file was found. Exiting ...\n")
        sys.exit(1)

    if len(tsr_file) != 1:
        sys.stderr.write("More than one tsrFinder file was found for this run of truQuant. Exiting ...\n")
        sys.exit(1)

    tsr_file = tsr_file[0]

    for file in [numerator_sequencing_filename_one, numerator_sequencing_filename_two,
                 denominator_sequencing_filename_one, denominator_sequencing_filename_two]:
        if not os.path.isfile(file):
            sys.stderr.write("File " + file + " was not found.\n")
            sys.exit(1)

    if bp_width % width != 0:
        sys.stderr.write("The width (bp) must be evenly divisible by the width (px). Exiting ...")
        sys.exit(1)

    if bp_width < width:
        sys.stderr.write("The width (bp) must be greater than width (px). Exiting ...")
        sys.exit(1)

    interval_size = int(bp_width / width)

    numerator_seq_files_data = [(numerator_sequencing_filename_one, numerator_spike_in_one),
                                (numerator_sequencing_filename_two, numerator_spike_in_two)]

    denominator_seq_files_data = [(denominator_sequencing_filename_one, denominator_spike_in_one),
                                (denominator_sequencing_filename_two, denominator_spike_in_two)]

    matrix_params = (upstream_distance, distance_past_tes, width, interval_size)
    filenames = (truQuant_output_file, tsr_file)

    return numerator_seq_files_data, denominator_seq_files_data, matrix_params, filenames, max_threads


def make_incremented_regions(regions_filename, interval_size, upstream_distance):
    # The difference between this function in the gene body heatmap one is that we don't go downstream because we don't want to quantify that

    # Using the regions provided, make incremented regions
    with open(regions_filename) as file:
        regions = []
        for i, line in enumerate(file):
            if i != 0:
                gene_name, chromosome, pause_left, pause_right, strand, total_reads, max_tss, max_tss_five_prime_reads, avg_tss, \
                                    std_tss, gene_body_left, gene_body_right, *_ = line.split()

                if strand == "+":
                    region_left = int(max_tss) - upstream_distance
                    region_right = int(gene_body_right)
                else:
                    region_left = int(gene_body_left)
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


def get_quant_list(truQuant_output_file, distances, sequencing_filename, norm_factor, blacklist_regions_file):
    distance_past_tes, interval_size, upstream_distance = distances

    # Step 1. Make regions to quantify
    intervals_filename = make_incremented_regions(truQuant_output_file, interval_size, upstream_distance)

    # Step 2. Quantify them
    quantified_regions_filename = quantify_intervals(sequencing_filename, blacklist_regions_file, intervals_filename)

    data = defaultdict(list)

    with open(quantified_regions_filename) as file:
        for line in file:
            chrom, left, right, gene_name, score, strand, counts, _, _, _ = line.split()
            data[gene_name].append(counts)

    # Go through the quantification and add each position up
    # This line will make the quant list equal to the length of the longest list in the data values
    quant_list = [0] * len(data[max(data, key=lambda x: len(data[x]))])

    for gene_name in data:
        for i, counts in enumerate(data[gene_name]):
            quant_list[i] += int(counts) * norm_factor

    return quant_list


def get_both_quant_lists(truQuant_output_file, distances, seq_files, blacklist_regions_file):
    # Get the quant lists and add them together
    pool = multiprocessing.Pool(2)

    # Get the arguments
    args = []
    for seq_file_data in seq_files:
        seq_file, norm_factor = seq_file_data
        args.append(
            (truQuant_output_file, distances, seq_file, norm_factor, blacklist_regions_file)
        )

    quant_lists = pool.starmap(get_quant_list, args)

    pool.close()

    # Add the quant lists together
    summed_quant_list = [x + y for x, y in zip(*quant_lists)]

    return summed_quant_list


def output_data(log_two_fold_change_list, upstream_distance, interval_size):
    # Want to print out the position and then the log2 fold change
    curr_position = -1 * upstream_distance

    for i, val in enumerate(log_two_fold_change_list):
        print(str(curr_position) + '\t' + str(val))
        curr_position += interval_size


def main(args):
    numerator_seq_files_data, denominator_seq_files_data, matrix_params, filenames, max_threads = get_args(args)

    # We want to make a list of the data for the numerators and denominators
    truQuant_output_file, tsr_file = filenames
    upstream_distance, distance_past_tes, width, interval_size = matrix_params
    distances = distance_past_tes, interval_size, upstream_distance

    blacklist_regions_file = blacklist_extended_gene_bodies(tsr_file, distance_past_tes)

    # Get the list for the numerators and denominators
    pool = NestedPool(2)

    args = [
        (truQuant_output_file, distances, numerator_seq_files_data, blacklist_regions_file),
        (truQuant_output_file, distances, denominator_seq_files_data, blacklist_regions_file)
    ]

    numerator_list, denominator_list = pool.starmap(get_both_quant_lists, args)

    pool.close()

    # Make the log2 fold change and then output it!
    log_two_fold_change = []

    for num, den in zip(numerator_list, denominator_list):
        if num == 0:
            num = 1
        if den == 0:
            den = 1
        log_two_fold_change.append(math.log2(num / den))

    output_data(log_two_fold_change, upstream_distance, interval_size)

if __name__ == '__main__':
    main(sys.argv[1:])