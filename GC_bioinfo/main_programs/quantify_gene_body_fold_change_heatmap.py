"""
Get the 2d list (which is not a matrix) and then sum at each position to get the combined then do fold difference
then get average of each position
"""

import sys
import glob
import statistics
import multiprocessing
import math

from collections import defaultdict

from main_programs import gene_body_heatmap, gene_body_fold_change_heatmap

from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.nested_multiprocessing_pool import NestedPool


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo quantify_gene_body_fold_change_heatmap <truQuant output file> <5' Buffer Distance>" +
                     " <Distance Past TES> <Width> <Height> <Gamma> <Max fold change> <Interval Size> <Spike in Correction> <Sequencing Filename>" +
                     "<Numerator Spike in Correction> <Numerator Sequencing Filename> <Numerator Spike in Correction> <Numerator Sequencing Filename>" +
                     "<Denomenator Spike in Correction> <Denomenator Sequencing Filename> <Denomenator Spike in Correction> <Denomenator Sequencing Filename> <Output Filename> \n")
    sys.stderr.write("\nMore information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/quantify_gene_body_fold_change_heatmap.rst\n")


def read_coverage_file(coverage_file):
    # Goal is to make a table with the gene name and then all of the values
    read_coverage_data = defaultdict(list)

    with open(coverage_file) as file:
        for line in file:
            chrom, left, right, gene_name, score, strand, counts, _, _, _ = line.split()
            read_coverage_data[gene_name].append(int(counts))

    return read_coverage_data


def normalize(coverage_data, spike_in):
    for gene_name in coverage_data:
        for i in range(len(coverage_data[gene_name])):
            coverage_data[gene_name][i] *= spike_in

    return coverage_data


def get_coverage_data(sequencing_filename, spike_in, blacklist_regions_file, intervals_filename):
    # Step 2. Quantify them
    quantified_regions_filename = gene_body_heatmap.quantify_intervals(sequencing_filename, blacklist_regions_file, intervals_filename)

    # Step 3. Read the coverage data and add them to a 2d list
    coverage_data_dict = read_coverage_file(quantified_regions_filename)

    remove_files(quantified_regions_filename)

    normalized_coverage_data = normalize(coverage_data_dict, spike_in)

    return normalized_coverage_data


def add_coverage_data(coverage_dict_one, coverage_dict_two):
    # Given two dicts just add them together
    # Will use coverage_dict_one as the final dict
    for gene_name in coverage_dict_one:
        for i in range(len(coverage_dict_one[gene_name])):
            coverage_dict_one[gene_name][i] += coverage_dict_two[gene_name][i]

    return coverage_dict_one


def get_combined_coverage_dict(filename_one, spike_in_one, filename_two, spike_in_two, blacklist_regions_file, intervals_filename):
    pool = multiprocessing.Pool(processes=2)
    all_coverage_data = pool.starmap(get_coverage_data, ((filename_one, spike_in_one, blacklist_regions_file, intervals_filename),
                                                         (filename_two, spike_in_two, blacklist_regions_file, intervals_filename)))
    return add_coverage_data(all_coverage_data[0], all_coverage_data[1])


def get_fold_change_dict(numerator_dict, denominator_dict):
    # Use the numerator dict as the one to return
    for gene_name in numerator_dict:
        for i in range(len(numerator_dict[gene_name])):
            numerator = numerator_dict[gene_name][i]
            denominator = denominator_dict[gene_name][i]

            if numerator == 0:
                numerator = 1

            if denominator == 0:
                denominator = 1

            fold_change = numerator / denominator

            # Get the log 2 fold change
            log_two_fold_change = math.log2(fold_change)

            numerator_dict[gene_name][i] = log_two_fold_change

    return numerator_dict


def get_averages_list(fold_change_dict):
    # This will go through each position in the lists and average them
    # This will return a list
    averages_list = []

    longest_gene_name = max(fold_change_dict, key=lambda x:len(fold_change_dict[x]))

    longest_gene_length = len(fold_change_dict[longest_gene_name])

    for i in range(longest_gene_length):
        all_values = []
        for gene_name in fold_change_dict:
            try:
                all_values.append(fold_change_dict[gene_name][i])
            except IndexError:
                # Just go to next gene
                continue

        # Once we are here we can append it to the averages list
        averages_list.append(statistics.mean(all_values))

    return averages_list


def main(args):
    gene_body_file, five_prime_buffer_distance, distance_past_tes, width, height, gamma, max_fold_change, \
    interval_size, numerator_spike_in_one, numerator_sequencing_filename_one, numerator_spike_in_two, numerator_sequencing_filename_two, \
    denominator_spike_in_one, denominator_sequencing_filename_one, denominator_spike_in_two, denominator_sequencing_filename_two, output_filename_prefix = gene_body_fold_change_heatmap.get_args(args)

    # Find all regions to blacklist
    tsr_file = glob.glob(gene_body_file.replace("gene_body_regions.bed", "") + "*TSR.tab")

    if not tsr_file:
        sys.stderr.write("No tsrFinder file was found. Exiting ...\n")
        sys.exit(1)

    if len(tsr_file) != 1:
        sys.stderr.write("More than one tsrFinder file was found for this run of truQuant. Exiting ...\n")
        sys.exit(1)

    tsr_file = tsr_file[0]

    blacklist_regions_file = gene_body_heatmap.blacklist_extended_gene_bodies(tsr_file, distance_past_tes)

    # Step 1. Make regions to quantify
    intervals_filename = gene_body_heatmap.make_incremented_regions(gene_body_file, 0, interval_size,
                                                                    five_prime_buffer_distance)

    # Now get all the coverage data
    # Use the nested multiprocessing from gene_body_fold_change_heatmap
    pool = NestedPool(2)

    numerator_args = (numerator_sequencing_filename_one, numerator_spike_in_one, numerator_sequencing_filename_two, numerator_spike_in_two, blacklist_regions_file, intervals_filename)
    denominator_args = (denominator_sequencing_filename_one, denominator_spike_in_one, denominator_sequencing_filename_two, denominator_spike_in_two, blacklist_regions_file, intervals_filename)

    numerator_coverage_dict, denominator_coverage_dict = pool.starmap(get_combined_coverage_dict, [numerator_args, denominator_args])

    fold_change_dict = get_fold_change_dict(numerator_coverage_dict, denominator_coverage_dict)

    # Now we need to get the average at each position
    averages_list = get_averages_list(fold_change_dict)

    # Print the position and the averages_list value
    start_location = five_prime_buffer_distance

    output_filename = output_filename_prefix + "_max_" + str(max_fold_change) + "_width_" + str(
        width * interval_size) + "bp_quantification.txt"

    with open(output_filename, 'w') as file:
        file.write("\t".join(["Location", "Average Fold Change"]) + "\n")

        for val in averages_list:
            file.write("\t".join([str(start_location), str(val)]) + "\n")
            start_location += interval_size

    remove_files(intervals_filename, blacklist_regions_file)


if __name__ == '__main__':
    main(sys.argv[1:])