import sys
import argparse
import multiprocessing

from GC_bioinfo.utils.get_region_length import determine_region_length
from GC_bioinfo.utils.verify_bed_file import verify_bed_files
from GC_bioinfo.utils.verify_region_length_is_even import verify_region_length_is_even
from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.heatmap_utils.make_log_two_fold_change_matrix import make_log_two_fold_change_matrix
from GC_bioinfo.utils.nested_multiprocessing_pool import NestedPool
from GC_bioinfo.main_programs.read_end_combined_heatmap import get_read_end_combined_heatmap
from GC_bioinfo.utils.heatmap_utils.generate_heatmap import generate_heatmap


def make_heatmap(log2_matrix, max_log2_fold_change, output_filename):
    # Define gamma as 2.2 even though it won't be used
    gamma = 2.2

    if max_log2_fold_change:
        min_value = -1 * max_log2_fold_change
    else:
        min_value = None

    max_value = max_log2_fold_change

    generate_heatmap(log2_matrix, 'red/blue', output_filename, gamma, min_value, max_value)


def run_read_end_fold_change_heatmap(args):
    end, filenames, max_log2_fc, repeat_amounts, numerator_seq_files_data, denominator_seq_files_data, threads = parse_input(args)

    regions_file, output_prefix = filenames

    if threads == 1:
        pool_threads = 1
        comb_heatmap_threads = 1
    elif threads == 2 or threads == 3:
        pool_threads = 2
        comb_heatmap_threads = 1
    else:
        pool_threads = 2
        comb_heatmap_threads = 2

    # Get the numerators and denominators matrix
    with NestedPool(pool_threads) as pool:
        args = [
            (regions_file, numerator_seq_files_data, end, repeat_amounts, comb_heatmap_threads),
            (regions_file, denominator_seq_files_data, end, repeat_amounts, comb_heatmap_threads)
        ]

        numerator_matrix, denominator_matrix = pool.starmap(get_read_end_combined_heatmap, args)

    # Do the log2 fold change for them
    log2_matrix = make_log_two_fold_change_matrix(numerator_matrix, denominator_matrix)
    remove_files(numerator_matrix, denominator_matrix)

    px_per_bp, vertical_averaging = repeat_amounts
    output_filename = output_prefix.replace(".tiff", "") + "_max_" + str(max_log2_fc) + \
                      "_vertical_averaging_" + str(vertical_averaging) +  "_px_per_bp_" + str(px_per_bp) + ".tiff"

    # Make the heatmap
    make_heatmap(log2_matrix, max_log2_fc, output_filename)
    remove_files(log2_matrix)


def parse_input(args):
    # TODO:

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

    parser = argparse.ArgumentParser(prog='end_fold_change_heatmap')

    parser.add_argument('read_type', metavar='read type', type=str, choices=["five", "three", "whole"],
                        help='either five, three, or whole')

    parser.add_argument('regions_file', metavar='regions_file', type=str,
                        help='Bed formatted file containing all the regions you want to average the sequences')

    parser.add_argument('output_prefix', metavar='output_prefix', type=str, help='Prefix for the output filename')

    parser.add_argument('numerator_filename_one', metavar='numerator_filename_one', type=str, help='First sequencing file in the numerator')
    parser.add_argument('numerator_norm_factor_one', metavar='numerator_norm_factor_one', type=positive_float,
                        help='Correction factor for the first numerator sequencing file')
    parser.add_argument('numerator_filename_two', metavar='numerator_filename_two', type=str, help='Second sequencing file in the numerator')
    parser.add_argument('numerator_norm_factor_two', metavar='numerator_norm_factor_two', type=positive_float,
                        help='Correction factor for the second numerator sequencing file')

    parser.add_argument('denominator_filename_one', metavar='denominator_filename_one', type=str, help='First sequencing file in the denominator')
    parser.add_argument('denominator_norm_factor_one', metavar='denominator_norm_factor_one', type=positive_float,
                        help='Correction factor for the first denominator sequencing file')
    parser.add_argument('denominator_filename_two', metavar='denominator_filename_two', type=str, help='Second sequencing file in the denominator')
    parser.add_argument('denominator_norm_factor_two', metavar='denominator_norm_factor_two', type=positive_float,
                        help='Correction factor for the second denominator sequencing file')

    parser.add_argument('-m', '--max_log2_fc', metavar='max_log2_fc', dest='max_log2_fc',
                        type=positive_float, default=None, help='Max log2 fold change of the heatmap')

    parser.add_argument('-r', '--repeat_amount', metavar='repeat_amount', dest='repeat_amount',
                        type=int, default=1,
                        help='Each base will be shown in this number of pixels. Default is 1')

    parser.add_argument('-v', '--vertical_averaging', metavar='vertical_averaging', dest='vertical_averaging',
                        type=int, default=1,
                        help='Average this number of rows into one row. Default is 1')

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    args = parser.parse_args(args)

    numerator_seq_files_data = ((args.numerator_filename_one, args.numerator_norm_factor_one), (args.numerator_filename_two, args.numerator_norm_factor_two))
    denominator_seq_files_data = ((args.denominator_filename_one, args.denominator_norm_factor_one), (args.denominator_filename_two, args.denominator_norm_factor_two))

    verify_bed_files(args.regions_file)
    region_length = determine_region_length(args.regions_file)
    verify_region_length_is_even(region_length)

    filenames = args.regions_file, args.output_prefix
    repeat_amounts = args.repeat_amount, args.vertical_averaging

    return args.read_type, filenames, args.max_log2_fc, repeat_amounts, numerator_seq_files_data, denominator_seq_files_data, args.threads


if __name__ == '__main__':
    run_read_end_fold_change_heatmap(sys.argv[1:])
