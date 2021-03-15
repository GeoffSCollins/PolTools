import argparse
import sys
import multiprocessing

from GC_bioinfo.utils.verify_bed_file import verify_bed_files
from GC_bioinfo.utils.verify_region_length_is_even import verify_region_length_is_even
from GC_bioinfo.utils.get_region_length import determine_region_length
from GC_bioinfo.utils.heatmap_utils.add_matrices import add_matrices
from GC_bioinfo.utils.remove_files import remove_files

from GC_bioinfo.main_programs.region_heatmap import get_heatmap_matrix, make_heatmap


def get_read_end_combined_heatmap(regions_filename, seq_files_data, end, repeat_amounts, threads):
    with multiprocessing.Pool(threads) as pool:
        args = []
        for seq_file_data in seq_files_data:
            args.append(
                (regions_filename, seq_file_data, end, repeat_amounts)
            )

        matrix_filenames = pool.starmap(get_heatmap_matrix, args)

    # Add the two together
    combined_matrix = add_matrices(matrix_filenames)

    # Remove the files used to combine
    remove_files(matrix_filenames)

    return combined_matrix


def run_read_end_combined_heatmap(args):
    end, filenames, seq_files_data, heatmap_parameters, repeat_amounts, threads, tick_parameters = parse_input(args)

    # Get the matrices for each of the samples
    regions_filename, output_prefix = filenames

    combined_matrix = get_read_end_combined_heatmap(regions_filename, seq_files_data, end, repeat_amounts, threads)

    # Make the heatmap
    make_heatmap(combined_matrix, output_prefix, heatmap_parameters, tick_parameters)


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

    parser = argparse.ArgumentParser(prog=' GC_bioinfo region_fold_change_heatmap')

    parser.add_argument('read_type', metavar='read type', type=str, choices=["five", "three", "whole"],
                        help='either five, three, or whole')

    parser.add_argument('regions_file', metavar='regions_file', type=str,
                        help='Bed formatted file containing all the regions you want to average the sequences')

    parser.add_argument('output_prefix', metavar='output_prefix', type=str, help='Prefix for the output filename')

    parser.add_argument('seq_file_one', metavar='seq_file_one', type=str, help='First sequencing file')
    parser.add_argument('norm_factor_one', metavar='norm_factor_one', type=positive_float, help='Correction factor for the first sequencing file')

    parser.add_argument('seq_file_two', metavar='seq_file_one', type=str, help='Second sequencing file')
    parser.add_argument('norm_factor_two', metavar='norm_factor_two', type=positive_float, help='Correction factor for the second sequencing file')

    parser.add_argument('-m', '--max_black', metavar='max_black', dest='max_black',
                        type=float, default=None,
                        help='Max black value of the heatmap. Default is the maximum possible value')

    parser.add_argument('-r', '--repeat_amount', metavar='repeat_amount', dest='repeat_amount',
                        type=int, default=1,
                        help='Each base will be shown in this number of pixels. Default is 1')

    parser.add_argument('-v', '--vertical_averaging', metavar='vertical_averaging', dest='vertical_averaging',
                        type=int, default=1,
                        help='Average this number of rows into one row. Default is 1')

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    parser.add_argument('-g', '--gamma', metavar='gamma', dest='gamma',
                        type=positive_float, default=2.2,
                        help='Gamma value of the heatmap. Default is 2.2, which is no gamma correction.')

    parser.add_argument('--minor_ticks', metavar='minor_ticks', dest='minor_ticks',
                        type=positive_int, default=None, help='Distance between minor ticks (bp). Default is 10 bp.')

    parser.add_argument('--major_ticks', metavar='major_ticks', dest='major_ticks',
                        type=positive_int, default=None, help='Distance between major ticks (bp). Default is 100 bp')

    args = parser.parse_args(args)

    verify_bed_files(args.regions_file)
    region_length = determine_region_length(args.regions_file)
    verify_region_length_is_even(region_length)

    filenames = args.regions_file, args.output_prefix
    seq_files_data = [(args.seq_file_one, args.norm_factor_one), (args.seq_file_two, args.norm_factor_two)]
    repeat_amounts = args.repeat_amount, args.vertical_averaging
    heatmap_parameters = args.max_value, args.gamma

    if args.minor_ticks != None and args.major_ticks != None:
        tick_parameters = (args.minor_ticks * args.repeat_amount, args.major_ticks * args.repeat_amount)
    else:
        tick_parameters = (None, None)

    return args.read_type, filenames, seq_files_data, heatmap_parameters, repeat_amounts, args.threads, tick_parameters


if __name__ == '__main__':
    run_read_end_combined_heatmap(sys.argv[1:])