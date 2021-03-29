import sys
import argparse
import multiprocessing
import os

from PIL import Image
from PolTools.main_programs import region_heatmap

from PolTools.utils.get_region_length import determine_region_length
from PolTools.utils.verify_bed_file import verify_bed_files
from PolTools.utils.verify_region_length_is_even import verify_region_length_is_even
from PolTools.utils.remove_files import remove_files
from PolTools.utils.heatmap_utils.make_log_two_fold_change_matrix import make_log_two_fold_change_matrix
from PolTools.utils.nested_multiprocessing_pool import NestedPool
from PolTools.utils.heatmap_utils.generate_heatmap import generate_heatmap, Ticks, make_ticks_matrix
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.constants import generate_heatmap_location

def make_ticks_image(width, tick_params, tick_width):
    # Make the tick marks
    t = Ticks(*tick_params, offset=1, width=tick_width)

    # Ticks matrix with a height of 50 px and a max black value of 1
    ticks_matrix = make_ticks_matrix(width, 50, 1, t)

    # Write to a file
    ticks_matrix_filename = generate_random_filename('.matrix')
    with open(ticks_matrix_filename, 'w') as file:
        for row in ticks_matrix:
            file.write("\t".join([str(val) for val in row]) + "\n")

    ticks_image_filename = generate_random_filename(".tiff")

    os.system("/usr/bin/Rscript " + generate_heatmap_location + " " +
              " ".join([ticks_matrix_filename, "gray", ticks_image_filename, "2.2"]))

    remove_files(ticks_matrix_filename)

    return ticks_image_filename


def combine_images(ticks_image_filename, only_heatmap_filename, output_filename):
    ticks_image = Image.open(ticks_image_filename)
    heatmap_image = Image.open(only_heatmap_filename)

    final_image = Image.new('RGB', (ticks_image.width, ticks_image.height + heatmap_image.height))

    final_image.paste(heatmap_image, (0, 0))
    final_image.paste(ticks_image, (0, heatmap_image.height))

    final_image.save(output_filename + ".tiff")
    remove_files(ticks_image_filename, only_heatmap_filename)


def make_heatmap(log2_matrix, max_log2_fold_change, tick_parameters, output_filename, px_per_bp):
    # Define gamma as 2.2 even though it won't be used
    gamma = 2.2
    min_value = -1 * max_log2_fold_change if max_log2_fold_change else None
    max_value = max_log2_fold_change
    red_blue_image = generate_random_filename('.tiff')

    if tick_parameters != (None, None):
        # Get the width of the matrix so we can make the ticks
        with open(log2_matrix) as file:
            width = len(file.readline().split())

        ticks_image = make_ticks_image(width, tick_parameters, px_per_bp)
        generate_heatmap(log2_matrix, 'red/blue', red_blue_image, gamma, min_value, max_value)
        combine_images(ticks_image, red_blue_image, output_filename)

    else:
        generate_heatmap(log2_matrix, 'red/blue', output_filename, gamma, min_value, max_value)


def run_read_end_fold_change_heatmap(args):
    end, filenames, max_log2_fc, repeat_amounts, numerator_seq_files_data, denominator_seq_files_data, \
        threads, tick_parameters = parse_input(args)

    regions_file, output_prefix = filenames

    # We use max_threads / 2 because we will be running two instances of the combined
    threads_per_heatmap = int(threads / 2)

    # Make sure that if the user only wants to run on one thread that it does not default to 0
    if threads_per_heatmap == 0:
        threads_per_heatmap = 1

    # Get the numerators and denominators matrix
    with NestedPool(threads) as pool:
        args = [
            (regions_file, numerator_seq_files_data, end, repeat_amounts, threads_per_heatmap),
            (regions_file, denominator_seq_files_data, end, repeat_amounts, threads_per_heatmap)
        ]

        numerator_matrix, denominator_matrix = pool.starmap(region_heatmap.get_matrix, args)

    # Do the log2 fold change for them
    log2_matrix = make_log_two_fold_change_matrix(numerator_matrix, denominator_matrix)
    remove_files(numerator_matrix, denominator_matrix)

    px_per_bp, vertical_averaging = repeat_amounts
    output_filename = output_prefix.replace(".tiff", "") + "_max_" + str(max_log2_fc) + "_vertical_averaging_" + \
                      str(vertical_averaging) +  "_px_per_bp_" + str(px_per_bp) + "_region_heatmap.tiff"

    # Make the heatmap
    make_heatmap(log2_matrix, max_log2_fc, tick_parameters, output_filename, px_per_bp)
    remove_files(log2_matrix)


def parse_input(args):
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

    parser = argparse.ArgumentParser(prog='PolTools region_fold_change_heatmap')

    parser.add_argument('read_type', metavar='read type', type=str, choices=["five", "three", "whole"],
                        help='either five, three, or whole')

    parser.add_argument('regions_file', metavar='regions_file', type=str,
                        help='Bed formatted file containing all the regions you want to average the sequences')

    parser.add_argument('output_prefix', metavar='output_prefix', type=str, help='Prefix for the output filename')

    parser.add_argument('--numerator', nargs=2, action='append', metavar=('seq_file', 'spike_in'),
                        required=True, help='Provide the sequencing file with its correction factor. You can supply '
                                            'more than one sequencing file by adding multiple --numerator arguments.')

    parser.add_argument('--denominator', nargs=2, action='append', metavar=('seq_file', 'spike_in'),
                        required=True, help='Provide the sequencing file with its correction factor. You can supply '
                                            'more than one sequencing file by adding multiple -denominator arguments.')

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

    parser.add_argument('--minor_ticks', metavar='minor_ticks', dest='minor_ticks',
                        type=positive_int, default=None, help='Distance between minor ticks (bp). Default is 10 bp.')

    parser.add_argument('--major_ticks', metavar='major_ticks', dest='major_ticks',
                        type=positive_int, default=None, help='Distance between major ticks (bp). Default is 100 bp')

    args = parser.parse_args(args)

    verify_bed_files(args.regions_file)
    region_length = determine_region_length(args.regions_file)
    verify_region_length_is_even(region_length)

    filenames = args.regions_file, args.output_prefix
    repeat_amounts = args.repeat_amount, args.vertical_averaging

    numerator_seq_files_data = []
    denominator_seq_files_data = []

    for dataset in args.numerator:
        seq_file, corr_factor = dataset

        corr_factor = positive_float(corr_factor)
        numerator_seq_files_data.append((seq_file, corr_factor))

        if not os.path.isfile(seq_file):
            sys.stderr.write("File " + seq_file + " was not found.\n")
            sys.exit(1)

    for dataset in args.denominator:
        seq_file, corr_factor = dataset

        corr_factor = positive_float(corr_factor)
        denominator_seq_files_data.append((seq_file, corr_factor))

        if not os.path.isfile(seq_file):
            sys.stderr.write("File " + seq_file + " was not found.\n")
            sys.exit(1)

    if args.minor_ticks != None and args.major_ticks != None:
        tick_parameters = (args.minor_ticks * args.repeat_amount, args.major_ticks * args.repeat_amount)
    else:
        tick_parameters = (None, None)

    return args.read_type, filenames, args.max_log2_fc, repeat_amounts, numerator_seq_files_data, denominator_seq_files_data, args.threads, tick_parameters


if __name__ == '__main__':
    run_read_end_fold_change_heatmap(sys.argv[1:])
