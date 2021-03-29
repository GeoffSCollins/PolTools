import glob
import os
import sys
import argparse
import multiprocessing

from PolTools.main_programs import TES_heatmap

from PolTools.utils.constants import generate_heatmap_location
from PolTools.utils.heatmap_utils.generate_heatmap import generate_heatmap, Ticks, make_ticks_matrix
from PolTools.utils.heatmap_utils.make_log_two_fold_change_matrix import make_log_two_fold_change_matrix
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.nested_multiprocessing_pool import NestedPool
from PolTools.utils.remove_files import remove_files
from PolTools.utils.heatmap_utils.set_matrix_bounds import set_matrix_bounds
from PolTools.main_programs.gene_body_fold_change_heatmap import combine_images


def get_fold_change_matrix(numerator_seq_files_data, denominator_seq_files_data, matrix_params, filenames, max_threads):
    # We use max_threads / 2 because we will be running two instances of the combined
    threads_per_heatmap = int(max_threads / 2)

    # Make sure that if the user only wants to run on one thread that it does not default to 0
    if threads_per_heatmap == 0:
        threads_per_heatmap = 1

    numerator_args = (numerator_seq_files_data, matrix_params, filenames, threads_per_heatmap)
    denominator_args = (denominator_seq_files_data, matrix_params, filenames, threads_per_heatmap)

    with NestedPool(max_threads) as pool:
        numerator_matrix_filename, denominator_matrix_filename = pool.starmap(TES_heatmap.get_matrix,
                                                                              [numerator_args, denominator_args])

    # Make the fold change matrix
    log_two_fold_change_matrix_filename = make_log_two_fold_change_matrix(numerator_matrix_filename,
                                                                          denominator_matrix_filename)

    remove_files(numerator_matrix_filename, denominator_matrix_filename)

    return log_two_fold_change_matrix_filename


def set_max_fold_change(fold_change_matrix_filename, max_fold_change):
    return set_matrix_bounds(fold_change_matrix_filename, -1 * max_fold_change, max_fold_change)


def make_ticks_image(width, interval_size, tick_params):
    minor_ticks_bp, major_ticks_bp = tick_params

    # Make the tick marks
    t = Ticks(minor_tick_mark_interval_size=(minor_ticks_bp / interval_size),
              major_tick_mark_interval_size=(major_ticks_bp / interval_size))

    ticks_matrix = make_ticks_matrix(width, 50, 1, t)

    # Write to a file
    ticks_matrix_filename = generate_random_filename()
    with open(ticks_matrix_filename, 'w') as file:
        for row in ticks_matrix:
            file.write("\t".join([str(val) for val in row]) + "\n")

    ticks_image_filename = generate_random_filename().replace(".bed", ".tiff")

    os.system("/usr/bin/Rscript " + generate_heatmap_location + " " +
              " ".join([ticks_matrix_filename, "gray", ticks_image_filename, "2.2"]))

    remove_files(ticks_matrix_filename)

    return ticks_image_filename


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

    parser = argparse.ArgumentParser(prog='PolTools TES_fold_change_heatmap',
                                     description="Generate a heatmap of 3' ends for each gene sorted by gene length " +
                                                 "aligned by the transcription end site\n" +
                                                 "More information can be found at " +
                                                 "https://geoffscollins.github.io/PolTools/TES_fold_change_heatmap.html")

    parser.add_argument('truQuant_output_file', metavar='truQuant_output_file', type=str,
                        help='truQuant output file which ends in -truQuant_output.txt')

    parser.add_argument('--numerator', nargs=2, action='append', metavar=('seq_file', 'spike_in'),
                        required=True, help='Provide the sequencing file with its correction factor. You can supply '
                                            'more than one sequencing file by adding multiple --numerator arguments.')

    parser.add_argument('--denominator', nargs=2, action='append', metavar=('seq_file', 'spike_in'),
                        required=True, help='Provide the sequencing file with its correction factor. You can supply '
                                            'more than one sequencing file by adding multiple -denominator arguments.')

    parser.add_argument('output_prefix', metavar='output_prefix', type=str, help='Prefix for the output filename')

    parser.add_argument('-w', '--width', metavar='width', dest='width',
                        type=positive_int, default=2_000, help='Width of the heatmap in pixels')

    parser.add_argument('-e', '--height', metavar='height', dest='height',
                        type=positive_int, default=2_000, help='Height of the heatmap in pixels')

    parser.add_argument('-d', '--downstream_distance', metavar='downstream_distance', dest='downstream_distance',
                        type=positive_int, default=50_000, help='Distance downstream from the transcription end site')

    parser.add_argument('-u', '--upstream_distance', metavar='upstream_distance', dest='upstream_distance',
                        type=positive_int, default=50_000, help='Distance upstream of the start of the gene body')

    parser.add_argument('-b', '--bp_width', metavar='bp_width', dest='bp_width', default=400_000, type=positive_int,
                        help='Total number of base pairs shown on the heatmap. This number must be greater than the ' +
                             'upstream distance + distance past TES.')

    parser.add_argument('-m', '--max_log2_fc', metavar='max_log2_fc', dest='max_log2_fc',
                        type=positive_float, default=None, help='Max log2 fold change of the heatmap')

    parser.add_argument('--minor_ticks', metavar='minor_ticks', dest='minor_ticks',
                        type=positive_int, default=10_000, help='Distance between minor ticks (bp)')

    parser.add_argument('--major_ticks', metavar='major_ticks', dest='major_ticks',
                        type=positive_int, default=50_000, help='Distance between major ticks (bp)')

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
    max_log2_fc = args.max_log2_fc
    minor_ticks = args.minor_ticks
    major_ticks = args.major_ticks
    max_threads = args.threads

    # Find all regions to blacklist
    tsr_file = glob.glob(truQuant_output_file.replace("-truQuant_output.txt", "") + "*TSR.tab")

    if not tsr_file:
        sys.stderr.write("No tsrFinder file was found. Exiting ...\n")
        sys.exit(1)

    if len(tsr_file) != 1:
        sys.stderr.write("More than one tsrFinder file was found for this run of truQuant. Exiting ...\n")
        sys.exit(1)

    tsr_file = tsr_file[0]

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


    # If the interval size is not an integer, then we can't use it
    if bp_width % width:
        sys.stderr.write(
            "The heatmap width in px must be a factor of the base pair width (bp width / px width must be an integer)")
        sys.exit(1)

    interval_size = int(bp_width / width)

    matrix_params = (upstream_distance, downstream_distance, bp_width, width, height, interval_size)
    heatmap_params = (bp_width, width, height, max_log2_fc, interval_size, minor_ticks, major_ticks)
    filenames = (truQuant_output_file, tsr_file, output_filename_prefix)

    return numerator_seq_files_data, denominator_seq_files_data, matrix_params, heatmap_params, filenames, max_threads


def main(args):
    numerator_seq_files_data, denominator_seq_files_data, matrix_params, heatmap_params, filenames, max_threads = get_args(args)

    # Get the fold change matrix
    fold_change_matrix = get_fold_change_matrix(numerator_seq_files_data, denominator_seq_files_data, matrix_params, filenames, max_threads)

    # Now plot!
    bp_width, width, height, max_log2_fc, interval_size, minor_ticks, major_ticks = heatmap_params
    output_prefix = filenames[-1]

    output_filename = output_prefix + "_max_" + str(max_log2_fc) +  "_width_" + str(bp_width) + \
                      "bp_fold_change_TES_heatmap"

    only_heatmap_filename = generate_random_filename(".tiff")

    negative_log2_value = -1 * max_log2_fc if max_log2_fc else None

    generate_heatmap(fold_change_matrix, 'red/blue', only_heatmap_filename, 2.2, negative_log2_value, max_log2_fc)

    tick_params = (minor_ticks, major_ticks)

    ticks_image_filename = make_ticks_image(width, interval_size, tick_params)

    combine_images(ticks_image_filename, only_heatmap_filename, output_filename)

    remove_files(fold_change_matrix, ticks_image_filename, only_heatmap_filename)


if __name__ == '__main__':
    main(sys.argv[1:])