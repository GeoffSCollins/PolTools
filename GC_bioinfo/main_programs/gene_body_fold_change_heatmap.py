"""
The goal of this is to make the plot that David described
"""
import glob
import os
import sys
import argparse
import multiprocessing

from PIL import Image

from GC_bioinfo.main_programs import gene_body_combined_heatmap
from GC_bioinfo.utils.constants import generate_heatmap_location
from GC_bioinfo.utils.heatmap_utils.generate_heatmap import generate_heatmap, Ticks, make_ticks_matrix
from GC_bioinfo.utils.heatmap_utils.make_log_two_fold_change_matrix import make_log_two_fold_change_matrix
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.nested_multiprocessing_pool import NestedPool
from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.heatmap_utils.set_matrix_bounds import set_matrix_bounds


def set_max_fold_change(fold_change_matrix_filename, max_fold_change):
    return set_matrix_bounds(fold_change_matrix_filename, -1 * max_fold_change, max_fold_change)


def make_ticks_image(width, interval_size, tick_params):
    minor_ticks_bp, major_ticks_bp = tick_params

    # Make the tick marks
    t = Ticks(minor_tick_mark_interval_size=(minor_ticks_bp / interval_size),
              major_tick_mark_interval_size=(major_ticks_bp / interval_size))

    # Ticks matrix with a height of 50 px and a max black value of 1
    ticks_matrix = make_ticks_matrix(width, 50, 1, t)

    # Write to a file
    ticks_matrix_filename = generate_random_filename()
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


def get_fold_change_matrix(numerator_seq_files_data, denominator_seq_files_data, matrix_params, filenames, max_threads):

    # If max threads is 1, we make the pool have 1 thread and each combined heatmap have one thread
    # If it is 2 or 3, the two combined heatmaps will be made concurrently
    # 4 or more means run at max speed
    if max_threads == 1:
        pool_threads = 1
        comb_heatmap_threads = 1
    elif max_threads == 2 or max_threads == 3:
        pool_threads = 2
        comb_heatmap_threads = 1
    else:
        pool_threads = 2
        comb_heatmap_threads = 2

    numerator_args = (numerator_seq_files_data, matrix_params, filenames, comb_heatmap_threads)
    denominator_args = (denominator_seq_files_data, matrix_params, filenames, comb_heatmap_threads)

    # We use max_threads / 2 because we will be running two instances of the combined
    pool = NestedPool(pool_threads)
    numerator_matrix_filename, denominator_matrix_filename = pool.starmap(gene_body_combined_heatmap.get_combined_matrix,
                                                                          [numerator_args, denominator_args])
    pool.close()

    # Make the fold change matrix
    log_two_fold_change_matrix_filename = make_log_two_fold_change_matrix(numerator_matrix_filename, denominator_matrix_filename)

    remove_files(numerator_matrix_filename, denominator_matrix_filename)

    return log_two_fold_change_matrix_filename


def make_rgb_heatmap(fold_change_matrix_filename, heatmap_params, output_filename_prefix):
    bp_width, width, height, gamma, max_fold_change, interval_size, minor_ticks_bp, major_ticks_bp = heatmap_params

    tick_params = minor_ticks_bp, major_ticks_bp

    only_heatmap_filename = generate_random_filename(extension=".tiff")

    if max_fold_change != None:
        negative_max_fold_change = -1 * max_fold_change
    else:
        negative_max_fold_change = None

    generate_heatmap(fold_change_matrix_filename, 'red/blue', only_heatmap_filename, gamma, negative_max_fold_change,
                     max_fold_change, ticks=None)

    ticks_image_filename = make_ticks_image(width, interval_size, tick_params)

    # Combine the two images together
    output_filename = output_filename_prefix + "_max_" + str(max_fold_change) + "_width_" + str(
        bp_width) + "bp_gene_body_fold_change_heatmap"

    combine_images(ticks_image_filename, only_heatmap_filename, output_filename)

    remove_files(fold_change_matrix_filename,
                 ticks_image_filename, only_heatmap_filename)


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

    parser.add_argument('output_prefix', metavar='output_prefix', type=str, help='Prefix for the output filename')

    parser.add_argument('-u', '--upstream_distance', metavar='upstream_distance', dest='upstream_distance',
                        type=positive_int, default=50_000, help='Distance upstream of the max TSS')

    parser.add_argument('-d', '--distance_past_tes', metavar='distance_past_tes', dest='distance_past_tes',
                        type=positive_int, default=50_000, help='Distance downstream of the transcription end site')

    parser.add_argument('-b', '--bp_width', metavar='bp_width', dest='bp_width', default=400_000, type=positive_int,
                        help='Total number of base pairs shown on the heatmap. This number must be greater than the ' +
                             'upstream distance + distance past TES.')

    parser.add_argument('-w', '--width', metavar='width', dest='width',
                        type=positive_int, default=2_000, help='Width of the heatmap in pixels')

    parser.add_argument('-e', '--height', metavar='height', dest='height',
                        type=positive_int, default=2_000, help='Height of the heatmap in pixels')

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
    upstream_distance = args.upstream_distance
    distance_past_tes = args.distance_past_tes
    bp_width = args.bp_width
    width = args.width
    height = args.height
    max_log2_fc = args.max_log2_fc

    numerator_spike_in_one = args.numerator_correction_factor_one
    numerator_sequencing_filename_one = args.numerator_seq_file_one
    numerator_spike_in_two = args.numerator_correction_factor_two
    numerator_sequencing_filename_two = args.numerator_seq_file_two

    denominator_spike_in_one = args.denominator_correction_factor_one
    denominator_sequencing_filename_one = args.denominator_seq_file_one
    denominator_spike_in_two = args.denominator_correction_factor_two
    denominator_sequencing_filename_two = args.denominator_seq_file_two

    output_filename_prefix = args.output_prefix
    minor_ticks = args.minor_ticks
    major_ticks = args.major_ticks

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

    # We make a gamma variable because the heatmap generation function needs it
    # However, the gamma will never be used
    gamma = 2.2

    numerator_seq_files_data = [(numerator_sequencing_filename_one, numerator_spike_in_one),
                                (numerator_sequencing_filename_two, numerator_spike_in_two)]

    denominator_seq_files_data = [(denominator_sequencing_filename_one, denominator_spike_in_one),
                                (denominator_sequencing_filename_two, denominator_spike_in_two)]

    matrix_params = (upstream_distance, distance_past_tes, width, height, interval_size)
    heatmap_params = (bp_width, width, height, gamma, max_log2_fc, interval_size, minor_ticks, major_ticks)
    filenames = (truQuant_output_file, tsr_file, output_filename_prefix)

    return numerator_seq_files_data, denominator_seq_files_data, matrix_params, heatmap_params, filenames, max_threads


def main(args):
    numerator_seq_files_data, denominator_seq_files_data, matrix_params, heatmap_params, filenames, max_threads = get_args(args)

    log_two_fold_change_matrix_filename = get_fold_change_matrix(
        numerator_seq_files_data, denominator_seq_files_data, matrix_params, filenames, max_threads)

    output_filename_prefix = filenames[-1]

    make_rgb_heatmap(log_two_fold_change_matrix_filename, heatmap_params, output_filename_prefix)


if __name__ == '__main__':
    main(sys.argv[1:])
