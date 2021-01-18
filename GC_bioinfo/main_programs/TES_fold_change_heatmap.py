import glob
import os
import sys
import argparse
import multiprocessing

from GC_bioinfo.utils.constants import generate_heatmap_location
from GC_bioinfo.utils.generate_heatmap import generate_heatmap, Ticks, make_ticks_matrix
from GC_bioinfo.utils.make_log_two_fold_change_matrix import make_log_two_fold_change_matrix
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.nested_multiprocessing_pool import NestedPool
from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.set_matrix_bounds import set_matrix_bounds
from GC_bioinfo.main_programs.TES_combined_heatmap import get_combined_matrix
from GC_bioinfo.main_programs.gene_body_fold_change_heatmap import combine_images


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

    pool = NestedPool(pool_threads)

    numerator_args = (numerator_seq_files_data, matrix_params, filenames, comb_heatmap_threads)
    denominator_args = (denominator_seq_files_data, matrix_params, filenames, comb_heatmap_threads)

    matrices = pool.starmap(get_combined_matrix, [numerator_args, denominator_args])

    pool.close()

    numerator_matrix, denominator_matrix = matrices

    fold_change_matrix = make_log_two_fold_change_matrix(numerator_matrix, denominator_matrix)

    remove_files(numerator_matrix, denominator_matrix)

    return fold_change_matrix


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

        return num

    def positive_float(num):
        try:
            val = float(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")

        return num

    parser = argparse.ArgumentParser(prog='GC_bioinfo TES_fold_change_heatmap',
                                     description="Generate a heatmap of 3' ends for each gene sorted by gene length " +
                                                 "aligned by the transcription end site\n" +
                                                 "More information can be found at " +
                                                 "https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/TES_fold_change_heatmap.rst")

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

    parser.add_argument('-g', '--gamma', metavar='gamma', dest='gamma',
                        type=positive_float, default=2.2, help='Gamma value of the heatmap')

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
    numerator_spike_in_one = args.numerator_spike_in_one
    numerator_sequencing_filename_one = args.numerator_sequencing_filename_one
    numerator_spike_in_two = args.numerator_spike_in_two
    numerator_sequencing_filename_two = args.numerator_sequencing_filename_two
    denominator_spike_in_one = args.denominator_spike_in_one
    denominator_sequencing_filename_one = args.denominator_sequencing_filename_one
    denominator_spike_in_two = args.denominator_spike_in_two
    denominator_sequencing_filename_two = args.denominator_sequencing_filename_two
    output_filename_prefix = args.output_filename_prefix
    width = args.width
    height = args.height
    downstream_distance = args.downstream_distance
    upstream_distance = args.upstream_distance
    bp_width = args.bp_width
    gamma = args.gamma
    max_log2_fc = args.max_log2_fc
    minor_ticks = args.minor_ticks
    major_ticks = args.major_ticks
    max_threads = args.threads

    # Find all regions to blacklist
    tsr_file = glob.glob(truQuant_output_file.replace("gene_body_regions.bed", "") + "*TSR.tab")

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


    # If the interval size is not an integer, then we can't use it
    if bp_width % width:
        sys.stderr.write(
            "The heatmap width in px must be a factor of the base pair width (bp width / px width must be an integer)")
        sys.exit(1)

    interval_size = int(bp_width / width)

    numerator_seq_file_data = [(numerator_sequencing_filename_one, numerator_spike_in_one),
                               (numerator_sequencing_filename_two, numerator_spike_in_two)]
    denominator_seq_file_data = [(denominator_sequencing_filename_one, denominator_spike_in_one),
                               (denominator_sequencing_filename_two, denominator_spike_in_two)]
    matrix_params = (upstream_distance, downstream_distance, bp_width, width, height, interval_size)
    heatmap_params = (bp_width, width, height, gamma, max_log2_fc, interval_size, minor_ticks, major_ticks)
    filenames = (truQuant_output_file, tsr_file, output_filename_prefix)

    return numerator_seq_file_data, denominator_seq_file_data, matrix_params, heatmap_params, filenames, max_threads


def main(args):
    numerator_seq_file_data, denominator_seq_file_data, matrix_params, heatmap_params, filenames, max_threads = get_args(args)

    # Get the fold change matrix
    fold_change_matrix = get_fold_change_matrix(numerator_seq_file_data, denominator_seq_file_data, matrix_params, filenames, max_threads)

    # Now plot!
    bp_width, width, height, gamma, max_log2_fc, interval_size, minor_ticks, major_ticks = heatmap_params
    output_prefix = filenames[-1]

    output_filename = output_prefix + "_max_" + str(max_log2_fc) +  "_width_" + str(bp_width) + \
                      "bp_fold_change_TES_heatmap.tiff"

    only_heatmap_filename = generate_random_filename().replace(".bed", ".tiff")

    generate_heatmap(fold_change_matrix, 'red/blue', only_heatmap_filename, gamma, (-1 * max_log2_fc),
                     max_log2_fc)

    tick_params = (minor_ticks, major_ticks)

    ticks_image_filename = make_ticks_image(width, interval_size, tick_params)

    combine_images(ticks_image_filename, only_heatmap_filename, output_filename)

    remove_files(fold_change_matrix, ticks_image_filename, only_heatmap_filename)


if __name__ == '__main__':
    main(sys.argv[1:])