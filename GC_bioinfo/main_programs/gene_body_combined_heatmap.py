import glob
import multiprocessing
import os
import sys
import argparse

from GC_bioinfo.main_programs import gene_body_heatmap
from GC_bioinfo.utils.add_matrices import add_matrices
from GC_bioinfo.utils.remove_files import remove_files


def get_combined_matrix(seq_files_data, matrix_params, filenames, max_threads):
    args = []
    for seq_data in seq_files_data:
        args.append( (seq_data, matrix_params, filenames) )

    pool = multiprocessing.Pool(processes=max_threads)
    matrix_filenames = pool.starmap(gene_body_heatmap.build_matrix, args)
    pool.close()

    # Add the matricies together
    combined_matrix = add_matrices(matrix_filenames)

    remove_files(matrix_filenames)

    return combined_matrix


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

    parser = argparse.ArgumentParser(prog='GC_bioinfo gene_body_combined_heatmap',
                                     description="Generate a heatmap of 3' ends for each gene sorted by gene length " +
                                                 "aligned by the TSS\n" +
                                                 "More information can be found at " +
                                                 "https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/gene_body_combined_heatmap.rst")

    parser.add_argument('truQuant_output_file', metavar='truQuant_output_file', type=str,
                        help='truQuant output file which ends in -truQuant_output.txt')

    parser.add_argument('correction_factor_one', metavar='correction_factor_one', type=positive_float,
                        help='Correction factor for the first dataset')

    parser.add_argument('seq_file_one', metavar='seq_file_one', type=str,
                        help='First bed formatted sequencing file')

    parser.add_argument('correction_factor_two', metavar='correction_factor_two', type=positive_float,
                        help='Correction factor for the second dataset')

    parser.add_argument('seq_file_two', metavar='seq_file_two', type=str,
                        help='Second bed formatted sequencing file')

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

    parser.add_argument('-g', '--gamma', metavar='gamma', dest='gamma',
                        type=positive_float, default=2.2, help='Gamma value of the heatmap')

    parser.add_argument('-m', '--max_black', metavar='max_black', dest='max_black',
                        type=positive_float, default=None, help='Max black value of the heatmap')

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
    gamma = args.gamma
    max_black_value = args.max_black
    spike_in_one = args.correction_factor_one
    sequencing_filename_one = args.seq_file_one
    spike_in_two = args.correction_factor_two
    sequencing_filename_two = args.seq_file_two
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
        sys.stderr.write("HERE: " + str(tsr_file))
        sys.exit(1)

    tsr_file = tsr_file[0]

    if not os.path.isfile(sequencing_filename_one):
        sys.stderr.write("File " + sequencing_filename_one + " was not found.\n")
        sys.exit(1)

    if not os.path.isfile(sequencing_filename_two):
        sys.stderr.write("File " + sequencing_filename_two + " was not found.\n")
        sys.exit(1)

    if not tsr_file:
        sys.stderr.write("No tsrFinder file was found. Exiting ...\n")
        sys.exit(1)

    if bp_width % width != 0:
        sys.stderr.write("The width (bp) must be evenly divisible by the width (px). Exiting ...")
        sys.exit(1)

    if bp_width < width:
        sys.stderr.write("The width (bp) must be greater than width (px). Exiting ...")
        sys.exit(1)

    interval_size = int(bp_width / width)

    seq_files_data = [(sequencing_filename_one, spike_in_one), (sequencing_filename_two, spike_in_two)]
    matrix_params = (upstream_distance, distance_past_tes, width, height, interval_size)
    heatmap_params = (bp_width, width, height, gamma, max_black_value, interval_size, minor_ticks, major_ticks)
    filenames = (truQuant_output_file, tsr_file, output_filename_prefix)

    return seq_files_data, matrix_params, heatmap_params, filenames, max_threads


def main(args):
    seq_files_data, matrix_params, heatmap_params, filenames, max_threads = get_args(args)

    combined_matrix = get_combined_matrix(seq_files_data, matrix_params, filenames, max_threads)

    output_filename_prefix = filenames[-1]

    # Make the heatmap of the combined matrix
    gene_body_heatmap.make_heatmap(combined_matrix, heatmap_params, output_filename_prefix)

    # Step 5. Remove the averaged_matrix file
    remove_files(combined_matrix)


if __name__ == '__main__':
    main(sys.argv[1:])
