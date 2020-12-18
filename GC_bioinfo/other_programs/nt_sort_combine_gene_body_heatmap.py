import glob
import multiprocessing
import os
import sys

from GC_bioinfo.other_programs import nt_sort_gene_body_heatmap
from GC_bioinfo.utils.add_matrices import add_matrices
from GC_bioinfo.utils.remove_files import remove_files


def get_combined_matrix(seq_files_data, matrix_params, filenames, nt, nt_content_distance, min_gene_length):

    args = []
    for seq_data in seq_files_data:
        args.append( (seq_data, matrix_params, filenames, nt, nt_content_distance, min_gene_length) )

    pool = multiprocessing.Pool(processes=2)
    matrix_filenames = pool.starmap(nt_sort_gene_body_heatmap.build_matrix, args)

    # Add the matricies together
    combined_matrix = add_matrices(matrix_filenames)

    remove_files(matrix_filenames)

    return combined_matrix


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo combine_gene_body_heatmap <truQuant output file> <Upstream Distance>" +
                     " <Distance Past TES> <Width (bp)> <Width (px)> <Height> <Gamma> <Max black value> <Spike in Correction> <Sequencing Filename> <Spike in Correction> <Sequencing Filename> <Output Filename> \n")
    sys.stderr.write \
        ("\nMore information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/combine_gene_body_heatmap.rst\n")


def get_args(args):
    if len(args) != 16:
        print_usage()
        sys.exit(1)

    truQuant_output_file, upstream_distance, distance_past_tes, bp_width, width, height, gamma, max_black_value, spike_in_one, \
    sequencing_filename_one, spike_in_two, sequencing_filename_two, output_filename_prefix, nt, nt_content_distance, min_gene_length = args

    tsr_file = glob.glob(truQuant_output_file.replace("-truQuant_output.txt", "") + "*TSR.tab")

    if not tsr_file:
        sys.stderr.write("No tsrFinder file was found. Exiting ...\n")
        sys.exit(1)

    if len(tsr_file) != 1:
        sys.stderr.write("More than one tsrFinder file was found for this run of truQuant. Exiting ...\n")
        sys.exit(1)

    tsr_file = tsr_file[0]

    if not os.path.isfile(sequencing_filename_one):
        sys.stderr.write("File " + sequencing_filename_one + " was not found.\n")
        sys.exit(1)

    if not os.path.isfile(sequencing_filename_two):
        sys.stderr.write("File " + sequencing_filename_two + " was not found.\n")
        sys.exit(1)


    def try_to_convert_to_int(var, var_name):
        try:
            int_var = int(var)
            return int_var
        except ValueError:
            sys.stderr.write("The " + var_name + " could not be converted to an integer")
            sys.exit(1)

    # Make sure the distance_past_tes, width, max_gene_length are all integers
    upstream_distance = try_to_convert_to_int(upstream_distance, "5' buffer distance")
    distance_past_tes = try_to_convert_to_int(distance_past_tes, "distance past the TES")
    width = try_to_convert_to_int(width, "width (px)")
    bp_width = try_to_convert_to_int(bp_width, "width (bp)")
    height = try_to_convert_to_int(height, "height")

    interval_size = bp_width / width
    interval_size = try_to_convert_to_int(interval_size, "interval size")

    if bp_width % width != 0:
        sys.stderr.write("The width (bp) must be evenly divisible by the width (px). Exiting ...")
        sys.exit(1)

    if bp_width < width:
        sys.stderr.write("The width (bp) must be greater than width (px). Exiting ...")
        sys.exit(1)

    try:
        gamma = float(gamma)
    except ValueError:
        sys.stderr.write("The gamma could not be converted to a float")

    try:
        max_black_value = float(max_black_value)
    except ValueError:
        sys.stderr.write("The gamma could not be converted to a float")

    try:
        spike_in_one = float(spike_in_one)
    except ValueError:
        sys.stderr.write("The spike in correction factor could not be converted to a float")

    try:
        spike_in_two = float(spike_in_two)
    except ValueError:
        sys.stderr.write("The spike in correction factor could not be converted to a float")

    interval_size = try_to_convert_to_int(interval_size, "interval size")

    nt_content_distance = int(nt_content_distance)
    min_gene_length = int(min_gene_length)

    seq_files_data = [(sequencing_filename_one, spike_in_one), (sequencing_filename_two, spike_in_two)]
    matrix_params = (upstream_distance, distance_past_tes, width, height, interval_size)
    heatmap_params = (bp_width, width, height, gamma, max_black_value, interval_size)
    filenames = (truQuant_output_file, tsr_file, output_filename_prefix)

    return seq_files_data, matrix_params, heatmap_params, filenames, nt, nt_content_distance, min_gene_length


def main(args):
    seq_files_data, matrix_params, heatmap_params, filenames, nt, nt_content_distance, min_gene_length = get_args(args)

    combined_matrix = get_combined_matrix(seq_files_data, matrix_params, filenames, nt, nt_content_distance, min_gene_length)

    output_filename_prefix = filenames[-1]

    # Make the heatmap of the combined matrix
    nt_sort_gene_body_heatmap.make_heatmap(combined_matrix, heatmap_params, output_filename_prefix)

    # Step 5. Remove the averaged_matrix file
    remove_files(combined_matrix)


if __name__ == '__main__':
    main(sys.argv[1:])
