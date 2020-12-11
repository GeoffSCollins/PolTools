import sys
import glob
import os
import multiprocessing

from main_programs import gene_body_heatmap

from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.add_matrices import add_matrices
from GC_bioinfo.utils.generate_heatmap import generate_heatmap, Ticks
from GC_bioinfo.utils.average_matrix import average_matrix
from GC_bioinfo.utils.scale_matrix import scale_matrix


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo combine_gene_body_heatmap <truQuant output file> <Upstream Distance>" +
                     " <Distance Past TES> <Width (bp)> <Width (px)> <Height> <Gamma> <Max black value> <Spike in Correction> <Sequencing Filename> <Spike in Correction> <Sequencing Filename> <Output Filename> \n")
    sys.stderr.write \
        ("\nMore information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/combine_gene_body_heatmap.rst\n")


def get_args(args):
    if len(args) != 13:
        print_usage()
        sys.exit(1)

    truQuant_output_file, upstream_distance, distance_past_tes, bp_width, width, height, gamma, max_black_value, spike_in_one, \
    sequencing_filename_one, spike_in_two, sequencing_filename_two, output_filename_prefix = args

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

    return truQuant_output_file, tsr_file, upstream_distance, distance_past_tes, bp_width, width, height, gamma, max_black_value, \
           interval_size, spike_in_one, sequencing_filename_one, spike_in_two, sequencing_filename_two, output_filename_prefix


def get_individual_matrices(truQuant_output_file, distance_past_tes, interval_size, five_prime_buffer_distance,
                            sequencing_filename, blacklist_regions_file, width, height, spike_in):

    # Step 1. Make regions to quantify
    intervals_filename = gene_body_heatmap.make_incremented_regions(truQuant_output_file, distance_past_tes, interval_size,
                                                                    five_prime_buffer_distance)

    # Step 2. Quantify them
    quantified_regions_filename = gene_body_heatmap.quantify_intervals(sequencing_filename, blacklist_regions_file,
                                                                       intervals_filename)

    # Step 3. Read the coverage data and add them to a 2d list
    sorted_matrix_filename, num_lines = gene_body_heatmap.read_coverage_file(quantified_regions_filename, width)

    remove_files(intervals_filename, quantified_regions_filename)

    # Make the dimensions correct
    num_lines_to_average = int(num_lines / height)
    averaged_matrix_filename = average_matrix(sorted_matrix_filename, num_lines_to_average)

    remove_files(sorted_matrix_filename)

    finalized_matrix_filename = scale_matrix(averaged_matrix_filename, spike_in)

    remove_files(averaged_matrix_filename)

    return finalized_matrix_filename


def get_combined_matrix(truQuant_output_file, tsr_file, distance_past_tes, sequencing_filename_one, spike_in_one, sequencing_filename_two, spike_in_two,
                        interval_size, upstream_distance, width, height):

    blacklist_regions_file = gene_body_heatmap.blacklist_extended_gene_bodies(tsr_file, distance_past_tes)

    get_individual_matrices_args = []
    for data_tuple in [(sequencing_filename_one, spike_in_one), (sequencing_filename_two, spike_in_two)]:
        sequencing_filename, spike_in = data_tuple

        get_individual_matrices_args.append(
            [truQuant_output_file, distance_past_tes, interval_size, upstream_distance,
             sequencing_filename, blacklist_regions_file, width, height, spike_in])

    pool = multiprocessing.Pool(processes=2)

    matrix_filenames = pool.starmap(get_individual_matrices, get_individual_matrices_args)

    # Add the matricies together
    combined_matrix = add_matrices(matrix_filenames)

    remove_files(matrix_filenames, blacklist_regions_file)

    return combined_matrix


def main(args):
    truQuant_output_file, tsr_file, upstream_distance, distance_past_tes, bp_width, width, height, gamma, max_black_value, \
        interval_size, spike_in_one, sequencing_filename_one, spike_in_two, sequencing_filename_two, \
        output_filename_prefix = get_args(args)

    combined_matrix = get_combined_matrix(truQuant_output_file, tsr_file, distance_past_tes, sequencing_filename_one, spike_in_one,
                                          sequencing_filename_two, spike_in_two, interval_size,
                                          upstream_distance, width, height)

    # Step 4. Make the heatmap using the resulting file

    # Minor tick marks every 10 kb and major tick marks every 50 kb
    t = Ticks(minor_tick_mark_interval_size=(10_000 / interval_size), major_tick_mark_interval_size=(50_000 / interval_size))

    output_filename = output_filename_prefix + "_gamma_" + str(gamma) + "_max_" + str (max_black_value) + \
                      "_width_" + str(bp_width) + "bp_gene_body_combined_heatmap.tiff"

    generate_heatmap(combined_matrix, 'gray', output_filename, gamma, 0, max_black_value, ticks=t)

    # Step 5. Remove the averaged_matrix file
    remove_files(combined_matrix)


if __name__ == '__main__':
    main(sys.argv[1:])
