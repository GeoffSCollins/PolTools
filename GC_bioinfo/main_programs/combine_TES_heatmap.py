import sys
import glob
import os
import multiprocessing

from main_programs.TES_heatmap import get_matrix

from GC_bioinfo.utils.add_matrices import add_matrices
from GC_bioinfo.utils.generate_heatmap import generate_heatmap, Ticks
from GC_bioinfo.utils.remove_files import remove_files


def get_matrices(tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gene_body_file, interval_size, width, height, sequencing_filename_one, spike_in_one, sequencing_filename_two, spike_in_two):
    pool = multiprocessing.Pool(processes=2)

    matrix_one_args = (tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gene_body_file, interval_size, sequencing_filename_one, width, height, spike_in_one)
    matrix_two_args = (tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gene_body_file, interval_size, sequencing_filename_two, width, height, spike_in_two)

    matrix_filenames = pool.starmap(get_matrix, [matrix_one_args, matrix_two_args])

    return matrix_filenames


def get_combined_matrix(tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gene_body_file, interval_size, width, height, sequencing_filename_one,
                        spike_in_one, sequencing_filename_two, spike_in_two):
    # First get the two matricies
    matrix_filename_one, matrix_filename_two = get_matrices(tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gene_body_file, interval_size, width,
                                                            height, sequencing_filename_one, spike_in_one,
                                                            sequencing_filename_two, spike_in_two)

    # Now combine them
    combined_matrix = add_matrices([matrix_filename_one, matrix_filename_two])

    remove_files(matrix_filename_one, matrix_filename_two)

    return combined_matrix


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo combine_TES_heatmap <truQuant output file> <Width> <Height> <Downstream Distance> <Upstream Distance> <Distance Upstream of Gene Body Start> "
                     "<Gamma> <Max black value> <Spike in Correction> <Sequencing Filename> <Spike in Correction> <Sequencing Filename> <Output Filename> \n")
    sys.stderr.write("\nMore information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/gene_body_heatmap.rst\n")


def get_args(args):
    if len(args) != 13:
        print_usage()
        sys.exit(1)

    truQuant_output_file, width, height, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gamma, max_black_value, spike_in_one, sequencing_filename_one, spike_in_two, sequencing_filename_two, output_prefix = args

    gene_body_file = glob.glob(truQuant_output_file.replace("-truQuant_output.txt", "") + "*gene_body_regions.bed")

    if not gene_body_file:
        sys.stderr.write("No gene body file was found for that run of truQuant. Exiting ...\n")
        sys.exit(1)

    if len(gene_body_file) != 1:
        sys.stderr.write("More than one gene body file was found. Exiting ...\n")
        sys.exit(1)

    gene_body_file = gene_body_file[0]

    # Find all regions to blacklist
    tsr_file = glob.glob(gene_body_file.replace("gene_body_regions.bed", "") + "*TSR.tab")

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


    width = try_to_convert_to_int(width, "width")
    height = try_to_convert_to_int(height, "height")
    downstream_distance = try_to_convert_to_int(downstream_distance, "downstream distance")
    upstream_distance = try_to_convert_to_int(upstream_distance, "upstream distance")
    distance_upstream_of_gene_body_start = try_to_convert_to_int(distance_upstream_of_gene_body_start, "distance upstream of gene body start")
    interval_size = int((downstream_distance + upstream_distance) / width)

    # If the interval size is not an integer, then we can't use it
    if (downstream_distance + upstream_distance) % width:
        sys.stderr.write(
            "The heatmap width in px must be a factor of the base pair width (bp width / px width must be an integer)")
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

    return truQuant_output_file, tsr_file, width, height, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gamma, max_black_value, interval_size, spike_in_one, sequencing_filename_one, spike_in_two, sequencing_filename_two, output_prefix


def main(args):
    truQuant_output_file, tsr_file, width, height, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gamma, max_black_value, interval_size, spike_in_one, sequencing_filename_one, spike_in_two, sequencing_filename_two, output_prefix = get_args(args)

    combined_matrix = get_combined_matrix(tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, truQuant_output_file, interval_size, width,
                                            height, sequencing_filename_one, spike_in_one, sequencing_filename_two, spike_in_two)

    # Now plot!
    output_filename = output_prefix + "_max_" + str(max_black_value) + "_gamma_" + str(gamma) + "_width_" + str(
        downstream_distance + upstream_distance) + "bp_combined_TES_heamap.tiff"

    t = Ticks(10_000 / interval_size, 50_000 / interval_size)
    generate_heatmap(combined_matrix, "gray", output_filename, gamma, 0, max_black_value, t)

    remove_files(combined_matrix)


if __name__ == '__main__':
    main(sys.argv[1:])