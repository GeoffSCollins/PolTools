import glob
import os
import sys

from GC_bioinfo.utils.constants import generate_heatmap_location
from GC_bioinfo.utils.generate_heatmap import generate_heatmap, Ticks, make_ticks_matrix
from GC_bioinfo.utils.make_fold_change_matrix import make_fold_change_matrix
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.nested_multiprocessing_pool import NestedPool
from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.set_matrix_bounds import set_matrix_bounds
from main_programs.combine_TES_heatmap import get_combined_matrix
from main_programs.gene_body_fold_change_heatmap import combine_images


def get_fold_change_matrix(tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gene_body_file, interval_size, width, height, numerator_sequencing_filename_one,
                           numerator_spike_in_one, numerator_sequencing_filename_two, numerator_spike_in_two,
                           denominator_sequencing_filename_one, denominator_spike_in_one, denominator_sequencing_filename_two, denominator_spike_in_two):

    pool = NestedPool(processes=2)

    numerator_args = (tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gene_body_file, interval_size, width, height,
                                           numerator_sequencing_filename_one, numerator_spike_in_one,
                                           numerator_sequencing_filename_two, numerator_spike_in_two)
    denominator_args = (tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gene_body_file, interval_size, width, height,
                                             denominator_sequencing_filename_one, denominator_spike_in_one,
                                             denominator_sequencing_filename_two, denominator_spike_in_two)

    matrices = pool.starmap(get_combined_matrix, [numerator_args, denominator_args])

    numerator_matrix, denominator_matrix = matrices

    fold_change_matrix = make_fold_change_matrix(numerator_matrix, denominator_matrix)

    remove_files(numerator_matrix, denominator_matrix)

    return fold_change_matrix


def set_max_fold_change(fold_change_matrix_filename, max_fold_change):
    return set_matrix_bounds(fold_change_matrix_filename, -1 * max_fold_change, max_fold_change)


def make_ticks_image(width, interval_size):
    # Make the tick marks
    # Minor tick marks every 10 kb and major tick marks every 50 kb
    t = Ticks(minor_tick_mark_interval_size=(10_000 / interval_size),
              major_tick_mark_interval_size=(50_000 / interval_size))

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


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo TES_fold_change_heatmap <truQuant output file> <Width> <Height> <Downstream Distance> <Upstream Distance> <Distance Upstream of Gene Body Start> <Gamma> <Max fold change> "
                     "<Numerator Spike in Correction> <Numerator Sequencing Filename> <Numerator Spike in Correction> <Numerator Sequencing Filename> "
                     " <Denominator Spike in Correction> <Denominator Sequencing Filename> <Denominator Spike in Correction> <Denominator Sequencing Filename> <Output Filename> \n")
    sys.stderr.write("\nMore information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/gene_body_heatmap.rst\n")


def get_args(args):
    if len(args) != 17:
        print_usage()
        sys.exit(1)

    truQuant_output_file, width, height, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gamma, max_fold_change, numerator_spike_in_one, \
        numerator_sequencing_filename_one, numerator_spike_in_two, numerator_sequencing_filename_two, \
        denominator_spike_in_one, denominator_sequencing_filename_one, denominator_spike_in_two, denominator_sequencing_filename_two, output_prefix = args

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

    if not os.path.isfile(numerator_sequencing_filename_one):
        sys.stderr.write("File " + numerator_sequencing_filename_one + " was not found.\n")
        sys.exit(1)

    if not os.path.isfile(numerator_sequencing_filename_two):
        sys.stderr.write("File " + numerator_sequencing_filename_two + " was not found.\n")
        sys.exit(1)

    if not os.path.isfile(denominator_sequencing_filename_one):
        sys.stderr.write("File " + denominator_sequencing_filename_one + " was not found.\n")
        sys.exit(1)

    if not os.path.isfile(denominator_sequencing_filename_two):
        sys.stderr.write("File " + denominator_sequencing_filename_two + " was not found.\n")
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
        sys.stderr.write("The heatmap width in px must be a factor of the base pair width (bp width / px width must be an integer)")
        sys.exit(1)

    try:
        gamma = float(gamma)
    except ValueError:
        sys.stderr.write("The gamma could not be converted to a float")

    try:
        max_fold_change = float(max_fold_change)
    except ValueError:
        sys.stderr.write("The gamma could not be converted to a float")

    try:
        numerator_spike_in_one = float(numerator_spike_in_one)
    except ValueError:
        sys.stderr.write("The spike in correction factor could not be converted to a float")

    try:
        numerator_spike_in_two = float(numerator_spike_in_two)
    except ValueError:
        sys.stderr.write("The spike in correction factor could not be converted to a float")

    try:
        denominator_spike_in_one = float(denominator_spike_in_one)
    except ValueError:
        sys.stderr.write("The spike in correction factor could not be converted to a float")

    try:
        denominator_spike_in_two = float(denominator_spike_in_two)
    except ValueError:
        sys.stderr.write("The spike in correction factor could not be converted to a float")

    return truQuant_output_file, tsr_file, width, height, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gamma, max_fold_change, interval_size, numerator_spike_in_one, numerator_sequencing_filename_one, numerator_spike_in_two, numerator_sequencing_filename_two, \
            denominator_spike_in_one, denominator_sequencing_filename_one, denominator_spike_in_two, denominator_sequencing_filename_two, output_prefix


def main(args):
    truQuant_output_file, tsr_file, width, height, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, gamma, max_fold_change, interval_size, numerator_spike_in_one, numerator_sequencing_filename_one, numerator_spike_in_two, numerator_sequencing_filename_two, \
    denominator_spike_in_one, denominator_sequencing_filename_one, denominator_spike_in_two, denominator_sequencing_filename_two, output_prefix = get_args(args)

    # Get the fold change matrix
    fold_change_matrix = get_fold_change_matrix(tsr_file, downstream_distance, upstream_distance, distance_upstream_of_gene_body_start, truQuant_output_file, interval_size, width, height, numerator_sequencing_filename_one,
                           numerator_spike_in_one, numerator_sequencing_filename_two, numerator_spike_in_two,
                           denominator_sequencing_filename_one, denominator_spike_in_one, denominator_sequencing_filename_two, denominator_spike_in_two)

    # Now plot!
    output_filename = output_prefix + "_max_" + str(max_fold_change) +  "_width_" + str(
        downstream_distance + upstream_distance) + "bp_fold_change_TES_heatmap.tiff"

    only_heatmap_filename = generate_random_filename().replace(".bed", ".tiff")

    generate_heatmap(fold_change_matrix, 'red/blue', only_heatmap_filename, gamma, (-1 * max_fold_change),
                     max_fold_change)

    ticks_image_filename = make_ticks_image(width, interval_size)

    combine_images(ticks_image_filename, only_heatmap_filename, output_filename)

    remove_files(fold_change_matrix, ticks_image_filename, only_heatmap_filename)


if __name__ == '__main__':
    main(sys.argv[1:])