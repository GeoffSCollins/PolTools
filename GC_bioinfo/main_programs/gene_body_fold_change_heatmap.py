"""
The goal of this is to make the plot that David described
"""
import sys
import os
import glob
import math

from PIL import Image

from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.generate_heatmap import generate_heatmap, Ticks, make_ticks_matrix
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.constants import generate_heatmap_location
from GC_bioinfo.utils.nested_multiprocessing_pool import NestedPool
from GC_bioinfo.utils.set_matrix_bounds import set_matrix_bounds

from main_programs import combine_gene_body_heatmap


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo gene_body_fold_change_heatmap <truQuant output file> <Upstream Distance>" +
                     " <Distance Past TES> <Width (bp)> <Width (px)> <Height> <Gamma> <Max fold change> <Spike in Correction> <Sequencing Filename>" +
                     "<Numerator Spike in Correction> <Numerator Sequencing Filename> <Numerator Spike in Correction> <Numerator Sequencing Filename>" +
                     "<Denomenator Spike in Correction> <Denomenator Sequencing Filename> <Denomenator Spike in Correction> <Denomenator Sequencing Filename> <Output Filename> \n")
    sys.stderr.write("\nMore information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/gene_body_fold_change_heatmap.rst\n")


def get_args(args):
    if len(args) != 17:
        print(len(args))
        print_usage()
        sys.exit(1)

    truQuant_output_file, upstream_distance, distance_past_tes, bp_width, width, height, gamma, max_fold_change, \
    numerator_spike_in_one, numerator_sequencing_filename_one, numerator_spike_in_two, numerator_sequencing_filename_two,\
    denominator_spike_in_one, denominator_sequencing_filename_one, denominator_spike_in_two, denominator_sequencing_filename_two, output_filename_prefix = args

    tsr_file = glob.glob(truQuant_output_file.replace("-truQuant_output.txt", "") + "*TSR.tab")

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

    def try_to_convert_to_float(var, var_name):
        try:
            int_var = float(var)
            return int_var
        except ValueError:
            sys.stderr.write("The " + var_name + " could not be converted to a float")
            sys.exit(1)


    gamma = try_to_convert_to_float(gamma, "gamma")
    max_fold_change = try_to_convert_to_float(max_fold_change, "max fold change")
    numerator_spike_in_one = try_to_convert_to_float(numerator_spike_in_one, "numerator spike in one")
    numerator_spike_in_two = try_to_convert_to_float(numerator_spike_in_two, "numerator spike in two")
    denominator_spike_in_one = try_to_convert_to_float(denominator_spike_in_one, "denominator spike in one")
    denominator_spike_in_two = try_to_convert_to_float(denominator_spike_in_two, "denominator spike in two")

    return truQuant_output_file, tsr_file, upstream_distance, distance_past_tes, bp_width, width, height, gamma, max_fold_change, \
        interval_size, numerator_spike_in_one, numerator_sequencing_filename_one, numerator_spike_in_two, numerator_sequencing_filename_two, \
        denominator_spike_in_one, denominator_sequencing_filename_one, denominator_spike_in_two, denominator_sequencing_filename_two, output_filename_prefix


def make_fold_change_matrix(numerator_filename, denominator_filename):
    # Going to divide the first by the second
    first_matrix = []
    with open(numerator_filename) as file:
        for line in file:
            first_matrix.append([float(val) for val in line.rstrip().split()])

    second_matrix = []
    with open(denominator_filename) as file:
        for line in file:
            second_matrix.append([float(val) for val in line.rstrip().split()])

    # Verify the matricies are the same size
    first_matrix_width = len(first_matrix[0])
    first_matrix_height = len(first_matrix)

    second_matrix_width = len(second_matrix[0])
    second_matrix_height = len(second_matrix)

    if not (first_matrix_width == second_matrix_width) and not (first_matrix_height == second_matrix_height):
        sys.stderr.write("The matricies are not the same size. Exiting ...")
        sys.exit(1)

    # Now divide the first by the second
    fold_change_matrix = []
    for i in range(first_matrix_height):
        fold_change_matrix.append([0 for _ in range(first_matrix_width)])

    for row in range(first_matrix_height):
        for col in range(first_matrix_width):
            numerator = first_matrix[row][col]
            denominator = second_matrix[row][col]

            if numerator == 0:
                numerator = 1

            if denominator == 0:
                denominator = 1

            value = numerator / denominator

            # Take the log2 of this value
            log_two_fold_change = math.log2(value)

            fold_change_matrix[row][col] = log_two_fold_change

    # Write to a file
    fold_change_matrix_filename = generate_random_filename()
    with open(fold_change_matrix_filename, 'w') as file:
        for row in fold_change_matrix:
            str_row = [str(val) for val in row]
            file.write("\t".join(str_row) + "\n")

    return fold_change_matrix_filename


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


def combine_images(ticks_image_filename, only_heatmap_filename, output_filename):
    ticks_image = Image.open(ticks_image_filename)
    heatmap_image = Image.open(only_heatmap_filename)

    final_image = Image.new('RGB', (ticks_image.width, ticks_image.height + heatmap_image.height))

    final_image.paste(heatmap_image, (0, 0))
    final_image.paste(ticks_image, (0, heatmap_image.height))

    final_image.save(output_filename + ".tiff")


def main(args):
    truQuant_output_file, tsr_file, upstream_distance, distance_past_tes, bp_width, width, height, gamma, max_fold_change, \
    interval_size, numerator_spike_in_one, numerator_sequencing_filename_one, numerator_spike_in_two, numerator_sequencing_filename_two, \
    denominator_spike_in_one, denominator_sequencing_filename_one, denominator_spike_in_two, denominator_sequencing_filename_two, output_filename_prefix = get_args(args)

    # Also run the quantification at the same time
    # p = subprocess.Popen("GC_bioinfo quantify_gene_body_fold_change_heatmap " + " ".join(args), shell=True)

    numerator_args = (truQuant_output_file, tsr_file, distance_past_tes, numerator_sequencing_filename_one, numerator_spike_in_one,
                                          numerator_sequencing_filename_two, numerator_spike_in_two, interval_size,
                                          upstream_distance, width, height)

    denominator_args = (truQuant_output_file, tsr_file, distance_past_tes, denominator_sequencing_filename_one, denominator_spike_in_one,
                                          denominator_sequencing_filename_two, denominator_spike_in_two, interval_size,
                                          upstream_distance, width, height)

    pool = NestedPool(2)
    numerator_matrix_filename, denominator_matrix_filename = pool.starmap(combine_gene_body_heatmap.get_combined_matrix, [numerator_args, denominator_args])

    # Make the fold change matrix

    fold_change_matrix_filename = make_fold_change_matrix(numerator_matrix_filename, denominator_matrix_filename)

    # Step 4. Make the heatmap using the resulting file

    only_heatmap_filename = generate_random_filename().replace(".bed", ".tiff")

    generate_heatmap(fold_change_matrix_filename, 'red/blue', only_heatmap_filename, gamma, (-1 * max_fold_change), max_fold_change, ticks=None)

    ticks_image_filename = make_ticks_image(width, interval_size)

    # Combine the two images together
    output_filename = output_filename_prefix + "_max_" + str(max_fold_change) + "_width_" + str(bp_width) + "bp_gene_body_fold_change_heatmap"

    combine_images(ticks_image_filename, only_heatmap_filename, output_filename)

    remove_files(numerator_matrix_filename, denominator_matrix_filename, fold_change_matrix_filename,
                 ticks_image_filename, only_heatmap_filename)


if __name__ == '__main__':
    main(sys.argv[1:])
