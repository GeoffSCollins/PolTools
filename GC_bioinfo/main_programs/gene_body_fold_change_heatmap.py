"""
The goal of this is to make the plot that David described
"""
import glob
import os
import sys

from PIL import Image

from GC_bioinfo.main_programs import combine_gene_body_heatmap
from GC_bioinfo.utils.constants import generate_heatmap_location
from GC_bioinfo.utils.generate_heatmap import generate_heatmap, Ticks, make_ticks_matrix
from GC_bioinfo.utils.make_fold_change_matrix import make_fold_change_matrix
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.nested_multiprocessing_pool import NestedPool
from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.set_matrix_bounds import set_matrix_bounds


def set_max_fold_change(fold_change_matrix_filename, max_fold_change):
    return set_matrix_bounds(fold_change_matrix_filename, -1 * max_fold_change, max_fold_change)


def make_ticks_image(width, interval_size, tick_params):
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


def get_fold_change_matrix(numerator_seq_files_data, denominator_seq_files_data, matrix_params, filenames):
    numerator_args = (numerator_seq_files_data, matrix_params, filenames)

    denominator_args = (denominator_seq_files_data, matrix_params, filenames)

    pool = NestedPool(2)
    numerator_matrix_filename, denominator_matrix_filename = pool.starmap(combine_gene_body_heatmap.get_combined_matrix,
                                                                          [numerator_args, denominator_args])

    # Make the fold change matrix
    fold_change_matrix_filename = make_fold_change_matrix(numerator_matrix_filename, denominator_matrix_filename)

    remove_files(numerator_matrix_filename, denominator_matrix_filename)

    return fold_change_matrix_filename


def make_rgb_heatmap(fold_change_matrix_filename, heatmap_params, tick_params, output_filename_prefix):
    bp_width, width, height, gamma, max_fold_change, interval_size = heatmap_params

    only_heatmap_filename = generate_random_filename(extension=".tiff")

    generate_heatmap(fold_change_matrix_filename, 'red/blue', only_heatmap_filename, gamma, (-1 * max_fold_change),
                     max_fold_change, ticks=None)

    ticks_image_filename = make_ticks_image(width, interval_size, tick_params)

    # Combine the two images together
    output_filename = output_filename_prefix + "_max_" + str(max_fold_change) + "_width_" + str(
        bp_width) + "bp_gene_body_fold_change_heatmap"

    combine_images(ticks_image_filename, only_heatmap_filename, output_filename)

    remove_files(fold_change_matrix_filename,
                 ticks_image_filename, only_heatmap_filename)


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

    for file in [numerator_sequencing_filename_one, numerator_sequencing_filename_two,
                 denominator_sequencing_filename_one, denominator_sequencing_filename_two]:
        if not os.path.isfile(file):
            sys.stderr.write("File " + file + " was not found.\n")
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



    numerator_seq_files_data = [(numerator_sequencing_filename_one, numerator_spike_in_one),
                                (numerator_sequencing_filename_two, numerator_spike_in_two)]

    denominator_seq_files_data = [(denominator_sequencing_filename_one, denominator_spike_in_one),
                                (denominator_sequencing_filename_two, denominator_spike_in_two)]

    matrix_params = (upstream_distance, distance_past_tes, width, height, interval_size)
    heatmap_params = (bp_width, width, height, gamma, max_fold_change, interval_size)
    filenames = (truQuant_output_file, tsr_file, output_filename_prefix)


    return numerator_seq_files_data, denominator_seq_files_data, matrix_params, heatmap_params, filenames


def main(args):
    numerator_seq_files_data, denominator_seq_files_data, matrix_params, heatmap_params, filenames = get_args(args)

    fold_change_matrix_filename = get_fold_change_matrix(
        numerator_seq_files_data, denominator_seq_files_data, matrix_params, filenames)

    output_filename_prefix = filenames[-1]

    minor_ticks = 10_000  # Minor ticks every 10kb
    major_ticks = 50_000  # Minor ticks every 10kb
    tick_params = (minor_ticks, major_ticks)

    make_rgb_heatmap(fold_change_matrix_filename, heatmap_params, tick_params, output_filename_prefix)


if __name__ == '__main__':
    main(sys.argv[1:])
