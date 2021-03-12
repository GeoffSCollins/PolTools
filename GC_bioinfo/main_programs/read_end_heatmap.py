import sys
import argparse

from GC_bioinfo.utils.get_region_length import determine_region_length
from GC_bioinfo.utils.build_counts_dict import build_counts_dict
from GC_bioinfo.utils.verify_bed_file import verify_bed_files
from GC_bioinfo.utils.verify_region_length_is_even import verify_region_length_is_even
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.heatmap_utils.generate_heatmap import generate_heatmap
from GC_bioinfo.utils.heatmap_utils.average_matrix import average_matrix


def get_original_matrix(regions_filename, sequencing_file, norm_factor, end):
    # 1. Load the five_prime_dict
    counts_dict = build_counts_dict(sequencing_file, end)

    matrix = []

    with open(regions_filename) as file:
        for line in file:
            if line.rstrip():
                chromosome, left, right, gene_name, _, strand = line.split()

                left, right = int(left), int(right)

                curr_gene_counts = [counts_dict[chromosome][strand][i] * norm_factor for i in range(left, right)]

                # If the gene is negative, we need to reverse the list
                if strand == "-":
                    curr_gene_counts = curr_gene_counts[::-1]

                # Now we loop through the region left to right
                matrix.append(curr_gene_counts)

    return matrix


def get_heatmap_matrix(regions_filename, seq_file_data, end, repeat_amounts):
    repeat_amount, vertical_averaging = repeat_amounts
    seq_file, norm_factor = seq_file_data

    # 2. Load 2D list containing the data to be outputted
    original_matrix = get_original_matrix(regions_filename, seq_file, norm_factor, end)

    # Expand the matrix using the repeat amounts and write it to a file
    matrix_filename = generate_random_filename(".matrix")

    with open(matrix_filename, 'w') as file:
        for row in original_matrix:
            # Make the row the correct size by repeating each element by repeat_amount
            output_list = []
            for val in row:
                for _ in range(repeat_amount):
                    output_list.append(str(val))

            file.write(
                "\t".join(output_list) + "\n"
            )

    # Do the vertical averaging
    heatmap_matrix = average_matrix(matrix_filename, vertical_averaging)
    remove_files(matrix_filename)

    return heatmap_matrix


def make_heatmap(matrix_filename, output_prefix, heatmap_parameters):
    max_value, gamma = heatmap_parameters

    # Now make the heatmap
    min_value = None
    output_filename = output_prefix + "_max_" + str(max_value) + ".tiff"

    generate_heatmap(matrix_filename, 'gray', output_filename, gamma, min_value, max_value)


def run_read_end_heatmap(end, filenames, seq_file_data, heatmap_parameters, repeat_amounts):
    regions_filename, output_prefix = filenames

    matrix = get_heatmap_matrix(regions_filename, seq_file_data, end, repeat_amounts)
    make_heatmap(matrix, output_prefix, heatmap_parameters)
    remove_files(matrix)


def parse_input(args):

    def positive_float(num):
        try:
            val = float(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")

        return val


    parser = argparse.ArgumentParser(prog='read_end_heatmap')

    parser.add_argument('read_type', metavar='read type', type=str, choices=["five", "three", "whole"],
                        help='either five, three, or whole')

    parser.add_argument('regions_file', metavar='regions_file', type=str,
                        help='Bed formatted file containing all the regions you want to average the sequences')

    parser.add_argument('seq_file', metavar='sequencing_file', type=str,
                        help='Bed formatted file from the sequencing experiment')

    parser.add_argument('norm_factor', metavar='norm_factor', type=positive_float,
                        help='Correction factor for the sequencing file')

    parser.add_argument('output_prefix', metavar='output_prefix', type=str, help='Prefix for the output filename')

    parser.add_argument('-m', '--max_black', metavar='max_black', dest='max_black',
                        type=int, default=None,
                        help='Max black value of the heatmap. Default is the maximum possible value')

    parser.add_argument('-r', '--repeat_amount', metavar='repeat_amount', dest='repeat_amount',
                        type=int, default=1,
                        help='Each base will be shown in this number of pixels. Default is 1')

    parser.add_argument('-v', '--vertical_averaging', metavar='vertical_averaging', dest='vertical_averaging',
                        type=int, default=1,
                        help='Average this number of rows into one row. Default is 1')

    parser.add_argument('-g', '--gamma', metavar='gamma', dest='gamma',
                        type=positive_float, default=2.2,
                        help='Gamma value of the heatmap. Default is 2.2, which is no gamma correction.')

    args = parser.parse_args(args)

    verify_bed_files(args.regions_file)
    region_length = determine_region_length(args.regions_file)
    verify_region_length_is_even(region_length)

    filenames = args.regions_file, args.seq_file, args.output_prefix
    repeat_amounts = args.repeat_amount, args.vertical_averaging
    heatmap_parameters = args.max_value, args.gamma
    seq_file_data = args.seq_file, args.norm_factor

    return args.read_type, filenames, seq_file_data, heatmap_parameters, repeat_amounts


if __name__ == '__main__':
    end, filenames, seq_file_data, heatmap_parameters, repeat_amounts = parse_input(sys.argv[1:])
    run_read_end_heatmap(end, filenames, seq_file_data, heatmap_parameters, repeat_amounts)
