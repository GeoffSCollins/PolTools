import sys
import argparse

from PolTools.utils.verify_bed_file import verify_bed_files
from PolTools.utils.get_region_length import determine_region_length
from PolTools.utils.bedtools_utils.run_bedtools_getfasta import run_getfasta
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.heatmap_utils.average_matrix import average_matrix
from PolTools.utils.heatmap_utils.generate_heatmap import generate_heatmap
from PolTools.utils.remove_files import remove_files


def expand_region(max_tss_file, region_width):
    expanded_regions_file = generate_random_filename()

    with open(max_tss_file) as file:
        with open(expanded_regions_file, 'w') as output_file:
            for line in file:
                chromosome, left, right, name, score, strand = line.split()

                # Add region_width / 2 to both sides
                region_left = str(int(left) - int(region_width / 2))
                region_right = str(int(left) + int(region_width / 2))

                # If the region is on the negative strand, we move right one because we use the right position
                if strand == "-":
                    region_left = str(int(region_left) + 1)
                    region_right = str(int(region_right) + 1)

                output_file.write(
                    "\t".join([chromosome, region_left, region_right, name, score, strand]) + "\n"
                )

    return expanded_regions_file


def get_sequences(regions_file):
    fasta_file = run_getfasta(regions_file)

    sequences = []

    with open(fasta_file) as file:
        for i, line in enumerate(file):
            if i % 2 == 0:
                # This line has the >
                pass
            else:
                sequences.append(line.rstrip())

    remove_files(fasta_file)

    return sequences


def convert_sequences_to_matrix_file(sequences, nucleotide, heatmap_width):
    matrix_file = generate_random_filename(".matrix")

    # We want to expand each nucleotide a certain amount of times depending on the desired pixel width of the heatmap
    repeat_amount = int(heatmap_width / len(sequences[0]))

    # The repeat amount must be at least one to generate a heatmap
    if repeat_amount < 1:
        repeat_amount = 1

    with open(matrix_file, 'w') as file:
        for seq in sequences:
            curr_seq_list = []

            for nt in seq:
                if nt == nucleotide:
                    for _ in range(repeat_amount):
                        curr_seq_list.append("1")
                else:
                    for _ in range(repeat_amount):
                        curr_seq_list.append("0")

            file.write(
                "\t".join(curr_seq_list) + "\n"
            )

    return matrix_file


def parse_args(args):
    def positive_int(value):
        try:
            val = int(value)
            if val < 0:
                raise argparse.ArgumentTypeError(value + " must be a postive integer")
        except:
            raise argparse.ArgumentTypeError(value + " must be a postive integer")

        return val

    parser = argparse.ArgumentParser(prog='PolTools nucleotide_heatmap',
                                     description='Create a grayscale heatmap for each nucleotide.\n' +
                                     "More information can be found at " +
                                     "https://geoffscollins.github.io/PolTools/nucleotide_heatmap.html")

    parser.add_argument('max_tss_file', metavar='max_tss_file', type=str,
                        help='Bed formatted file which has the base of the max TSS')

    parser.add_argument('region_width', metavar='region_width', type=positive_int, help='Base pair width of the heatmap')

    parser.add_argument('heatmap_width', metavar='heatmap_width', type=positive_int,
                        help='Pixel width of the heatmap')

    parser.add_argument('vertical_average', metavar='vertical_average', type=positive_int,
                        help='Number of lines to average vertically')

    args = parser.parse_args(args)

    max_tss_file = args.max_tss_file
    region_width = args.region_width
    heatmap_width = args.heatmap_width
    vertical_average = args.vertical_average

    # Verify maxTSS file exists
    verify_bed_files(max_tss_file)

    # Make sure the regions are 1 bp long
    region_length = determine_region_length(max_tss_file)

    if region_length != 1:
        sys.stderr.write("The maxTSS bed file must have regions of 1 bp.\n")
        sys.exit(1)

    return max_tss_file, region_width, heatmap_width, vertical_average


def main(args):
    max_tss_file, region_width, heatmap_width, vertical_average = parse_args(args)

    # Expand the maxTSS file to the desired width
    expanded_regions = expand_region(max_tss_file, region_width)

    # Get the sequences from the maxTSS file
    sequences = get_sequences(expanded_regions)

    # Convert the sequences into a list of 1's and 0's
    nt_binary_lists = {
        "A": convert_sequences_to_matrix_file(sequences, "A", heatmap_width),
        "T": convert_sequences_to_matrix_file(sequences, "T", heatmap_width),
        "G": convert_sequences_to_matrix_file(sequences, "G", heatmap_width),
        "C": convert_sequences_to_matrix_file(sequences, "C", heatmap_width),
    }

    # Average vertically
    for nt in nt_binary_lists:
        matrix_filename = nt_binary_lists[nt]
        nt_binary_lists[nt] = average_matrix(nt_binary_lists[nt], vertical_average)
        remove_files(matrix_filename)

    # Make the heatmaps for each nucleotide
    output_prefix = max_tss_file.split("/")[-1].replace(".bed", "") + "_average_" + str(vertical_average) + "_"

    for nt in nt_binary_lists:
        matrix_filename = nt_binary_lists[nt]
        output_filename = output_prefix + nt + ".tiff"

        generate_heatmap(matrix_filename, "gray", output_filename, 2.2, 0, 1)
        remove_files(matrix_filename)


if __name__ == '__main__':
    main(sys.argv[1:])