# This will make heatmaps for A, T, G, C for a certain distance surrounding the +1 nt

import sys

from GC_bioinfo.utils.verify_bed_file import verify_bed_files
from GC_bioinfo.utils.get_region_length import determine_region_length
from GC_bioinfo.utils.run_bedtools_getfasta import run_getfasta
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.average_matrix import average_matrix
from GC_bioinfo.utils.generate_heatmap import generate_heatmap
from GC_bioinfo.utils.remove_files import remove_files


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


def convert_sequences_to_matrix_file(sequences, nucleotide):
    matrix_file = generate_random_filename(".matrix")

    with open(matrix_file, 'w') as file:
        for seq in sequences:
            curr_seq_list = []

            for nt in seq:
                if nt == nucleotide:
                    curr_seq_list.append("1")
                else:
                    curr_seq_list.append("0")

            file.write(
                "\t".join(curr_seq_list) + "\n"
            )

    return matrix_file


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo nucleotide_heatmap <maxTSS file> <region width> <vertical average>\n")
    # TODO: add the github link
    sys.stderr.write(
        "More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/multicoverage.rst\n")


def parse_args(args):
    if len(args) != 3:
        print_usage()
        sys.exit(1)

    max_tss_file, region_width, vertical_average = args

    # Verify maxTSS file exists
    verify_bed_files(max_tss_file)

    # Make sure the regions are 1 bp long
    region_length = determine_region_length(max_tss_file)

    if region_length != 1:
        sys.stderr.write("The maxTSS bed file must have regions of 1 bp.\n")
        print_usage()
        sys.exit(1)

    # Verify the region width is an integer
    try:
        region_width = int(region_width)
    except:
        sys.stderr.write("The region width must be an integer.\n")
        print_usage()
        sys.exit(1)

    # Verify the vertical average is an integer
    try:
        vertical_average = int(vertical_average)
    except:
        sys.stderr.write("The vertical average must be an integer.\n")
        print_usage()
        sys.exit(1)

    return max_tss_file, region_width, vertical_average


def main(args):
    max_tss_file, region_width, vertical_average = parse_args(args)

    # Expand the maxTSS file to the desired width
    expanded_regions = expand_region(max_tss_file, region_width)

    # Get the sequences from the maxTSS file
    sequences = get_sequences(expanded_regions)

    # Convert the sequences into a list of 1's and 0's
    nt_binary_lists = {
        "A": convert_sequences_to_matrix_file(sequences, "A"),
        "T": convert_sequences_to_matrix_file(sequences, "T"),
        "G": convert_sequences_to_matrix_file(sequences, "G"),
        "C": convert_sequences_to_matrix_file(sequences, "C"),
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