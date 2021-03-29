import sys
import argparse

from PolTools.utils.get_region_length import determine_region_length
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.bedtools_utils.run_bedtools_getfasta import run_getfasta
from PolTools.utils.remove_files import remove_files
from PolTools.utils.constants import hg38_chrom_sizes_random_file

from collections import defaultdict


def get_regions_file(max_tss_file, search, chrom_sizes):
    search_left, search_right = search

    if search_left[0] == "+":
        distance_to_move_upstream = search_left[1] - 1
    else:
        distance_to_move_upstream = -1 * search_left[1]

    if search_right[0] == "+":
        distance_to_move_downstream = search_right[1] - 1
    else:
        distance_to_move_downstream = -1 * search_right[1]

    regions_filename = generate_random_filename()

    gene_names = []

    with open(max_tss_file) as file:
        with open(regions_filename, 'w') as outfile:
            for line in file:
                chromosome, left, right, name, score, strand = line.split()

                left = int(left)
                right = int(right)

                # Move the region to the search
                if strand == "+":
                    region_left = str(left + distance_to_move_upstream)
                    region_right = str(right + distance_to_move_downstream)
                else:
                    region_left = str(left - distance_to_move_downstream)
                    region_right = str(right - distance_to_move_upstream)

                # Make sure the region is possible
                if int(region_left) < 0 or int(region_right) > chrom_sizes[chromosome]:
                    continue

                outfile.write(
                    "\t".join([chromosome, region_left, region_right, name, score, strand]) + "\n"
                )

                gene_names.append(name)

    return regions_filename, gene_names


def get_chrom_sizes():
    chrom_sizes = {}

    with open(hg38_chrom_sizes_random_file) as file:
        for line in file:
            chrom, size = line.split()
            chrom_sizes[chrom] = int(size)

    return chrom_sizes


def map_sequences_to_gene_names(gene_names, fasta_file):
    sequence_dict = defaultdict(str)

    with open(fasta_file) as file:
        for i, line in enumerate(file):
            if i == 0:
                # This line has the >
                pass
            else:
                # This line has the sequence
                sequence = line.rstrip().upper()

                curr_gene_name = gene_names[int(i / 2)]

                sequence_dict[curr_gene_name] = sequence

    return sequence_dict


def output_sequence_dict(sequence_dict, underscore):
    end = "_" if underscore else "\n"
    for gene_name in sequence_dict:
        print(">" + gene_name, end=end)
        print(sequence_dict[gene_name])


def parse_args(args):
    # Argparse does not allow the input numbers to start with a -. So, we add a space to the input so it works.
    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    parser = argparse.ArgumentParser(prog='PolTools sequence_from_region_around_max_tss',
                                     description='Get the sequence from a certain region around the max TSS\n' +
                                                 "More information can be found at " +
                                                 "https://geoffscollins.github.io/PolTools/sequence_from_region_around_max_tss.html"
                                     )

    parser.add_argument('max_tss_file', metavar='max_tss_file', type=str,
                        help='Bed formatted file which has the base of the max TSS')

    parser.add_argument('left', metavar='left', type=int,
                        help='Left end of the region to get the sequence.')

    parser.add_argument('right', metavar='right', type=int,
                        help='Left end of the region to get the sequence.')

    parser.add_argument('-u', '--underscore', action="store_true", dest='underscore',
                        default=False, help='Split the sequences by an underscore instead of a new line.')

    args = parser.parse_args(args)
    max_tss_file = args.max_tss_file
    left = args.left
    right = args.right

    if 0 in [left, right]:
        sys.stderr.write("There is no 0 nucleotide. Please try again.\n")
        sys.exit(1)

    region_length = determine_region_length(max_tss_file)

    if region_length != 1:
        sys.stderr.write("The maxTSS bed file must have regions of 1 bp.\n")
        sys.exit(1)

    # Convert search region syntax to a list
    search_left = [
        "+" if left > 0 else "-",
        abs(left)
    ]

    search_right = [
        "+" if right > 0 else "-",
        abs(right)
    ]

    search = [search_left, search_right]

    return args.underscore, max_tss_file, search


def main(args):
    underscore, max_tss_file, search = parse_args(args)

    chrom_sizes = get_chrom_sizes()

    # Make a bed file which has the regions of interest
    regions_file, gene_names = get_regions_file(max_tss_file, search, chrom_sizes)

    # Get the sequences for each region
    fasta_file = run_getfasta(regions_file)

    sequence_dict = map_sequences_to_gene_names(gene_names, fasta_file)

    output_sequence_dict(sequence_dict, underscore)

    remove_files(regions_file, fasta_file)


if __name__ == '__main__':
    main(sys.argv[1:])
