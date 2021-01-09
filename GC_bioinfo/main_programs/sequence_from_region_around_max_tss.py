import sys

# Input is a maxTSS file and the position
# Output is the maxTSS file with the sequence at that position


from GC_bioinfo.utils.get_region_length import determine_region_length
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.run_bedtools_getfasta import run_getfasta
from GC_bioinfo.utils.remove_files import remove_files

from collections import defaultdict


def get_regions_file(max_tss_file, search):
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

                outfile.write(
                    "\t".join([chromosome, region_left, region_right, name, score, strand]) + "\n"
                )

                gene_names.append(name)

    return regions_filename, gene_names


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


def output_sequence_dict(sequence_dict):
    for gene_name in sequence_dict:
        print(">" + gene_name)
        print(sequence_dict[gene_name])


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo sequence_from_region_around_max_tss <maxTSS file> <region>\n")
    sys.stderr.write(
        "More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/sequence_from_region_around_max_tss.rst\n")


def parse_args(args):
    if len(args) != 2:
        print_usage()
        sys.exit(1)

    max_tss_file, search_region = args

    region_length = determine_region_length(max_tss_file)

    if region_length != 1:
        sys.stderr.write("The maxTSS bed file must have regions of 1 bp.\n")
        print_usage()
        sys.exit(1)

    # Convert search region syntax to a list
    search_left = [
        search_region[0],
        int(search_region.split("_")[0][1:])
    ]

    search_right = [
        search_region.split("_")[1][0],
        int(search_region.split("_")[1][1:])
    ]

    search = [search_left, search_right]

    return max_tss_file, search


def main(args):
    max_tss_file, search = parse_args(args)

    # Make a bed file which has the regions of interest
    regions_file, gene_names = get_regions_file(max_tss_file, search)

    # Get the sequences for each region
    fasta_file = run_getfasta(regions_file)

    sequence_dict = map_sequences_to_gene_names(gene_names, fasta_file)

    output_sequence_dict(sequence_dict)

    remove_files(regions_file, fasta_file)


if __name__ == '__main__':
    main(sys.argv[1:])