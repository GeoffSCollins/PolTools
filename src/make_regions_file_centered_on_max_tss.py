
import sys

from utils.print_tab_delimited import print_tab_delimited
from utils.verify_region_length_is_even import verify_region_length_is_even

def get_max_tsss(truQuant_filename):
    max_tsss = []

    with open(truQuant_filename) as file:
        for i, line in enumerate(file):
            if i != 0:
                gene_name, chromosome, pause_left, pause_right, strand, total_reads, max_tss, max_tss_five_prime_reads, *_ = line.split()
                max_tsss.append([gene_name, chromosome, max_tss, strand, max_tss_five_prime_reads])

    return max_tsss


def expand_max_tss(max_tsss, region_size):
    expanded_regions = []

    for max_tss in max_tsss:
        gene_name, chromosome, position, strand, five_prime_reads = max_tss

        position = int(position)

        left = position - int(region_size / 2)
        right = position + int(region_size / 2)

        if strand == "-":
            # If the strand is negative, we shift the right 1
            left += 1
            right += 1

        expanded_regions.append( [chromosome, left, right, gene_name, five_prime_reads, strand] )

    return expanded_regions


def output_data(expanded_regions):
    for region in expanded_regions:
        print_tab_delimited(region)


def parse_args(args):
    if len(args) != 2:
        print_usage()
        sys.exit(1)

    truQuant_filename, region_size = args

    region_size = int(region_size)
    verify_region_length_is_even(region_size)

    return truQuant_filename, region_size


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("python3 make_regions_file_centered_on_max_tss.py <truQuant output file> <region size>\n")
    sys.stderr.write("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/make_regions_file_centered_on_max_tss.rst\n")


def main(args):
    """
    python3 make_regions_file_centered_on_max_tss.py <truQuant output file> <region size>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/make_regions_file_centered_on_max_tss.rst

    :param args:
    :return:
    """
    truQuant_filename, region_size = parse_args(args)

    max_tsss = get_max_tsss(truQuant_filename)
    expanded_regions = expand_max_tss(max_tsss, region_size)
    output_data(expanded_regions)


if __name__ == '__main__':
    main(sys.argv[1:])