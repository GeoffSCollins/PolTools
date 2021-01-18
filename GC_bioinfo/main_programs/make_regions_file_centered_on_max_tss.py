import sys
import argparse

from GC_bioinfo.utils.print_tab_delimited import print_tab_delimited
from GC_bioinfo.utils.verify_region_length_is_even import verify_region_length_is_even


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

        if region_size == 1:
            # Only the maxTSS
            left = position
            right = position + 1
        else:
            left = position - int(region_size / 2)
            right = position + int(region_size / 2)

            if strand == "-":
                # We shift the right by 1
                left = int(left) + 1
                right = int(right) + 1

        expanded_regions.append( [chromosome, left, right, gene_name, five_prime_reads, strand] )

    return expanded_regions


def output_data(expanded_regions):
    for region in expanded_regions:
        print_tab_delimited(region)


def parse_args(args):
    parser = argparse.ArgumentParser(prog='GC_bioinfo make_regions_file_centered_on_max_tss',
        description='Make a region file centered on the max TSS from truQuant\n' +
                     "More information can be found at " +
                     "https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/make_regions_file_centered_on_max_tss.rst")

    parser.add_argument('truQuant_file', metavar='truQuant_file', type=str,
                        help='truQuant output file ending in -truQuant_output.txt')

    parser.add_argument('region_size', metavar='region_size', type=int,
                        help='size of the region to be generated. This must be an even integer or 1')

    args = parser.parse_args(args)

    truQuant_filename = args.truQuant_file
    region_size = args.region_size

    if region_size != 1:
        verify_region_length_is_even(region_size)

    return truQuant_filename, region_size


def main(args):
    """
    GC_bioinfo make_regions_file_centered_on_max_tss <truQuant output file> <region size>
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