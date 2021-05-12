import sys
import argparse

from PolTools.utils.print_tab_delimited import print_tab_delimited
from PolTools.utils.verify_region_length_is_even import verify_region_length_is_even


def get_max_tsss(truQuant_filename, tsr_finder):
    max_tsss = []

    # If it is a tsrFinder file, parse it. Otherwise, it must be a truQuant file
    if tsr_finder:
        with open(truQuant_filename) as file:
            for i, line in enumerate(file):
                chrom, left, right, read_sum, five_ends, strand, tss_left, tss_right, tss_strength, avg_tss = line.split()

                # Add one to tss_left because truQuant uses genome browser positions (starts at 1)
                tss_left = str(int(tss_left) + 1)

                max_tsss.append(
                    ["TSR" + str(i), chrom, tss_left, strand, tss_strength]
                )

    with open(truQuant_filename) as file:
        for i, line in enumerate(file):
            if i != 0:
                gene_name, chromosome, pause_left, pause_right, strand, total_reads, max_tss, max_tss_five_prime_reads, *_ = line.split()
                max_tsss.append(
                    [gene_name, chromosome, max_tss, strand, max_tss_five_prime_reads]
                )

    return max_tsss


def expand_max_tss(max_tsss, region_size):
    expanded_regions = []

    for max_tss in max_tsss:
        gene_name, chromosome, position, strand, five_prime_reads = max_tss

        # We subtract one from the position because the genome browser starts at 1 and bioinformatics programs start at 0
        position = int(position) - 1

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

        # Verify the region is not negative before adding it
        if left >= 0:
            expanded_regions.append( [chromosome, left, right, gene_name, five_prime_reads, strand] )

    return expanded_regions


def output_data(expanded_regions):
    for region in expanded_regions:
        print_tab_delimited(region)


def parse_args(args):
    parser = argparse.ArgumentParser(prog='PolTools make_regions_file_centered_on_max_tss',
        description='Make a region file centered on the max TSS from truQuant\n' +
                     "More information can be found at " +
                     "https://geoffscollins.github.io/PolTools/make_regions_file_centered_on_max_tss.html")

    parser.add_argument('truQuant_file', metavar='truQuant_file', type=str,
                        help='truQuant output file ending in -truQuant_output.txt')

    parser.add_argument('-f', '--tsrFinder', action="store_true", dest='tsr_finder',
                        default=False, help='File is a tsrFinder output file.')

    parser.add_argument('region_size', metavar='region_size', type=int,
                        help='size of the region to be generated. This must be an even integer or 1')

    args = parser.parse_args(args)

    region_size = args.region_size

    if region_size != 1:
        verify_region_length_is_even(region_size)

    return args.truQuant_file, region_size, args.tsr_finder


def main(args):
    """
    PolTools make_regions_file_centered_on_max_tss <truQuant output file> <region size>
    More information can be found at https://geoffscollins.github.io/PolTools/make_regions_file_centered_on_max_tss.html

    :param args:
    :return:
    """
    truQuant_filename, region_size, tsr_finder = parse_args(args)

    max_tsss = get_max_tsss(truQuant_filename, tsr_finder)

    expanded_regions = expand_max_tss(max_tsss, region_size)
    output_data(expanded_regions)


if __name__ == '__main__':
    main(sys.argv[1:])