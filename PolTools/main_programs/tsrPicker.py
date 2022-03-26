import sys
import argparse

from collections import defaultdict

from PolTools.utils.build_counts_dict import build_counts_dict
from PolTools.utils.constants import hg38_chrom_sizes_file


def get_possible_tss_dict(five_prime_counts_dict, count_threshold):
    # Filter all positions that are below the threshold
    filtered_dict = defaultdict(
        lambda: {
            "+": [],
            "-": []
        }
    )

    for chromosome in five_prime_counts_dict:
        for strand in five_prime_counts_dict[chromosome]:
            for position, counts in five_prime_counts_dict[chromosome][strand].items():
                if counts >= count_threshold:
                    filtered_dict[chromosome][strand].append(
                        [chromosome, position, position + 1, "TSS", counts, strand]
                    )

            # Sort by the number of counts
            filtered_dict[chromosome][strand].sort(key=lambda x: x[-2], reverse=True)

    return filtered_dict


def find_tsrs_per_chrom(possible_tss_dict, radius):
    # Go through each chromosome and pick the largest. Go +- radius bp to make a region size of 1 + (2*radius).
    # Remove any TSSs within that that region and then pick the next max
    identified_tsrs = []
    tsr_number = 0

    for chrom in possible_tss_dict:
        for strand in possible_tss_dict[chrom]:

            tss_heights = possible_tss_dict[chrom][strand]
            while tss_heights:
                # Choose the first one since it is the max
                curr_max_tss = tss_heights[0]
                chrom, left, right, name, height, strand = curr_max_tss

                new_tsr = [chrom, left - radius, right + radius, "TSR" + str(tsr_number), height, strand]
                new_tsr_chrom, new_tsr_left, new_tsr_right, new_tsr_name, height, new_tsr_strand = new_tsr
                identified_tsrs.append(new_tsr)

                # Eliminate overlapping TSSs
                non_overlapping_tss = []
                for tss in tss_heights:
                    tss_chrom, tss_left, tss_right, tss_name, tss_score, tss_strand = tss

                    if tss_left < new_tsr_left or tss_left >= new_tsr_right:
                        # If there is no overlap between the newly defined TSR and this TSS, add it to the non overlapping
                        non_overlapping_tss.append(tss)

                tss_heights = non_overlapping_tss
                tsr_number += 1

    return identified_tsrs


def remove_invalid_tsrs(found_tsrs):
    # First load the chromosome sizes
    valid_tsrs = []

    chrom_sizes_dict = {}
    with open(hg38_chrom_sizes_file) as file:
        for line in file:
            chrom, size = line.split()
            size = int(size)
            chrom_sizes_dict[chrom] = size

    for tsr in found_tsrs:
        chrom, left, right, name, score, strand = tsr

        if chrom in chrom_sizes_dict and int(right) < chrom_sizes_dict[chrom] and int(left) >= 0:
            valid_tsrs.append(tsr)

    return valid_tsrs


def output_tsrs(found_tsrs, count_threshold, dataset):
    # Write the TSRs to a file
    output_filename = dataset.rsplit(".")[0] + "_min_" + str(count_threshold) + "-TSR.bed"

    with open(output_filename, 'w') as file:
        for tsr in found_tsrs:
            file.write(
                "\t".join(
                    [str(val) for val in tsr]
                ) + "\n"
            )


def positive_int(num):
    try:
        val = int(num)
        if val <= 0:
            raise Exception("Go to the except")
    except:
        raise argparse.ArgumentTypeError(num + " must be positive")

    return val


def parse_args(args):
    parser = argparse.ArgumentParser(prog='PolTools tsrPicker',
                                     description="Generate transcription start regions\n" +
                                                 "More information can be found at " +
                                                 "https://geoffscollins.github.io/PolTools/tsrPicker.html")

    parser.add_argument('seq_file', metavar='seq_file', type=str, help='Bed formatted sequencing file to find the TSRs')
    parser.add_argument('min_seq_depth', metavar='min_seq_depth', type=positive_int,
                        help="Minimum number of 5' ends to be considered a TSR")
    parser.add_argument('-r', '--radius', dest='radius', metavar='radius', type=positive_int, nargs='?', default=5)

    parsed_args = parser.parse_args(args)
    return parsed_args.seq_file, parsed_args.min_seq_depth, parsed_args.radius


def main(args):
    dataset, tss_strength_threshold, radius = parse_args(args)

    # Convert the thresholds to integers
    tss_strength_threshold = int(tss_strength_threshold)

    five_prime_counts_dict = build_counts_dict(dataset, 'five')

    # Get the counts for each of those positions in each of the tbp datasets
    possible_tss_dict = get_possible_tss_dict(five_prime_counts_dict, tss_strength_threshold)

    # Go through each chromosome and find tsrs
    found_tsrs = find_tsrs_per_chrom(possible_tss_dict, radius)

    # Remove tsrs that go past the genome
    final_tsrs = remove_invalid_tsrs(found_tsrs)

    output_tsrs(final_tsrs, tss_strength_threshold, dataset)


if __name__ == '__main__':
    main(sys.argv[1:])