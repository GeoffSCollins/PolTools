import csv
import multiprocessing
import sys
import argparse

from PolTools.utils.check_dependencies import check_dependencies
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.make_read_end_file import make_read_end_file
from PolTools.utils.print_tab_delimited import print_tab_delimited
from PolTools.utils.remove_files import remove_files
from PolTools.utils.bedtools_utils.run_bedtools_coverage import run_coverage
from PolTools.utils.bedtools_utils.run_bedtools_subtract import run_subtract


def make_incremented_regions(regions_filename, upstream_distance, downstream_distance, interval_size):
    # Using the regions provided, make incremented regions
    with open(regions_filename) as file:
        regions = []
        for line in file:
            chromosome, left, right, gene_name, score, strand = line.split()

            if strand == "+":
                # We use the right position because that is the TES
                region_left = int(right) - upstream_distance
                region_right = int(right) + downstream_distance
            else:
                # We use the left position because that is the TES
                region_left = int(left) - downstream_distance
                region_right = int(left) + upstream_distance

            # Add the region to regions
            regions.append([chromosome, region_left, region_right, gene_name, score, strand])

    # Go through all the regions and make the incremented ones
    incremented_regions = []
    for region in regions:
        chromosome, left, right, gene_name, score, strand = region

        if strand == "+":
            # We work from left to right
            for i in range(left + interval_size, right + 1, interval_size):
                # Looping through each interval region
                incremented_regions.append([chromosome, i - interval_size, i, gene_name, score, strand])

        else:
            # We work from right to left
            if left + (interval_size - 1) > 0:
                for i in range(right, left + (interval_size - 1), (-1* interval_size)):
                    # Looping through each interval region
                    incremented_regions.append([chromosome, i - interval_size, i, gene_name, score, strand])

    region_intervals_filename = generate_random_filename()

    with open(region_intervals_filename, 'w') as tmp_region_file:
        output_writer = csv.writer(tmp_region_file, delimiter='\t', lineterminator='\n')
        for region in incremented_regions:
            output_writer.writerow(region)

    return region_intervals_filename


def blacklist_tsrs(sequencing_files, tsr_file):
    # Go through all of the sequencing files and blacklist TSRs
    blacklisted_filenames = []
    jobs = []
    for seq_filename in sequencing_files:
        curr_filename = generate_random_filename()
        blacklisted_filenames.append(curr_filename)
        p = multiprocessing.Process(target=run_subtract, args=[seq_filename, tsr_file],
                             kwargs={'output_filename': curr_filename})
        jobs.append(p)
        p.start()

    for job in jobs:
        job.join()

    return blacklisted_filenames


def get_coverage_files_helper(filename, region_intervals_file):
    # First makes the three bed file
    three_prime_end_file = make_read_end_file(filename, 'three')

    # Run coverage on the three bed file
    coverage_file = generate_random_filename()
    run_coverage(region_intervals_file, three_prime_end_file, output_filename=coverage_file)

    remove_files(three_prime_end_file)
    return coverage_file


def get_coverage_files(blacklisted_filenames, region_intervals_file, max_threads):
    # Go through all of the blacklisted sequencing files and make 3' end files
    with multiprocessing.Pool(processes=max_threads) as pool:
        coverage_files = pool.starmap(get_coverage_files_helper, [(filename, region_intervals_file) for filename in blacklisted_filenames])

    return coverage_files


def coverage_files_to_dictionary(coverage_files, sequencing_files):
    # Has keys of the region number and values of a list containing the counts for each file
    combined_dict = {}

    for i, filename in enumerate(coverage_files):
        with open(filename) as file:
            # We have open the bedtools coverage results file
            region_number = 0
            previous_gene_name = ""

            for j, line in enumerate(file):
                chromosome, left, right, gene_name, score, strand, counts, *_ = line.split()

                if previous_gene_name != gene_name:
                    region_number = 0

                if region_number not in combined_dict:
                    combined_dict[region_number] = [0] * len(sequencing_files)

                combined_dict[region_number][i] += int(counts)

                region_number += 1
                previous_gene_name = gene_name

    return combined_dict


def output_data(combined_dict, sequencing_files, upstream_distance, interval_size):
    # Now the data is in the combined_dict, we need to reduce it back down to positions again

    # First print out the headers
    print_tab_delimited(["Position"] + [seq_file.split("/")[-1] for seq_file in sequencing_files])

    real_position = upstream_distance * -1
    # Now output all of the data
    for position in combined_dict:
        print_tab_delimited([real_position] + combined_dict[position])
        real_position += interval_size


def parse_input(args):

    def positive_int(value):
        try:
            val = int(value)
            if val < 0:
                raise argparse.ArgumentTypeError(value + " must be a postive integer")

        except:
            raise argparse.ArgumentTypeError(value + " must be a postive integer")

        return val

    parser = argparse.ArgumentParser(prog='PolTools read_through_transcription',
                                     description='Quantify three prime ends around the TES of genes\n' +
                                     "More information can be found at " +
                                     "https://geoffscollins.github.io/PolTools/read_through_transcription.html")

    parser.add_argument('regions_filename', metavar='regions_filename', type=str,
                        help='Bed formatted regions file with an even region length or a region length of one.')

    parser.add_argument('tsr_file', metavar='tsr_file', type=str,
                        help='tsrFinder output file')

    parser.add_argument('upstream_distance', metavar='upstream_distance', type=positive_int,
                        help='Distance upstream of the TES')

    parser.add_argument('downstream_distance', metavar='downstream_distance', type=positive_int,
                        help='Distance downstream of the TES')

    parser.add_argument('interval_size', metavar='interval_size', type=positive_int,
                        help='Size of the intervals to split the regions to quantify into')

    parser.add_argument('seq_files', metavar='sequencing_files', nargs='+', type=str,
                        help='Bed formatted files from the sequencing experiment')

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    args = parser.parse_args(args)

    regions_filename = args.regions_filename
    tsr_file = args.tsr_file
    upstream_distance = args.upstream_distance
    downstream_distance = args.downstream_distance
    interval_size = args.interval_size
    sequencing_files = args.seq_files
    max_threads = args.threads

    return regions_filename, tsr_file, upstream_distance, downstream_distance, interval_size, sequencing_files, max_threads


def run_read_through_transcription(regions_filename, tsr_file, upstream_distance, downstream_distance, interval_size, sequencing_files, max_threads):
    # 1. Make the region intervals file from upstream distance to downstream distance in intervals
    incremented_regions_filename = make_incremented_regions(regions_filename, upstream_distance, downstream_distance, interval_size)

    # Blacklist the TSRs
    if tsr_file != 'no':
        blacklisted_filenames = blacklist_tsrs(sequencing_files, tsr_file)
        coverage_files = get_coverage_files(blacklisted_filenames, incremented_regions_filename, max_threads)
        remove_files(blacklisted_filenames)
    else:
        coverage_files = get_coverage_files(sequencing_files, incremented_regions_filename, max_threads)

    combined_dict = coverage_files_to_dictionary(coverage_files, sequencing_files)

    output_data(combined_dict, sequencing_files, upstream_distance, interval_size)

    # Remove all of the temporary files
    remove_files(incremented_regions_filename, coverage_files)


def main(args):
    """
    read_through_transcription <Regions Filename> <TSR Filename> <Upstream Distance> <Downstream Distance> <Interval Size>
    <Sequencing Files>
    More information can be found at https://geoffscollins.github.io/PolTools/read_through_transcription.html

    :param args: arguments provided to the program
    :type args: list
    :return:
    """
    check_dependencies("bedtools")
    # The user must give us a bed file with the regions and a list of the sequencing files
    regions_filename, tsr_file, upstream_distance, downstream_distance, interval_size, sequencing_files, max_threads = parse_input(args)

    run_read_through_transcription(regions_filename, tsr_file, upstream_distance, downstream_distance, interval_size,
                                   sequencing_files, max_threads)


if __name__ == '__main__':
    main(sys.argv[1:])
