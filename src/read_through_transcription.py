import csv
import sys
import multiprocessing

from utils.run_bedtools_coverage import run_coverage
from utils.run_bedtools_subtract import run_subtract
from utils.make_three_prime_bed_file import make_three_bed_file
from utils.make_random_filename import generate_random_filename
from utils.remove_files import remove_files
from utils.verify_bed_file import verify_bed_files
from utils.check_dependencies import check_dependencies
from utils.print_tab_delimited import print_tab_delimited


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
    three_prime_end_file = make_three_bed_file(filename)

    # Run coverage on the three bed file
    coverage_file = generate_random_filename()
    run_coverage(region_intervals_file, three_prime_end_file, output_filename=coverage_file)

    remove_files(three_prime_end_file)
    return coverage_file


def get_coverage_files(blacklisted_filenames, region_intervals_file):
    # Go through all of the blacklisted sequencing files and make 3' end files
    pool = multiprocessing.Pool(processes=len(blacklisted_filenames))
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


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo read_through_transcription.py <Regions Filename> <TSR Filename> " +
          "<Upstream Distance> <Downstream Distance> <Interval Distance> <Sequencing Files>\n")
    sys.stderr.write("\nMore information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/read_through_transcription.rst\n")


def parse_input(args):
    if len(args) == 0:
        print_usage()
        sys.exit(1)

    if len(args) < 6:
        sys.stderr.write("You did not provide all of the necessary arguments. Please try again.\n")
        print_usage()
        sys.exit(1)

    regions_filename, tsr_file, upstream_distance, downstream_distance, interval_size = args[:5]
    sequencing_files = args[5:]

    verify_bed_files(regions_filename, sequencing_files)

    try:
        upstream_distance = int(upstream_distance)

    except ValueError as _:
        # The values are not integers
        raise ValueError("The upstream distance you provided is not an integer.")

    try:
        downstream_distance = int(downstream_distance)

    except ValueError as _:
        # The values are not integers
        raise ValueError("The downstream distance you provided is not an integer.")

    try:
        interval_size = int(interval_size)

    except ValueError as _:
        # The values are not integers
        raise ValueError("The downstream distance you provided is not an integer.")


    if upstream_distance < 0 or downstream_distance < 0 or interval_size < 0:
        sys.stderr.write("Distances cannot be negative. Please enter a positive number.\n")
        sys.exit(1)

    return regions_filename, tsr_file, upstream_distance, downstream_distance, interval_size, sequencing_files


def run_read_through_transcription(regions_filename, tsr_file, upstream_distance, downstream_distance, interval_size, sequencing_files):
    # 1. Make the region intervals file from upstream distance to downstream distance in intervals
    incremented_regions_filename = make_incremented_regions(regions_filename, upstream_distance, downstream_distance, interval_size)

    # Blacklist the TSRs
    if tsr_file != 'no':
        blacklisted_filenames = blacklist_tsrs(sequencing_files, tsr_file)
        coverage_files = get_coverage_files(blacklisted_filenames, incremented_regions_filename)
        remove_files(blacklisted_filenames)
    else:
        coverage_files = get_coverage_files(sequencing_files, incremented_regions_filename)

    combined_dict = coverage_files_to_dictionary(coverage_files, sequencing_files)

    output_data(combined_dict, sequencing_files, upstream_distance, interval_size)

    # Remove all of the temporary files
    remove_files(incremented_regions_filename, coverage_files)


def main(args):
    """
    read_through_transcription <Regions Filename> <TSR Filename> <Upstream Distance> <Downstream Distance> <Interval Size>
    <Sequencing Files>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/read_through_transcription.rst

    :param args: arguments provided to the program
    :type args: list
    :return:
    """
    check_dependencies("bedtools")
    # The user must give us a bed file with the regions and a list of the sequencing files
    regions_filename, tsr_file, upstream_distance, downstream_distance, interval_size, sequencing_files = parse_input(args)

    run_read_through_transcription(regions_filename, tsr_file, upstream_distance, downstream_distance, interval_size,
                                   sequencing_files)


if __name__ == '__main__':
    main(sys.argv[1:])
