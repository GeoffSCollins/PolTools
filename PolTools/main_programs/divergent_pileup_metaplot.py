import multiprocessing
import os
import sys
import argparse

from PolTools.utils.check_dependencies import check_dependencies
from PolTools.utils.get_region_length import determine_region_length
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.output_metaplot_data import output_metaplot_data
from PolTools.utils.remove_files import remove_files
from utils.bedtools_utils.run_bedtools_coverage import run_coverage
from PolTools.utils.verify_region_length_is_even import verify_region_length_is_even


def split_bed_file(bed_file):
    fw_filename = generate_random_filename()
    rv_filename = generate_random_filename()

    with open(fw_filename, 'w') as fw_file:
        with open(rv_filename, 'w') as rv_file:
            with open(bed_file) as file:
                for line in file:
                    if "+" in line:
                        fw_file.write(line)
                    else:
                        rv_file.write(line)

    return fw_filename, rv_filename


def get_pileups_helper(region_files, seq_files, region_length, opposite_strand):
    fw_region_file, rv_region_file = region_files
    fw_seq_file, rv_seq_file = seq_files

    # Run bedtools coverage on all of them
    fw_output = run_coverage(fw_region_file, fw_seq_file, flags=["-d"])
    rv_output = run_coverage(rv_region_file, rv_seq_file, flags=["-d"])

    # Combine the two files together
    combined_file = generate_random_filename()
    os.system("cat " + fw_output + " " + rv_output + " > " + combined_file)

    counts_list = [0] * region_length

    with open(combined_file) as file:
        for line in file:
            if not opposite_strand:
                if "+" in line:
                    position = int(line.split()[-2]) - 1  # Subtract 1 because position starts at 1
                else:
                    position = region_length - int(line.split()[-2])

                counts_list[position] += int(line.split()[-1])
            else:
                # If it is the opposite strand, we reverse the position and make the counts negative
                if "-" in line:
                    position = int(line.split()[-2]) - 1  # Subtract 1 because position starts at 1
                else:
                    position = region_length - int(line.split()[-2])

                counts_list[position] -= int(line.split()[-1])

    remove_files(fw_output, rv_output, combined_file)

    return counts_list


def make_rev_region_file(region_filename):
    rev_region_filename = generate_random_filename()

    with open(region_filename) as file:
        with open(rev_region_filename, 'w') as outfile:
            for line in file:
                if "+" in line:
                    outfile.write(line.replace("+", "-"))
                else:
                    outfile.write(line.replace("-", "+"))

    return rev_region_filename


def get_pileups(region_files, rev_region_files, seq_files, region_length, orig_seq_filename):
    pileups = get_pileups_helper(region_files, seq_files, region_length, False)
    rev_pileups = get_pileups_helper(rev_region_files, seq_files, region_length, True)

    fw_region_filename, rv_region_filename = region_files

    with open(fw_region_filename) as file:
        number_of_regions_in_fw = sum(1 for _ in file)

    with open(rv_region_filename) as file:
        number_of_regions_in_rv = sum(1 for _ in file)

    number_of_regions = number_of_regions_in_fw + number_of_regions_in_rv

    same_averages = [x / number_of_regions for x in pileups]
    rev_averages = [x / number_of_regions for x in rev_pileups]

    return list(zip(same_averages, rev_averages)), orig_seq_filename


def parse_args(args):
    def positive_int(num):
        try:
            val = int(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")
        return val

    parser = argparse.ArgumentParser(prog='PolTools divergent_pileup_metaplot',
        description='Generate a metaplot of the pileup reads.\n'+
        'More information can be found at https://geoffscollins.github.io/PolTools/divergent_pileup_metaplot.html')

    parser.add_argument('regions_file', metavar='regions_file', type=str,
                        help='Bed formatted file containing all the regions you want to average the sequences')

    parser.add_argument('seq_files', metavar='sequencing_files', nargs='+', type=str,
                        help='Bed formatted files from the sequencing experiment')

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    args = parser.parse_args(args)
    regions_filename = args.regions_file
    max_threads = args.threads

    region_length = determine_region_length(regions_filename)
    verify_region_length_is_even(region_length)

    sequencing_files_list = args.seq_files

    return regions_filename, sequencing_files_list, region_length, max_threads


def run_divergent_pileup_plots(regions_filename, sequencing_files_list, region_length, max_threads):
    # Make the opposite stranded regions file
    rev_region_filename = make_rev_region_file(regions_filename)

    # Split the files by strand
    region_files = split_bed_file(regions_filename)
    rev_region_files = split_bed_file(rev_region_filename)

    split_seq_filenames = [(split_bed_file(filename), filename) for filename in sequencing_files_list]

    with multiprocessing.Pool(processes=max_threads) as pool:
        averages = pool.starmap(get_pileups, [(region_files, rev_region_files, seq_files, region_length, filename)
                                          for (seq_files, filename) in split_seq_filenames])

    split_seq_filenames_to_delete = [tup[0] for tup in split_seq_filenames]

    # Remove all the stranded region files
    remove_files(region_files, rev_region_files, rev_region_filename, split_seq_filenames_to_delete)

    output_metaplot_data(averages, region_length, "")


def main(args):
    """
    python3 divergent_pileup_metaplot <Regions Filename> <Sequencing Files>
    More information can be found at https://geoffscollins.github.io/PolTools/divergent_pileup_metaplot.html

    :param args: arguments provided to the program
    :type args: list
    :return:
    """
    check_dependencies("bedtools")
    regions_filename, sequencing_files_list, region_length, max_threads = parse_args(args)
    run_divergent_pileup_plots(regions_filename, sequencing_files_list, region_length, max_threads)


if __name__ == '__main__':
    main(sys.argv[1:])
