import multiprocessing
import sys
import argparse
from itertools import chain

from PolTools.utils.print_tab_delimited import print_tab_delimited
from PolTools.utils.get_region_length import determine_region_length
from PolTools.utils.build_counts_dict import build_counts_dict
from PolTools.utils.verify_bed_file import verify_bed_files
from PolTools.utils.verify_region_length_is_even import verify_region_length_is_even

def average_vertically(input_2d_list):
    """
    Takes in a 2d list and gets the average of each x position

    :param input_2d_list:
    :return: a list of the averaged values
    :rtype: list
    """
    height = len(input_2d_list)

    if height == 0:
        return []

    width = len(input_2d_list[0])

    averages_list = [0] * width

    # We loop through each base in the region
    for curr_position in range(width):
        current_sum = 0

        # Loop through all of the regions
        for region in input_2d_list:
            current_sum += region[curr_position]

        avg = current_sum / height
        averages_list[curr_position] = avg

    return averages_list


def get_primes_data_helper(regions_filename, counts_dict, region_length):
    with open(regions_filename) as file:
        output_list = []
        rev_output_list = []

        for line in file:
            chromosome, left, right, gene_name, _, strand = line.split()

            curr_list = [-1] * region_length
            rev_current_list = [-1] * region_length

            if strand == "+":
                current_position = 0
                add_next = 1
                opp_strand = "-"
            else:
                current_position = region_length - 1
                add_next = -1
                opp_strand = "+"

            # Now we loop through the region left to right
            for i in range(int(left), int(right)):
                curr_list[current_position] = counts_dict[chromosome][strand][i]
                rev_current_list[current_position] = -1 * counts_dict[chromosome][opp_strand][i]

                current_position += add_next

            output_list.append(curr_list)
            rev_output_list.append(rev_current_list)

    return output_list, rev_output_list


def get_primes_data(regions_filename, sequencing_file, region_length, read_type):
    counts_dict = build_counts_dict(sequencing_file, read_type)

    # 2. Load 2D list containing the data to be outputted
    sense_regions_data, divergent_regions_data = get_primes_data_helper(regions_filename, counts_dict, region_length)

    # 3. Get averages data
    same_averages = average_vertically(sense_regions_data)
    rev_averages = average_vertically(divergent_regions_data)

    return list(zip(same_averages, rev_averages)), sequencing_file


def output_metaplot_data(averages, region_length, prime_name):
    """

    :param averages: averages list from the metaplots programs
    :type averages: list
    :param region_length: length of the region
    :type region_length: int
    :param prime_name: either "five prime" or "three prime"
    :type prime_name: str
    :return:
    """
    avgs_data, files = [x for x in zip(*averages)]

    # Merge all of the lists together
    merged_list = [list(chain.from_iterable(x)) for x in zip(*avgs_data)]

    # 5. Put the data into a file
    header = ["Position"]
    # Write the header first
    for file in files:
        if prime_name:
            # Include a space before print the prime name
            header.append(file.split("/")[-1] + " " + prime_name + " sense strand")
            header.append(file.split("/")[-1] + " " + prime_name + " divergent strand")
        else:
            header.append(file.split("/")[-1] + " sense strand")
            header.append(file.split("/")[-1] + " divergent strand")

    print_tab_delimited(header)

    for i, base_list in enumerate(merged_list):
        position = i - (region_length / 2)
        if position >= 0:
            position += 1

        print_tab_delimited([position] + base_list)


def run_metaplot(regions_filename, sequencing_files_list, region_length, read_type, max_threads):
    with multiprocessing.Pool(processes=max_threads) as pool:
        averages = pool.starmap(get_primes_data, [(regions_filename, sequencing_file, region_length, read_type)
                                                  for sequencing_file in sequencing_files_list])

    if read_type == "five":
        output_prefix = "5'"
    elif read_type == "three":
        output_prefix = "3'"
    else:
        output_prefix = "whole"

    output_metaplot_data(averages, region_length, output_prefix)


def parse_input(args):
    def positive_int(num):
        try:
            val = int(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")

        return val

    parser = argparse.ArgumentParser(prog="PolTools metaplot",
                                     description="Compute the number of 5'/3'/whole prime reads at each position of the given region.\n" +
                                                 "More information can be found at " +
                                                 "https://github.com/GeoffSCollins/PolTools/blob/master/docs/metaplot.rst")

    parser.add_argument('read_type', metavar='read type', type=str, choices=["five", "three", "whole"],
                        help='either five, three, or whole')

    parser.add_argument('regions_file', metavar='regions_file', type=str,
                        help='Bed formatted file containing all the regions you want to average the sequences')

    parser.add_argument('seq_files', metavar='sequencing_files', nargs='+', type=str,
                        help='Bed formatted files from the sequencing experiment')

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    args = parser.parse_args(args)
    regions_filename = args.regions_file
    sequencing_files_list = args.seq_files
    max_threads = args.threads

    verify_bed_files(regions_filename)
    region_length = determine_region_length(regions_filename)
    verify_region_length_is_even(region_length)

    return args.read_type, regions_filename, sequencing_files_list, region_length, max_threads


def main(args):
    read_type, regions_filename, sequencing_files_list, region_length, max_threads = parse_input(args)
    run_metaplot(regions_filename, sequencing_files_list, region_length, read_type, max_threads)


if __name__ == '__main__':
    main(sys.argv[1:])
