import sys
import multiprocessing

from utils.make_five_and_three_dict import build_counts_dict
from utils.get_region_length import determine_region_length
from utils.verify_bed_file import verify_bed_files
from utils.output_metaplot_data import output_metaplot_data
from utils.get_metaplot_averages import get_averages
from utils.verify_region_length_is_even import verify_region_length_is_even


def get_primes_data_helper(regions_filename, three_prime_counts_dict):
    with open(regions_filename) as file:
        three_output_list = []
        three_rev_output_list = []

        for line in file:
            chromosome, left, right, gene_name, _, strand = line.split()

            three_current_list = [-1] * region_length
            three_rev_current_list = [-1] * region_length

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
                three_current_list[current_position] = three_prime_counts_dict[chromosome][strand][i]
                three_rev_current_list[current_position] = three_prime_counts_dict[chromosome][opp_strand][i]

                current_position += add_next

            three_output_list.append(three_current_list)
            three_rev_output_list.append(three_rev_current_list)

    return three_output_list, three_rev_output_list


def get_primes_data(regions_filename, sequencing_file):
    # 1. Load the three_prime_counts_dict
    _, three_prime_counts_dict = build_counts_dict(sequencing_file)

    # 2. Load 2D list containing the data to be outputted
    three_same_regions_data, three_rev_regions_data = get_primes_data_helper(regions_filename, three_prime_counts_dict)

    # 3. Get averages data
    three_same_averages = get_averages(three_same_regions_data, region_length)
    three_rev_averages = get_averages(three_rev_regions_data, region_length)

    return list(zip(three_same_averages, three_rev_averages)), sequencing_file


def parse_input(args):
    if len(args) < 2:
        print_usage()
        sys.exit(1)

    regions_filename = args[0]
    sequencing_files_list = args[1:]

    global region_length
    region_length = determine_region_length(regions_filename)

    verify_region_length_is_even(region_length)

    verify_bed_files([regions_filename] + sequencing_files_list)

    return regions_filename, sequencing_files_list


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("python3 three_prime_metaplot.py <Regions Filename> <Sequencing Files>\n")
    sys.stderr.write("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/three_prime_metaplot.rst\n")


def run_three_prime_metaplots(regions_filename, sequencing_files_list):
    pool = multiprocessing.Pool(processes=len(sequencing_files_list))
    averages = pool.starmap(get_primes_data, [(regions_filename, sequencing_file) for sequencing_file in sequencing_files_list])
    output_metaplot_data(averages, region_length, "three prime")


def main(args):
    """
    python3 three_prime_metaplot.py <Regions Filename> <Sequencing Files>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/three_prime_metaplot.rst

    :param args: arguments provided to the program
    :type args: list
    :return:
    """

    regions_filename, sequencing_files_list = parse_input(args)
    run_three_prime_metaplots(regions_filename, sequencing_files_list)

if __name__ == '__main__':
    main(sys.argv[1:])