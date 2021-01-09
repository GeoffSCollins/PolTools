import multiprocessing
import sys

from GC_bioinfo.utils.get_metaplot_averages import average_vertically
from GC_bioinfo.utils.get_region_length import determine_region_length
from GC_bioinfo.utils.make_five_and_three_dict import build_counts_dict
from GC_bioinfo.utils.output_metaplot_data import output_metaplot_data
from GC_bioinfo.utils.verify_bed_file import verify_bed_files
from GC_bioinfo.utils.verify_region_length_is_even import verify_region_length_is_even


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


def get_primes_data(regions_filename, sequencing_file, region_length, end):
    # 1. Load the five_prime_dict
    five_prime_counts_dict, three_prime_counts_dict = build_counts_dict(sequencing_file)

    if end == "five":
        counts_dict = five_prime_counts_dict
    elif end == "three":
        counts_dict = three_prime_counts_dict
    else:
        sys.stderr.write("End provided is not five or three. Exiting ...")
        sys.exit(1)

    # 2. Load 2D list containing the data to be outputted
    sense_regions_data, divergent_regions_data = get_primes_data_helper(regions_filename, counts_dict, region_length)

    # 3. Get averages data
    same_averages = average_vertically(sense_regions_data)
    rev_averages = average_vertically(divergent_regions_data)

    return list(zip(same_averages, rev_averages)), sequencing_file


def run_metaplot(regions_filename, sequencing_files_list, region_length, end):
    pool = multiprocessing.Pool(processes=len(sequencing_files_list))
    averages = pool.starmap(get_primes_data, [(regions_filename, sequencing_file, region_length, end)
                                              for sequencing_file in sequencing_files_list])

    if end == "five":
        output_prefix = "5'"
    else:
        output_prefix = "3'"

    output_metaplot_data(averages, region_length, output_prefix)


def print_five_prime_metaplot_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo five_prime_metaplot <Regions Filename> <Sequencing Files>\n")
    sys.stderr.write("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/five_prime_metaplot.rst\n")

def print_three_prime_metaplot_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo three_prime_metaplot <Regions Filename> <Sequencing Files>\n")
    sys.stderr.write("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/three_prime_metaplot.rst\n")


def parse_input(args, end):
    if end not in ["five", "three"]:
        sys.stderr.write("End provided is not five or three. Exiting ...")
        sys.exit(1)

    if len(args) < 2:
        if end == "five":
            print_five_prime_metaplot_usage()
            sys.exit(1)
        else:
            print_three_prime_metaplot_usage()
            sys.exit(1)

    regions_filename = args[0]
    sequencing_files_list = args[1:]

    region_length = determine_region_length(regions_filename)

    verify_region_length_is_even(region_length)

    verify_bed_files([regions_filename] + sequencing_files_list)

    return regions_filename, sequencing_files_list, region_length
