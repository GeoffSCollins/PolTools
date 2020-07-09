import sys

from itertools import chain

from utils.make_five_and_three_dict import build_counts_dict
from utils.get_region_length import determine_region_length
from utils.verify_bed_file import verify_bed_files
from utils.print_tab_delimited import print_tab_delimited


def get_primes_data(regions_filename, five_prime_counts_dict):
    with open(regions_filename) as file:
        five_output_list = []
        five_rev_output_list = []

        for line in file:
            chromosome, left, right, gene_name, _, strand = line.split()

            five_current_list = [0] * region_length
            five_rev_current_list = [0] * region_length

            if strand == "+":
                current_position = -1 * region_length
                add_next = 1
                opp_strand = "-"
            else:
                current_position = region_length
                add_next = -1
                opp_strand = "+"

            # Now we loop through the region left to right
            for i in range(int(left), int(right) + 1):  # right+1 to include right
                # If the position is in the dictionary, add it to the 2D list

                if i in five_prime_counts_dict[chromosome][strand]:
                    five_current_list[current_position] = five_prime_counts_dict[chromosome][strand][i]

                if i in five_prime_counts_dict[chromosome][opp_strand]:
                    five_rev_current_list[current_position] = five_prime_counts_dict[chromosome][opp_strand][i]

                current_position += add_next

            five_output_list.append(five_current_list)
            five_rev_output_list.append(five_rev_current_list)

    return five_output_list, five_rev_output_list


def get_averages(input_2d_list):
    # This will take in a 2d list containing the regions and counts at that specific base
    # and output a list of the averages at each base
    averages_list = [0] * region_length

    # We loop through each base in the region
    for current_base_position, counts in enumerate(input_2d_list[0]):
        current_sum = 0

        # The total is the number of regions
        total = len(input_2d_list)

        # Loop through all of the regions
        for region in input_2d_list:
            current_sum += region[current_base_position]

        avg = current_sum / total
        averages_list[current_base_position] = avg

    return averages_list


def output_data(files, merged_list):
    header = ["Position"]
    # Write the header first
    for file in files:
        header.append(file + " 5' same strand")
        header.append(file + " 5' reverse strand")

    print_tab_delimited(header)

    for i, base_list in enumerate(merged_list):
        position = i - region_length / 2

        if position >= 0:
            position += 1

        print_tab_delimited([position] + base_list)



def parse_input(args):
    if len(args) < 2:
        print("You did not provide any sequencing files! Exiting ...")
        sys.exit(1)

    regions_filename = args[0]
    sequencing_files_list = args[1:]

    global region_length
    region_length = determine_region_length(regions_filename)

    verify_bed_files([regions_filename] + sequencing_files_list)

    return regions_filename, sequencing_files_list


def print_usage():
    print("Usage: ")
    print("python3 five_prime_metaplot.py <Regions Filename> <Sequencing Files>")
    print("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/five_prime_metaplot.rst")


def run_five_prime_metaplot(regions_filename, sequencing_files_list):
    # Output dictionary has keys of sequencing files and values of a list containing average lists
    # The order of the list goes as follows: 5'same, 5'rev
    output_dictionary = {}

    for sequencing_file in sequencing_files_list:
        # 1. Load the five_prime_dict
        five_prime_counts_dict, _ = build_counts_dict(sequencing_file)

        # 2. Load 2D list containing the data to be outputted
        five_same_regions_data, five_rev_regions_data = get_primes_data(regions_filename, five_prime_counts_dict)

        # 3. Get averages data
        five_same_averages = get_averages(five_same_regions_data)
        five_rev_averages = get_averages(five_rev_regions_data)

        # The transposed_l looks like [ [], [] ] where the brackets contain an individual nt
        output_dictionary[sequencing_file] = list(zip(five_same_averages, five_rev_averages))

    files = []
    all_data = []
    for file in output_dictionary:
        all_data.append(output_dictionary[file])
        files.append(file.split("/")[-1])

    # Merge all of the lists together
    merged_list = [list(chain.from_iterable(x)) for x in zip(*all_data)]
    output_data(files, merged_list)

def main(args):
    """
    python3 five_prime_metaplot.py <Regions Filename> <Sequencing Files>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/five_prime_metaplot.rst

    :param args: arguments provided to the program
    :type args: list
    :return:
    """
    regions_filename, sequencing_files_list = parse_input(args)
    run_five_prime_metaplot(regions_filename, sequencing_files_list)


if __name__ == '__main__':
    main(sys.argv[1:])