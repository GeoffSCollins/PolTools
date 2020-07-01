'''
This python script takes in two files: a regions file (in bed format) and a sequencing file.
It will output a single file which has the number of 5' ends at a distance away from the annotated start site for each strand
'''

# Todo: Test this program


import sys
import csv

from itertools import chain

from utils import make_five_and_three_dict, get_region_length


def get_primes_data(regions_filename, region_length):
    with open(regions_filename) as file:
        five_output_list = []
        five_rev_output_list = []

        three_output_list = []
        three_rev_output_list = []

        for line in file:
            chromosome, left, right, gene_name, _, strand = line.split()

            five_current_list = [0] * region_length
            five_rev_current_list = [0] * region_length

            three_current_list = [0] * region_length
            three_rev_current_list = [0] * region_length

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

                if i in three_prime_counts_dict[chromosome][strand]:
                    three_current_list[current_position] = three_prime_counts_dict[chromosome][strand][i]

                if i in three_prime_counts_dict[chromosome][opp_strand]:
                    three_rev_current_list[current_position] = three_prime_counts_dict[chromosome][opp_strand][i]

                current_position += add_next

            five_output_list.append(five_current_list)
            five_rev_output_list.append(five_rev_current_list)

            three_output_list.append(three_current_list)
            three_rev_output_list.append(three_rev_current_list)

    return five_output_list, five_rev_output_list, three_output_list, three_rev_output_list


def get_averages(input_2d_list):
    # This will take in a 2d list containing the regions and counts at that specific base
    # and output a list of the averages at each base
    averages_list = [0] * 1001

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


if __name__ == '__main__':
    regions_filename = sys.argv[1]
    sequencing_files_list = sys.argv[2:]

    output_filename = regions_filename.replace(".bed", "_metaplots.tsv")

    if len(sequencing_files_list) == 0:
        print("You did not provide any sequencing files! Exiting ...")
        sys.exit(1)

    region_length = get_region_length.determine_region_length(regions_filename)

    # Output dictionary has keys of sequencing files and values of a list containing average lists
    # The order of the list goes as follows: 5'same, 5'rev, 3'same, 3'rev, pileup
    output_dictionary = {}

    for sequencing_file in sequencing_files_list:
        # 1. Load the five_prime_dict
        five_prime_counts_dict, three_prime_counts_dict = make_five_and_three_dict.build_counts_dict(sequencing_file)

        # 2. Load 2D list containing the data to be outputted
        five_same_regions_data, five_rev_regions_data, three_same_regions_data, three_rev_regions_data = \
            get_primes_data(regions_filename, region_length)

        # 3. Get averages data
        five_same_averages = get_averages(five_same_regions_data)
        five_rev_averages = get_averages(five_rev_regions_data)
        three_same_averages = get_averages(three_same_regions_data)
        three_rev_averages = get_averages(three_rev_regions_data)

        l = [five_same_averages, five_rev_averages, three_same_averages, three_rev_averages]
        transposed_l = list(map(list, zip(*l)))

        # The transposed_l looks like [ [], [], [], [] ] where the brackets contain an individual nt
        output_dictionary[sequencing_file] = transposed_l

    files = []
    all_data = []
    for file in output_dictionary:
        all_data.append(output_dictionary[file])
        files.append(file)

    # Merge all of the lists together
    merged_list = [list(chain.from_iterable(x)) for x in zip(*all_data)]

    # 5. Put the data into a file
    with open(output_filename, 'w') as file:
        output_writer = csv.writer(file, delimiter='\t', lineterminator='\n')

        header = ["Position"]
        # Write the header first
        for file in files:
            header.append(file + " 5' same strand")
            header.append(file + " 5' reverse strand")
            header.append(file + " 3' same strand")
            header.append(file + " 3' reverse strand")

        output_writer.writerow(header)

        for i, base_list in enumerate(merged_list):
            output_writer.writerow([i - region_length / 2] + base_list)
