'''
This python script takes in two files: a regions file (in bed format) and a sequencing file.
It will output a single file which has the number of 5' ends at a distance away from the annotated start site for each strand
'''

# Todo

import sys
import os
import csv

from itertools import chain

five_prime_counts_dict = {}
three_prime_counts_dict = {}


def split_sequencing_file(sequencing_file):
    with open(sequencing_file.replace(".bed", "-FW.bed"), 'w') as fw_file:
        with open(sequencing_file.replace(".bed", "-RV.bed"), 'w') as rv_file:
            with open(sequencing_file) as file:
                for line in file:
                    if "+" in line:
                        fw_file.write(line)
                    else:
                        rv_file.write(line)


def split_region_file(region_file):
    with open(region_file) as file:
        with open(region_file.replace(".bed", "-FW.bed"), 'w') as fw_file:
            with open(region_file.replace(".bed", "-RV.bed"), 'w') as rv_file:
                for line in file:
                    if "+" in line:
                        fw_file.write(line)
                    else:
                        rv_file.write(line)


def remove_split_files(sequencing_file, region_file, rev_region_file):
    os.system("rm -rf " + sequencing_file.replace(".bed", "-FW.bed") + sequencing_file.replace(".bed", "-RV.bed"))
    os.system("rm -rf " + region_file.replace(".bed", "-FW.bed") + region_file.replace(".bed", "-RV.bed"))
    os.system("rm -rf " + rev_region_file.replace(".bed", "-FW.bed") + rev_region_file.replace(".bed", "-RV.bed"))


def get_pileups(region_filename, sequencing_file):
    split_sequencing_file(sequencing_file)
    split_region_file(region_filename)

    # Run bedtools coverage on all of them
    os.system(
        "bedtools coverage -d -a " + region_filename.replace(".bed", "-FW.bed") + " -b " + sequencing_file.replace(
            ".bed", "-FW.bed") + " > ./tmp_pileup_data_FW.bed")

    os.system(
        "bedtools coverage -d -a " + region_filename.replace(".bed", "-RV.bed") + " -b " + sequencing_file.replace(
            ".bed", "-RV.bed") + " > ./tmp_pileup_data_RV.bed")

    # Combine the two files together
    os.system("cat ./tmp_pileup_data_FW.bed ./tmp_pileup_data_RV.bed > ./tmp_combined_pileup_data.bed")

    # Remove the separate output files
    os.system("rm -rf ./tmp_pileup_data_FW.bed ./tmp_pileup_data_RV.bed")

    # Counts list will start with the -500 nt and go to +500 nt
    counts_list = [0] * 1000

    with open("./tmp_combined_pileup_data.bed") as file:
        for i, line in enumerate(file):
            if "+" in line:
                position = int(line.split()[-2]) - 1  # Subtract 1 because position starts at 1
            else:
                position = 999 - int(line.split()[-2])  # Subtract 1 because position starts at 1

            counts_list[position] += int(line.split()[-1])

    avg_list = [0] * 1000

    # Get the averages quick
    for i, val in enumerate(counts_list):
        avg_list[i] = val / len(counts_list)

    # Clean up the regions file

    # We are wanting to get the opposite strand as well

    rev_region_filename = region_filename.replace(".bed", "rev_strands.bed")

    # Make a new version of the regions file
    with open(region_filename) as file:
        with open(rev_region_filename, 'w') as outfile:
            for line in file:
                if "+" in line:
                    outfile.write(line.replace("+", "-"))
                else:
                    outfile.write(line.replace("-", "+"))

    split_region_file(rev_region_filename)

    # Run bedtools coverage on all of them
    os.system(
        "bedtools coverage -d -a " + rev_region_filename.replace(".bed", "-FW.bed") + " -b " + sequencing_file.replace(
            ".bed", "-FW.bed") + " > ./tmp_pileup_data_FW.bed")

    os.system(
        "bedtools coverage -d -a " + rev_region_filename.replace(".bed", "-RV.bed") + " -b " + sequencing_file.replace(
            ".bed", "-RV.bed") + " > ./tmp_pileup_data_RV.bed")

    # Combine the two files together
    os.system("cat ./tmp_pileup_data_FW.bed ./tmp_pileup_data_RV.bed > ./tmp_combined_pileup_data.bed")

    # Remove the separate output files
    os.system("rm -rf ./tmp_pileup_data_FW.bed ./tmp_pileup_data_RV.bed")

    # Counts list will start with the -500 nt and go to +500 nt
    rev_counts_list = [0] * 1000

    with open("./tmp_combined_pileup_data.bed") as file:
        for i, line in enumerate(file):
            if "-" in line:
                position = int(line.split()[-2]) - 1  # Subtract 1 because position starts at 1
            else:
                position = 999 - int(line.split()[-2])  # Subtract 1 because position starts at 1

            rev_counts_list[position] += int(line.split()[-1])

    rev_avg_list = [0] * 1000

    # Get the averages quick
    for i, val in enumerate(rev_counts_list):
        rev_avg_list[i] = val / len(rev_counts_list)

    # We finally clean up the split files and the tmp_combined_pileup_data.bed
    remove_split_files(sequencing_file, region_filename, rev_region_filename)
    os.system("rm -rf ./tmp_combined_pileup_data.bed")

    return avg_list, rev_avg_list


if __name__ == '__main__':
    regions_filename = sys.argv[1]
    sequencing_files_list = sys.argv[2:]

    output_filename = regions_filename.replace(".bed", "_pileup_metaplots.tsv")

    if len(sequencing_files_list) == 0:
        print("You did not provide any sequencing files! Exiting ...")
        sys.exit(1)

    # Output dictionary has keys of sequencing files and values of a list containing average lists
    # The order of the list goes as follows: 5'same, 5'rev, 3'same, 3'rev, pileup
    output_dictionary = {}

    for i, sequencing_file in enumerate(sequencing_files_list):
        print("Working on file " + str(i + 1) + " out of " + str(len(sequencing_files_list)))
        # 4. Now we do the pileup data
        pileup_averages, rev_pileup_averages = get_pileups(regions_filename, sequencing_file)

        l = [pileup_averages, rev_pileup_averages]
        transposed_l = list(map(list, zip(*l)))
        output_dictionary[sequencing_file] = transposed_l

    files = []
    all_data = []
    for file in output_dictionary:
        all_data.append(output_dictionary[file])
        files.append(file)

    # Merge all of the lists together
    merged_list = [list(chain.from_iterable(x)) for x in zip(*all_data)]

    output_filename = regions_filename.replace(".bed", "_pileup_metaplots.tsv")

    # 5. Put the data into a file
    with open(output_filename, 'w') as file:
        output_writer = csv.writer(file, delimiter='\t', lineterminator='\n')

        header = ["Position"]
        # Write the header first
        for file in files:
            header.append(file + " pileup same strand")
            header.append(file + " pileup reverse strand")

        output_writer.writerow(header)

        for i, base_list in enumerate(merged_list):
            output_writer.writerow([i - 500] + base_list)
