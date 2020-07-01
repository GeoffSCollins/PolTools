'''
This will take in a a region file and a sequencing file to output the average pausing distance for each gene
'''

# Todo: Test this program

import sys
import csv

from utils import make_transcripts_dict, get_region_length

# Has keys of chromosomes and values of a dictionary with keys of strand and values of a list of counts
counts_dict = {}
output_dict = {}


def most_common(lst):
    return max(set(lst), key=lst.count)


def get_pausing_distances(region_filename, region_length):
    # Go through each gene and get the distances from each one
    with open(region_filename) as file:
        for line in file:
            chromosome, left, right, gene_name, fold_change, strand = line.rstrip().split()

            five_prime_position = int(left) + region_length / 2

            # Get all of the transcript lengths at this position
            if five_prime_position in counts_dict[chromosome][strand]:
                most_common_pausing_position = most_common(counts_dict[chromosome][strand][five_prime_position])
                most_common_pausing_distance = abs(five_prime_position - most_common_pausing_position)

                # Add the gene to the output dict
                if gene_name not in output_dict:
                    output_dict[gene_name] = {}

                if sequencing_filename not in output_dict[gene_name]:
                    output_dict[gene_name][sequencing_filename] = -1

                output_dict[gene_name][sequencing_filename] = most_common_pausing_distance

            else:
                # Add the gene to the output dict
                if gene_name not in output_dict:
                    output_dict[gene_name] = {}

                if sequencing_filename not in output_dict[gene_name]:
                    output_dict[gene_name][sequencing_filename] = "N/A"


if __name__ == '__main__':
    regions_filename = sys.argv[1]
    sequencing_files_list = sys.argv[2:]

    if len(sequencing_files_list) == 0:
        print("You did not provide any sequencing files! Exiting ...")
        sys.exit(1)

    region_length = get_region_length.determine_region_length()

    for i, sequencing_filename in enumerate(sequencing_files_list):
        counts_dict = {}
        make_transcripts_dict.build_transcripts_dict(sequencing_filename)
        get_pausing_distances(regions_filename, region_length)

    # Output to file
    with open(regions_filename.replace(".bed", "maxTPS.txt"), 'w') as file:
        output_writer = csv.writer(file, delimiter='\t', lineterminator='\n')
        for i, gene in enumerate(output_dict):
            if i == 0:
                # We print the headers
                output_writer.writerow(["Gene"] + [x.split("/")[-1] for x in output_dict[gene].keys()])

            output_writer.writerow([gene] + [output_dict[gene][sequencing_filename] for sequencing_filename in output_dict[gene]])
