'''
This file will run truQuant annotation and quantitation together and make the huge output file
June 20, 2020
Geoff Collins
'''

import csv
import sys
import os
import math
import threading
from pathlib import Path

from utils.make_five_prime_bed_file import make_five_bed_file
from utils.make_three_prime_bed_file import make_three_bed_file
from utils.make_five_and_three_dict import build_counts_dict
from utils.make_random_filename import generate_random_filename
from utils.run_bedtools_coverage import run_coverage
from utils.run_bedtools_subtract import run_subtract
from utils.remove_files import remove_files
from utils.verify_bed_file import verify_bed_files
from utils.check_dependencies import check_dependencies
from utils.constants import rna_blacklist_file, hg38_chrom_sizes_file, annotation_file


def make_search_regions_list(regions_filename):
    search_regions_dict = {}
    annotations_dict = {}
    # This makes the search regions
    with open(regions_filename, 'r') as file:
        for line in file:
            chromosome, left, right, strand, gene_name, met_left, met_right = line.split()

            if strand == "+":
                # We want to go from the left position to the met_left position
                if chromosome in search_regions_dict:
                    search_regions_dict[chromosome].append(
                        [chromosome, int(left) - annotation_extension, met_left, gene_name, 0, strand])
                else:
                    search_regions_dict[chromosome] = [
                        [chromosome, int(left) - annotation_extension, met_left, gene_name, 0, strand]]
            else:
                if chromosome in search_regions_dict:
                    search_regions_dict[chromosome].append(
                        [chromosome, met_right, int(right) + annotation_extension, gene_name, 0, strand])
                else:
                    search_regions_dict[chromosome] = [
                        [chromosome, met_right, int(right) + annotation_extension, gene_name, 0, strand]]

            annotations_dict[gene_name] = [chromosome, left, right, gene_name, 0, strand]

    return search_regions_dict, annotations_dict


def map_tsrs_to_search_regions(tsr_filename, search_regions_dict):
    gene_tsr_dict = {}
    flow_through_tsrs = []

    with open(tsr_filename, 'r') as file:
        # Loop through the TSR file
        for i, line in enumerate(file):
            if i != 0:
                tsr_chromosome, tsr_left, tsr_right, tsr_read_sum, tsr_strength, tsr_strand, tss_left, tss_right, \
                max_tss, tss_strength, avg_tss, max_tss_mins_avg_tss, std_dev_avg_tss = line.split()

                has_mapped = False

                if tsr_chromosome in search_regions_dict:
                    for region in search_regions_dict[tsr_chromosome]:
                        chromosome, left, right, gene_name, _, strand = region

                        # If there is a TSR as that base pair, add it to the gene_tsr dict
                        if tsr_strand == strand and not (int(tsr_right) < int(left) or int(tsr_left) > int(right)):
                            has_mapped = True
                            if gene_name not in gene_tsr_dict:
                                gene_tsr_dict[gene_name] = [line.split()[0:6] + [avg_tss]]
                            else:
                                gene_tsr_dict[gene_name].append(line.split()[0:6] + [avg_tss])

                        # If the left of the region is past the right side of the TSR,
                        # we don't need to search the rest of the file
                        if int(left) > int(tsr_right):
                            break

                if not has_mapped:
                    flow_through_tsrs.append(line.split()[0:6] + [avg_tss])

    return gene_tsr_dict, flow_through_tsrs


def find_max_tsr_in_search_region(gene_tsr_dict):
    max_tsrs_dict = {}
    non_max_tsrs_dict = {}

    for gene_name in gene_tsr_dict:
        for tsr in gene_tsr_dict[gene_name]:
            tsr_chromosome, tsr_left, tsr_right, _, tsr_counts, tsr_strand, avg_tss = tsr

            if gene_name not in max_tsrs_dict:
                max_tsrs_dict[gene_name] = tsr
            elif int(tsr_counts) > int(max_tsrs_dict[gene_name][4]):
                # If the current max is larger than the previous max, add the old tsr to the non_max
                if gene_name not in non_max_tsrs_dict:
                    non_max_tsrs_dict[gene_name] = [max_tsrs_dict[gene_name]]
                else:
                    non_max_tsrs_dict[gene_name].append(max_tsrs_dict[gene_name])

                # Put in the new one
                max_tsrs_dict[gene_name] = tsr
            else:
                if gene_name not in non_max_tsrs_dict:
                    non_max_tsrs_dict[gene_name] = [tsr]
                else:
                    non_max_tsrs_dict[gene_name].append(tsr)

    return max_tsrs_dict, non_max_tsrs_dict


def define_pause_regions_and_gene_bodies(paused_region_filename, gene_body_region_filename, max_tsrs_dict,
                                         annotations_dict):
    # This function will define the pause region as the max TSR for each gene and gene bodies as the end of
    # the pause region to the TES of the annotation

    with open(paused_region_filename, 'w') as paused_file:
        with open(gene_body_region_filename, 'w') as gene_bodies_file:
            gene_bodies_writer = csv.writer(gene_bodies_file, delimiter='\t', lineterminator='\n')
            paused_writer = csv.writer(paused_file, delimiter='\t', lineterminator='\n')

            # Loop through each gene and print out the positions
            for gene in max_tsrs_dict:
                if max_tsrs_dict[gene][0] == annotations_dict[gene][0]:
                    # If the annotation is on the same chromosome as the TSR
                    # I include this is because the VAMP7 gene is annotated on both the X and Y chromosome

                    # The paused region is the expansion of the avgTSS
                    if max_tsrs_dict[gene][-2] == "+":
                        pause_left = int(round(float(max_tsrs_dict[gene][-1]))) - 75
                        pause_right = int(round(float(max_tsrs_dict[gene][-1]))) + 75

                        gene_left = pause_right + 1
                        gene_right = annotations_dict[gene][2]  # The right position of the annotation
                    else:
                        pause_left = int(round(float(max_tsrs_dict[gene][-1]))) - 75
                        pause_right = int(round(float(max_tsrs_dict[gene][-1]))) + 75

                        gene_left = annotations_dict[gene][1]  # The left position of the annotation
                        gene_right = pause_left - 1

                    paused_writer.writerow(
                        [max_tsrs_dict[gene][0]] + [pause_left, pause_right] + [gene] + max_tsrs_dict[gene][4:-1])

                    # Now we just write the gene body to the file
                    gene_bodies_writer.writerow(
                        [max_tsrs_dict[gene][0]] + [gene_left, gene_right] + [gene] + max_tsrs_dict[gene][4:-1])

                    # Fill the truQuant_regions_dict

                    if gene not in truQuant_regions_dict:
                        truQuant_regions_dict[gene] = {"Pause": [], "Body": []}

                    truQuant_regions_dict[gene]["Pause"] = [max_tsrs_dict[gene][0]] + [pause_left, pause_right] + [
                        max_tsrs_dict[gene][-2]]
                    truQuant_regions_dict[gene]["Body"] = [gene_left, gene_right, int(gene_right) - int(gene_left)]


def map_flow_through_tsrs(annotations_dict, flow_through_tsrs):
    # To make this function faster, we split the annotations dict by chromosome and strand
    temp_annotations_dict = {}
    mapped_flow_through_tsrs_dict = {}

    for gene_name in annotations_dict:
        chromosome, left, right, gene_name, _, strand = annotations_dict[gene_name]

        if chromosome not in temp_annotations_dict:
            temp_annotations_dict[chromosome] = {"+": [], "-": []}

        temp_annotations_dict[chromosome][strand].append([left, right, gene_name])

    for tsr in flow_through_tsrs:
        tsr_chromosome, tsr_left, tsr_right, _, tsr_counts, tsr_strand, avg_tss = tsr

        for region in temp_annotations_dict[tsr_chromosome][tsr_strand]:
            left, right, gene_name = region

            # If there is a TSR as that base pair, add it to the gene_tsr dict
            if not (int(tsr_right) < int(left) or int(tsr_left) > int(right)):
                if gene_name not in mapped_flow_through_tsrs_dict:
                    mapped_flow_through_tsrs_dict[gene_name] = [
                        [tsr_chromosome, tsr_left, tsr_right, gene_name, tsr_counts, tsr_strand]]
                else:
                    mapped_flow_through_tsrs_dict[gene_name].append(
                        [tsr_chromosome, tsr_left, tsr_right, gene_name, tsr_counts, tsr_strand])

    return mapped_flow_through_tsrs_dict


def make_blacklisted_regions(blacklist_filename, annotations_dict, max_tsrs_dict, non_max_tsrs_dict, flow_through_tsrs):
    # We are blacklisting all non max TSRs that are not inside ANY paused region

    blacklisted_tsrs = []

    # Map all of the TSRs that did not get mapped to search regions
    mapped_flow_through_tsrs_dict = map_flow_through_tsrs(annotations_dict, flow_through_tsrs)

    # Go through each gene and add TSRs to the blacklist that are at least 30% of the max TSR
    for gene_name in mapped_flow_through_tsrs_dict:
        for tsr in mapped_flow_through_tsrs_dict[gene_name]:
            tsr_chromosome, tsr_left, tsr_right, gene_name, tsr_counts, tsr_strand = tsr

            if gene_name in max_tsrs_dict:
                max_tsr_counts = int(max_tsrs_dict[gene_name][-3])

                if int(tsr_counts) >= 0.3 * max_tsr_counts:
                    # It should be blacklisted
                    blacklisted_tsrs.append(tsr)

    # Loop through the non_max_tsrs and find the ones we need to blacklist
    for gene_name in non_max_tsrs_dict:
        for tsr in non_max_tsrs_dict[gene_name]:
            tsr_chromosome, tsr_left, tsr_right, _, tsr_counts, tsr_strand, avg_tss = tsr

            if gene_name in max_tsrs_dict:
                max_tsr_counts = int(max_tsrs_dict[gene_name][-3])

                if int(tsr_counts) >= 0.3 * max_tsr_counts:
                    blacklisted_tsrs.append(tsr[:-1])

    with open(blacklist_filename, 'w') as output_file:
        tsv_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

        for tsr in blacklisted_tsrs:
            tsv_writer.writerow(tsr)


# ----------------------------------------------     Quantitation     -------------------------------------------------#

def get_counts_in_paused_region(regions_filename, blacklisted_sequencing_file, sequencing_file):
    five_bed_filename = make_five_bed_file(blacklisted_sequencing_file)

    # Run bedtools coverage on the 5' bed file
    random_filename = run_coverage(regions_filename, five_bed_filename)

    with open(random_filename) as file:
        for line in file:
            counts = int(line.split()[-4])
            gene_name = line.split()[3]

            if gene_name not in gene_counts_dict:
                gene_counts_dict[gene_name] = {"Pause": [-1] * len(sequencing_files),
                                               "Body": [-1] * len(sequencing_files)}

            # We are going to place the current data in the position that it is found in the sequencing_files list
            for i, seq_file in enumerate(sequencing_files):
                if sequencing_file == seq_file:
                    gene_counts_dict[gene_name]["Pause"][i] = counts
                    break

    remove_files(five_bed_filename, random_filename)


def get_counts_in_gene_bodies(regions_filename, blacklisted_sequencing_file, sequencing_file):
    three_bed_filename = make_three_bed_file(blacklisted_sequencing_file)

    # Generate a random filename
    random_filename = run_coverage(regions_filename, three_bed_filename)

    with open(random_filename) as file:
        for line in file:
            counts = line.split()[-4]
            gene_name = line.split()[3]

            # We are going to place the current data in the position that it is found in the sequencing_files list
            for i, seq_file in enumerate(sequencing_files):
                if sequencing_file == seq_file:
                    gene_counts_dict[gene_name]["Body"][i] = counts
                    break

    remove_files(three_bed_filename, random_filename)


def get_region_data(region, five_prime_counts_dict):
    chromosome, left, right, strand = region

    five_prime_sum = 0
    tsr_position_sum = 0
    stdev_weighted_pause_region_center = 0

    max_tss_position = int(left)
    max_tss_counts = 0

    # Loop through each base in the region
    for i in range(int(left), int(right) + 1):
        if i in five_prime_counts_dict[chromosome][strand]:
            height = five_prime_counts_dict[chromosome][strand][i]
        else:
            height = 0

        position = i - int(left)

        tsr_position_sum += position * height
        five_prime_sum += height

        if height > max_tss_counts:
            max_tss_position = i
            max_tss_counts = height

    weighted_pause_region_center = int(left) + (tsr_position_sum / five_prime_sum)

    # Round it and make it an integer
    weighted_pause_region_center = int(round(weighted_pause_region_center))

    for i in range(int(left), int(right)):
        if i in five_prime_counts_dict[chromosome][strand]:
            height = five_prime_counts_dict[chromosome][strand][i]
        else:
            height = 0
        position = i - int(left)

        stdev_weighted_pause_region_center += ((position - (weighted_pause_region_center - int(left))) ** 2) * height

    stdev_weighted_pause_region_center = math.sqrt(stdev_weighted_pause_region_center / five_prime_sum)

    return [five_prime_sum, max_tss_position, max_tss_counts, weighted_pause_region_center,
            stdev_weighted_pause_region_center]


def gather_data(sequencing_file, blacklist_filename, annotated_dataset, paused_region_filename,
                gene_body_region_filename):
    # We need to blacklist the data before running the program
    blacklisted_sequencing_filename = generate_random_filename()

    run_subtract(sequencing_file, rna_blacklist_file, blacklist_filename,
                 output_filename=blacklisted_sequencing_filename)

    get_counts_in_paused_region(paused_region_filename, blacklisted_sequencing_filename, sequencing_file)
    get_counts_in_gene_bodies(gene_body_region_filename, blacklisted_sequencing_filename, sequencing_file)

    # Only get the region data from the dataset which was annotated
    if annotated_dataset:
        five_prime_counts_dict, _ = build_counts_dict(sequencing_file)

        for gene in truQuant_regions_dict:
            region_data_dict[gene] = get_region_data(truQuant_regions_dict[gene]["Pause"], five_prime_counts_dict)

    remove_files(blacklisted_sequencing_filename)


def output_data(output_filename):
    pause_region_headers = [file.split("/")[-1] + " Pause Region" for file in sequencing_files]
    gene_body_headers = [file.split("/")[-1] + " Gene Body" for file in sequencing_files]

    with open(output_filename, 'w') as output_file:
        output_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

        output_writer.writerow(["Gene", "Chromosome", "Pause Region Left", "Pause Region Right", "Strand",
                                "Total 5' Reads", "MaxTSS", "MaxTSS 5' Reads", "Weighted Pause Region Center",
                                "STDEV of TSSs", "Gene Body Left", "Gene Body Right", "Gene Body Distance"]
                               + pause_region_headers + gene_body_headers)

        for gene in truQuant_regions_dict:
            # Gene, chrom, left, right, strand
            output_list = [gene] + truQuant_regions_dict[gene]["Pause"]

            # Add the total reads, max TSS, maxTSS 5' reads, avgTSS, stdev
            output_list += region_data_dict[gene]

            # Add the gene body positions and length
            output_list += truQuant_regions_dict[gene]["Body"]

            # Add the data from each of the sequencing files
            output_list += gene_counts_dict[gene]["Pause"] + gene_counts_dict[gene]["Body"]

            output_writer.writerow(output_list)


def print_usage():
    print("Usage: ")
    print("python3 truQuant <Sequencing Files>")
    print("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/truQuant.rst")


def parse_input():
    if len(sys.argv) == 1:
        print_usage()
        sys.exit(1)

    seq_files = sys.argv[1:]
    verify_bed_files(seq_files)

    return seq_files


def run_truQuant():
    # Blacklist the file first
    blacklisted_first_sequencing_file = sequencing_files[0].replace(".bed", "-blacklisted.bed")
    run_subtract(sequencing_files[0], rna_blacklist_file, output_filename=blacklisted_first_sequencing_file)

    # Run tsrFinder on the first file
    os.system("tsrFinderPARALLEL " + blacklisted_first_sequencing_file + " 150 20 30 600 " + hg38_chrom_sizes_file)

    tsr_file = blacklisted_first_sequencing_file.replace(".bed", "_150_20_30_600-TSR.tab")

    paused_region_filename = output_directory + os.path.basename(tsr_file).replace('-TSR.tab', '-paused_regions.bed')
    gene_body_region_filename = output_directory + os.path.basename(tsr_file).replace('-TSR.tab',
                                                                                      '-gene_body_regions.bed')
    blacklist_filename = output_directory + os.path.basename(tsr_file).replace('-TSR.tab', '-blacklisted_regions.bed')
    output_filename = output_directory + os.path.basename(sequencing_files[0].replace('.bed', '-truQuant_output.txt'))

    # 1: Make the regions we are going to be searching for max TSSs in max TSRs
    search_regions_dict, annotations_dict = make_search_regions_list(annotation_file)

    # 2: Make the pause regions and gene bodies
    gene_tsr_dict, flow_through_tsrs = map_tsrs_to_search_regions(tsr_file, search_regions_dict)
    max_tsrs_dict, non_max_tsrs_dict = find_max_tsr_in_search_region(gene_tsr_dict)
    define_pause_regions_and_gene_bodies(paused_region_filename, gene_body_region_filename,
                                         max_tsrs_dict, annotations_dict)
    make_blacklisted_regions(blacklist_filename, annotations_dict, max_tsrs_dict, non_max_tsrs_dict, flow_through_tsrs)

    # 3. Get the number of counts (multithreaded)
    threads_list = []
    for i, sequencing_file in enumerate(sequencing_files):
        annotated_dataset = i == 0
        t = threading.Thread(target=gather_data, args=[sequencing_file, blacklist_filename, annotated_dataset,
                                                       paused_region_filename, gene_body_region_filename])
        t.start()
        threads_list.append(t)

    for thread in threads_list:
        thread.join()

    # Output the data
    output_data(output_filename)


if __name__ == "__main__":
    check_dependencies("bedtools", "tsrFinderPARALLEL")
    sequencing_files = parse_input()

    annotation_extension = 1000
    blacklist_extension = 0
    percent_for_blacklisting = 0.3
    output_directory = str(Path(sequencing_files[0]).parent) + "/"

    gene_counts_dict = {}
    truQuant_regions_dict = {}
    region_data_dict = {}
    run_truQuant()
