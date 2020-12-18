import csv
import math
import multiprocessing
import os
import sys
from collections import defaultdict
from pathlib import Path

from GC_bioinfo.utils.check_dependencies import check_dependencies
from GC_bioinfo.utils.constants import rna_blacklist_file, hg38_chrom_sizes_file, annotation_file
from GC_bioinfo.utils.make_five_and_three_dict import build_counts_dict
from GC_bioinfo.utils.make_five_prime_bed_file import make_five_bed_file
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.make_three_prime_bed_file import make_three_bed_file
from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.run_bedtools_coverage import run_coverage
from GC_bioinfo.utils.run_bedtools_subtract import run_subtract
from GC_bioinfo.utils.verify_bed_file import verify_bed_files


def make_search_regions_list(regions_filename, annotation_extension):
    search_regions_dict = defaultdict(list)
    annotations_dict = {}
    # This makes the search regions
    with open(regions_filename, 'r') as file:
        for line in file:
            chromosome, left, right, strand, gene_name, met_left, met_right = line.split()

            if strand == "+":
                # We want to go from the left position to the met_left position
                search_regions_dict[chromosome].append([chromosome, int(left) - annotation_extension, met_left, gene_name, 0, strand])
            else:
                search_regions_dict[chromosome].append([chromosome, met_right, int(right) + annotation_extension, gene_name, 0, strand])

            annotations_dict[gene_name] = [chromosome, left, right, gene_name, 0, strand]

    return search_regions_dict, annotations_dict


def map_tsrs_to_search_regions(tsr_filename, search_regions_dict):
    gene_tsr_dict = defaultdict(list)
    flow_through_tsrs = []

    search_regions_filename = generate_random_filename()

    with open(search_regions_filename, 'w') as file:
        for chrom in search_regions_dict:
            for region in search_regions_dict[chrom]:
                file.write("\t".join([str(val) for val in region]) + "\n")

    mapped_tsrs_filename = generate_random_filename()
    flow_through_tsrs_filename = generate_random_filename()

    os.system("bedtools intersect -s -a " + search_regions_filename + " -b " + tsr_filename + " -wa -wb > " + mapped_tsrs_filename)
    os.system("bedtools intersect -s -v -a " + tsr_filename + " -b " + search_regions_filename + " > " + flow_through_tsrs_filename)

    with open(mapped_tsrs_filename) as file:
        for line in file:
            gene_name = line.split()[3]
            mapped_tsr = line.split()[6:]
            tsr_chromosome, tsr_left, tsr_right, tsr_read_sum, tsr_strength, tsr_strand, tss_left, tss_right, \
            tss_strength, avg_tss, max_tss_mins_avg_tss = mapped_tsr

            gene_tsr_dict[gene_name].append(
                [tsr_chromosome, tsr_left, tsr_right, tsr_read_sum, tsr_strength, tsr_strand, avg_tss])

    with open(flow_through_tsrs_filename) as file:
        for line in file:
            tsr_chromosome, tsr_left, tsr_right, tsr_read_sum, tsr_strength, tsr_strand, tss_left, tss_right, \
            tss_strength, avg_tss, max_tss_mins_avg_tss = line.split()
            flow_through_tsrs.append(
                [tsr_chromosome, tsr_left, tsr_right, tsr_read_sum, tsr_strength, tsr_strand, avg_tss])

    remove_files(search_regions_filename, mapped_tsrs_filename, flow_through_tsrs_filename)

    return gene_tsr_dict, flow_through_tsrs


def find_max_tsr_in_search_region(gene_tsr_dict):
    max_tsrs_dict = {}
    non_max_tsrs_dict = defaultdict(list)

    for gene_name in gene_tsr_dict:
        for tsr in gene_tsr_dict[gene_name]:
            tsr_chromosome, tsr_left, tsr_right, _, tsr_counts, tsr_strand, avg_tss = tsr

            if gene_name not in max_tsrs_dict:
                max_tsrs_dict[gene_name] = tsr
            elif int(tsr_counts) > int(max_tsrs_dict[gene_name][4]):
                # If the current max is larger than the previous max, add the old tsr to the non_max
                non_max_tsrs_dict[gene_name].append(max_tsrs_dict[gene_name])

                # Put in the new one
                max_tsrs_dict[gene_name] = tsr
            else:
                non_max_tsrs_dict[gene_name].append(tsr)

    return max_tsrs_dict, non_max_tsrs_dict


def define_pause_regions_and_gene_bodies(paused_region_filename, gene_body_region_filename, max_tsrs_dict,
                                         annotations_dict, truQuant_regions_dict):
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
    temp_annotations_dict = defaultdict(lambda: {"+": [], "-": []})
    mapped_flow_through_tsrs_dict = {}

    for gene_name in annotations_dict:
        chromosome, left, right, gene_name, _, strand = annotations_dict[gene_name]

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


def make_blacklisted_regions(blacklist_filename, annotations_dict, max_tsrs_dict, non_max_tsrs_dict, flow_through_tsrs,
                             percent_for_blacklisting):
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

                if int(tsr_counts) >= percent_for_blacklisting * max_tsr_counts:
                    # It should be blacklisted
                    blacklisted_tsrs.append(tsr)

    # Loop through the non_max_tsrs and find the ones we need to blacklist
    for gene_name in non_max_tsrs_dict:
        for tsr in non_max_tsrs_dict[gene_name]:
            tsr_chromosome, tsr_left, tsr_right, _, tsr_counts, tsr_strand, avg_tss = tsr

            if gene_name in max_tsrs_dict:
                max_tsr_counts = int(max_tsrs_dict[gene_name][-3])

                if int(tsr_counts) >= percent_for_blacklisting * max_tsr_counts:
                    blacklisted_tsrs.append(tsr[:-1])

    with open(blacklist_filename, 'w') as output_file:
        tsv_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

        for tsr in blacklisted_tsrs:
            tsv_writer.writerow(tsr)


# ----------------------------------------------     Quantitation     -------------------------------------------------#

def get_counts_in_paused_region(regions_filename, blacklisted_sequencing_file):
    five_bed_filename = make_five_bed_file(blacklisted_sequencing_file)

    # Run bedtools coverage on the 5' bed file
    random_filename = run_coverage(regions_filename, five_bed_filename)

    # Not using a defaultdict because multiprocessing does not like it
    indv_gene_counts_dict = {}

    with open(random_filename) as file:
        for line in file:
            counts = int(line.split()[-4])
            gene_name = line.split()[3]

            if gene_name not in indv_gene_counts_dict:
                indv_gene_counts_dict[gene_name] = {"Pause": -1, "Body": -1}

            indv_gene_counts_dict[gene_name]["Pause"] = counts

    remove_files(five_bed_filename, random_filename)
    return indv_gene_counts_dict


def get_counts_in_gene_bodies(regions_filename, blacklisted_sequencing_file, indv_gene_counts_dict):
    three_bed_filename = make_three_bed_file(blacklisted_sequencing_file)

    # Generate a random filename
    random_filename = run_coverage(regions_filename, three_bed_filename)

    with open(random_filename) as file:
        for line in file:
            counts = line.split()[-4]
            gene_name = line.split()[3]

            indv_gene_counts_dict[gene_name]["Body"] = counts

    remove_files(three_bed_filename, random_filename)
    return indv_gene_counts_dict


def get_region_data(region, five_prime_counts_dict):
    chromosome, left, right, strand = region

    five_prime_sum = 0
    tsr_position_sum = 0
    stdev_weighted_pause_region_center = 0

    max_tss_position = int(left)
    max_tss_counts = 0

    # Loop through each base in the region
    for i in range(int(left), int(right) + 1):
        height = five_prime_counts_dict[chromosome][strand][i]

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
        height = five_prime_counts_dict[chromosome][strand][i]
        position = i - int(left)

        stdev_weighted_pause_region_center += ((position - (weighted_pause_region_center - int(left))) ** 2) * height

    stdev_weighted_pause_region_center = math.sqrt(stdev_weighted_pause_region_center / five_prime_sum)

    return [five_prime_sum, max_tss_position, max_tss_counts, weighted_pause_region_center,
            stdev_weighted_pause_region_center]


def gather_data(sequencing_file, blacklist_filename, annotated_dataset, paused_region_filename,
                gene_body_region_filename, truQuant_regions_dict):
    region_data_dict = {}
    # We need to blacklist the data before running the program
    blacklisted_sequencing_filename = generate_random_filename()

    run_subtract(sequencing_file, rna_blacklist_file, blacklist_filename,
                 output_filename=blacklisted_sequencing_filename)

    indv_gene_counts_dict = get_counts_in_paused_region(paused_region_filename, blacklisted_sequencing_filename)
    get_counts_in_gene_bodies(gene_body_region_filename, blacklisted_sequencing_filename, indv_gene_counts_dict)

    # Only get the region data from the dataset which was annotated
    if annotated_dataset:
        five_prime_counts_dict, _ = build_counts_dict(sequencing_file)

        for gene in truQuant_regions_dict:
            region_data_dict[gene] = get_region_data(truQuant_regions_dict[gene]["Pause"], five_prime_counts_dict)

    remove_files(blacklisted_sequencing_filename)

    return sequencing_file, indv_gene_counts_dict, region_data_dict


def combine_indv_gene_counts_dict(count_information_list, sequencing_files):
    # Need to combine the dictionaries to fit the gene_counts dict format
    ret_region_data_dict = {}
    gene_counts_dict = defaultdict(lambda: {"Pause": [-1] * len(sequencing_files), "Body": [-1] * len(sequencing_files)})

    for i, tup in enumerate(count_information_list):
        # Tup has the sequencing_file, indv_gene_counts_dict, and region_data_dict
        sequencing_file, indv_gene_counts_dict, region_data_dict = tup

        for gene_name in indv_gene_counts_dict:
            gene_counts_dict[gene_name]["Pause"][i] = indv_gene_counts_dict[gene_name]["Pause"]
            gene_counts_dict[gene_name]["Body"][i] = indv_gene_counts_dict[gene_name]["Body"]

        if region_data_dict != {}:
            ret_region_data_dict = region_data_dict

    return ret_region_data_dict, gene_counts_dict


def output_data(output_filename, region_data_dict, gene_counts_dict, truQuant_regions_dict, sequencing_files):
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
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo truQuant <Sequencing Files>\n")
    sys.stderr.write("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/truQuant.rst\n")


def parse_input(args):
    if len(args) == 0:
        print_usage()
        sys.exit(1)

    seq_files = args
    verify_bed_files(seq_files)

    return seq_files


def run_truQuant(annotation_extension, percent_for_blacklisting, truQuant_regions_dict, output_directory,
                 sequencing_files):

    # Blacklist the file first
    blacklisted_first_sequencing_file = sequencing_files[0].replace(".bed", "-blacklisted.bed")
    run_subtract(sequencing_files[0], rna_blacklist_file, output_filename=blacklisted_first_sequencing_file)


    # Define necessary files
    tsr_file = blacklisted_first_sequencing_file.replace(".bed", "_150_20_30_600-TSR.tab")

    paused_region_filename = output_directory + os.path.basename(tsr_file).replace('-TSR.tab', '-paused_regions.bed')
    gene_body_region_filename = output_directory + os.path.basename(tsr_file).replace('-TSR.tab',
                                                                                      '-gene_body_regions.bed')
    blacklist_filename = output_directory + os.path.basename(tsr_file).replace('-TSR.tab', '-blacklisted_regions.bed')
    output_filename = output_directory + os.path.basename(sequencing_files[0].replace('.bed', '-truQuant_output.txt'))


    # Run tsrFinder on the first file
    os.system("GC_bioinfo tsrFinder " + blacklisted_first_sequencing_file + " 150 20 30 600 " + hg38_chrom_sizes_file)

    # 1: Make the regions we are going to be searching for max TSSs in max TSRs
    search_regions_dict, annotations_dict = make_search_regions_list(annotation_file, annotation_extension)

    # 2: Make the pause regions and gene bodies
    gene_tsr_dict, flow_through_tsrs = map_tsrs_to_search_regions(tsr_file, search_regions_dict)
    max_tsrs_dict, non_max_tsrs_dict = find_max_tsr_in_search_region(gene_tsr_dict)
    define_pause_regions_and_gene_bodies(paused_region_filename, gene_body_region_filename,
                                         max_tsrs_dict, annotations_dict, truQuant_regions_dict)
    make_blacklisted_regions(blacklist_filename, annotations_dict, max_tsrs_dict, non_max_tsrs_dict, flow_through_tsrs,
                             percent_for_blacklisting)

    # 3. Get the number of counts (multithreaded)

    pool = multiprocessing.Pool(processes=len(sequencing_files))

    count_information_list = pool.starmap(gather_data,
                                          [(sequencing_file, blacklist_filename, i == 0, paused_region_filename,
                                            gene_body_region_filename, truQuant_regions_dict) for i, sequencing_file in enumerate(sequencing_files)])

    region_data_dict, gene_counts_dict = combine_indv_gene_counts_dict(count_information_list, sequencing_files)

    # Output the data
    output_data(output_filename, region_data_dict, gene_counts_dict, truQuant_regions_dict, sequencing_files)


def main(args):
    """
    truQuant <Sequencing Files>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/truQuant.rst

    :param args: list of the sequencing files
    :type args: list
    :return:
    """
    check_dependencies("bedtools")

    sequencing_files = parse_input(args)

    annotation_extension = 1000
    percent_for_blacklisting = 0.3
    output_directory = str(Path(sequencing_files[0]).parent) + "/"

    truQuant_regions_dict = {}

    run_truQuant(annotation_extension, percent_for_blacklisting, truQuant_regions_dict, output_directory,
                 sequencing_files)


if __name__ == "__main__":
    main(sys.argv[1:])
