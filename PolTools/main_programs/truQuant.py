import csv
import math
import multiprocessing
import os
import sys
import argparse

from collections import defaultdict
from pathlib import Path

from PolTools.utils.check_dependencies import check_dependencies
from PolTools.utils.constants import rna_blacklist_file, hg38_chrom_sizes_file, annotation_file
from PolTools.utils.build_counts_dict import build_counts_dict
from PolTools.utils.make_read_end_file import make_read_end_file
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.remove_files import remove_files
from PolTools.utils.bedtools_utils.run_bedtools_coverage import run_coverage
from PolTools.utils.bedtools_utils.run_bedtools_subtract import run_subtract
from PolTools.utils.verify_bed_file import verify_bed_files


def make_search_regions(regions_filename, annotation_extension):
    search_regions_dict = defaultdict(list)
    annotations_dict = {}
    # Make the search regions from the 5' end extended by annotation_extension to the most downstream methionine
    with open(regions_filename, 'r') as file:
        for line in file:
            chromosome, left, right, strand, gene_name, met_left, met_right = line.split()

            if strand == "+":
                # We want to go from the left position to the met_left position
                search_regions_dict[chromosome].append(
                    [chromosome, str(int(left) - annotation_extension), met_left, gene_name, "0", strand]
                )
            else:
                search_regions_dict[chromosome].append(
                    [chromosome, met_right, str(int(right) + annotation_extension), gene_name, "0", strand]
                )

            annotations_dict[gene_name] = [chromosome, left, right, gene_name, "0", strand]

    return search_regions_dict, annotations_dict


def map_tsrs_to_search_regions(tsr_filename, search_regions_dict):
    # Maps the TSRs found from tsrFinder to the search regions
    # (5' end of gene extended upstream to the most downstream methionine)

    gene_tsr_dict = defaultdict(list)
    flow_through_tsrs = []

    search_regions_filename = generate_random_filename()

    with open(search_regions_filename, 'w') as file:
        for chrom in search_regions_dict:
            for region in search_regions_dict[chrom]:
                file.write("\t".join([str(val) for val in region]) + "\n")

    mapped_tsrs_filename = generate_random_filename()
    flow_through_tsrs_filename = generate_random_filename()

    os.system("bedtools intersect -s -a " + search_regions_filename +
              " -b " + tsr_filename + " -wa -wb > " + mapped_tsrs_filename)
    os.system("bedtools intersect -s -v -a " + tsr_filename +
              " -b " + search_regions_filename + " -wa > " + flow_through_tsrs_filename)

    with open(mapped_tsrs_filename) as file:
        for line in file:
            gene_name = line.split()[3]
            mapped_tsr = line.split()[6:]
            tsr_chromosome, tsr_left, tsr_right, tsr_read_sum, tsr_strength, tsr_strand, tss_left, tss_right, \
            tss_strength, avg_tss = mapped_tsr

            gene_tsr_dict[gene_name].append(
                [tsr_chromosome, tsr_left, tsr_right, tsr_read_sum, tsr_strength, tsr_strand, avg_tss])

    with open(flow_through_tsrs_filename) as file:
        for line in file:
            tsr_chromosome, tsr_left, tsr_right, tsr_read_sum, tsr_strength, tsr_strand, tss_left, tss_right, \
            tss_strength, avg_tss = line.split()
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


def define_pause_regions_and_gene_bodies(region_filenames, max_tsrs_dict, annotations_dict, pause_region_radius):
    # This function will define the pause region as the max TSR for each gene and gene bodies as the end of
    # the pause region to the TES of the annotation

    paused_region_filename, gene_body_region_filename = region_filenames

    truQuant_regions_dict = {}

    with open(paused_region_filename, 'w') as paused_file:
        with open(gene_body_region_filename, 'w') as gene_bodies_file:
            paused_writer = csv.writer(paused_file, delimiter='\t', lineterminator='\n')
            gene_bodies_writer = csv.writer(gene_bodies_file, delimiter='\t', lineterminator='\n')

            # Loop through each gene and print out the positions
            for gene in max_tsrs_dict:
                tsr_chromosome, tsr_left, tsr_right, tsr_read_sum, tsr_counts, tsr_strand, avg_tss = max_tsrs_dict[gene]
                ann_chromosome, ann_left, ann_right, ann_gene_name, ann_score, ann_strand = annotations_dict[gene]

                if tsr_chromosome == ann_chromosome:
                    # I include this is because the VAMP7 gene is annotated on both the X and Y chromosome

                    # The paused region is the expansion of the avgTSS
                    if tsr_strand == "+":
                        pause_left = int(round(float(avg_tss))) - pause_region_radius
                        pause_right = int(round(float(avg_tss))) + pause_region_radius

                        gene_left = pause_right
                        gene_right = int(ann_right)
                    else:
                        pause_left = int(round(float(avg_tss))) - pause_region_radius + 1
                        pause_right = int(round(float(avg_tss))) + pause_region_radius + 1

                        gene_left = int(ann_left)
                        gene_right = pause_left

                    paused_writer.writerow([tsr_chromosome, pause_left, pause_right, gene, tsr_counts, tsr_strand])

                    # Now we just write the gene body to the file
                    gene_bodies_writer.writerow([tsr_chromosome, gene_left, gene_right, gene, tsr_counts, tsr_strand])

                    # Fill the truQuant_regions_dict

                    if gene not in truQuant_regions_dict:
                        truQuant_regions_dict[gene] = {"Pause": [], "Body": []}

                    truQuant_regions_dict[gene]["Pause"] = [tsr_chromosome, pause_left, pause_right, tsr_strand]
                    truQuant_regions_dict[gene]["Body"] = [gene_left, gene_right, int(gene_right) - int(gene_left)]

    return truQuant_regions_dict


def map_flow_through_tsrs(annotations_dict, flow_through_tsrs):
    # Write the flow through tsrs to a temporary file
    tsr_filename = generate_random_filename()

    with open(tsr_filename, 'w') as file:
        for tsr in flow_through_tsrs:
            file.write("\t".join(tsr) + "\n")

    # Write the annotations dict to a temporary file
    annotation_file = generate_random_filename()

    with open(annotation_file, 'w') as file:
        for gene in annotations_dict:
            # Don't write last field because that is the average TSS
            file.write("\t".join(annotations_dict[gene]) + "\n")

    mapped_flow_through_tsrs_dict = defaultdict(list)

    mapped_tsrs_filename = generate_random_filename()

    os.system("bedtools intersect -s -a " + annotation_file + " -b " + tsr_filename + " -wa -wb > " + mapped_tsrs_filename)

    with open(mapped_tsrs_filename) as file:
        for line in file:
            gene_name = line.split()[3]
            mapped_tsr = line.split()[6:]
            tsr_chromosome, tsr_left, tsr_right, _, tsr_counts, tsr_strand, avg_tss = mapped_tsr

            mapped_flow_through_tsrs_dict[gene_name].append(
                [tsr_chromosome, tsr_left, tsr_right, gene_name, tsr_counts, tsr_strand])

    remove_files(tsr_filename, annotation_file, mapped_tsrs_filename)

    return mapped_flow_through_tsrs_dict


def make_blacklisted_regions(blacklist_filename, mapped_tsrs, percent_for_blacklisting):
    max_tsrs_dict, non_max_tsrs_dict, mapped_flow_through_tsrs_dict = mapped_tsrs

    # We are blacklisting all non max TSRs that are not inside ANY paused region
    blacklisted_tsrs = []

    # Go through each gene and add TSRs to the blacklist that are at least 30% of the max TSR
    for gene_name in mapped_flow_through_tsrs_dict:
        for tsr in mapped_flow_through_tsrs_dict[gene_name]:
            tsr_chromosome, tsr_left, tsr_right, tsr_read_len_sum, tsr_counts, tsr_strand = tsr

            if gene_name in max_tsrs_dict:
                max_tsr_counts = int(max_tsrs_dict[gene_name][-3])

                if int(tsr_counts) >= percent_for_blacklisting * max_tsr_counts:
                    # Replace the read_len_sum with the gene name for output
                    tsr[3] = gene_name
                    blacklisted_tsrs.append(tsr)

    # Loop through the non_max_tsrs and find the ones we need to blacklist
    for gene_name in non_max_tsrs_dict:
        for tsr in non_max_tsrs_dict[gene_name]:
            tsr_chromosome, tsr_left, tsr_right, tsr_read_len_sum, tsr_counts, tsr_strand, avg_tss = tsr

            if gene_name in max_tsrs_dict:
                max_tsr_counts = int(max_tsrs_dict[gene_name][-3])

                if int(tsr_counts) >= percent_for_blacklisting * max_tsr_counts:
                    # Replace the read_len_sum with the gene name for output
                    tsr[3] = gene_name
                    blacklisted_tsrs.append(tsr[:-1])

    with open(blacklist_filename, 'w') as output_file:
        tsv_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

        for tsr in blacklisted_tsrs:
            tsv_writer.writerow(tsr)


# ----------------------------------------------     Quantitation     -------------------------------------------------#

def get_counts_in_paused_region(pause_region_filename, blacklisted_sequencing_file):
    five_bed_filename = make_read_end_file(blacklisted_sequencing_file, 'five')

    # Run bedtools coverage on the 5' bed file
    random_filename = run_coverage(pause_region_filename, five_bed_filename)

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
    three_bed_filename = make_read_end_file(blacklisted_sequencing_file, 'three')

    random_filename = run_coverage(regions_filename, three_bed_filename)

    with open(random_filename) as file:
        for line in file:
            counts = int(line.split()[-4])
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

    # We need to add 1 to the max_tss_counts because the genome browser starts at 1
    max_tss_position += 1

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


def gather_data(sequencing_file, blacklist_filename, annotated_dataset, region_filenames, truQuant_regions_dict):
    paused_region_filename, gene_body_region_filename = region_filenames

    region_data_dict = {}
    # We need to blacklist the data before running the program
    blacklisted_sequencing_filename = generate_random_filename()

    run_subtract(sequencing_file, rna_blacklist_file, blacklist_filename, strand_specific=False,
                 output_filename=blacklisted_sequencing_filename)

    indv_gene_counts_dict = get_counts_in_paused_region(paused_region_filename, blacklisted_sequencing_filename)
    get_counts_in_gene_bodies(gene_body_region_filename, blacklisted_sequencing_filename, indv_gene_counts_dict)

    # Only get the region data from the dataset which was annotated
    if annotated_dataset:
        five_prime_counts_dict = build_counts_dict(sequencing_file, "five")

        for gene in truQuant_regions_dict:
            region_data_dict[gene] = get_region_data(truQuant_regions_dict[gene]["Pause"], five_prime_counts_dict)

    remove_files(blacklisted_sequencing_filename)

    return sequencing_file, indv_gene_counts_dict, region_data_dict


def combine_indv_gene_counts_dict(count_information_list, sequencing_files):
    # Need to combine the dictionaries to fit the gene_counts dict format
    ret_region_data_dict = {}
    gene_counts_dict = defaultdict(lambda: {"Pause": [-1] * len(sequencing_files), "Body": [-1] * len(sequencing_files)})

    for i, tup in enumerate(count_information_list):
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


def run_tsrFinder(first_seq_file, max_threads, tsrFinder_parameters):
    tsrFinder_param_string = " ".join([str(val) for val in tsrFinder_parameters])
    tsrFinder_filename_string = "_" + tsrFinder_param_string.replace(" ", "_") + "-TSR.tab"

    blacklisted_first_sequencing_file = first_seq_file.replace(".bed", "-blacklisted.bed")

    tsr_file = blacklisted_first_sequencing_file.replace(".bed", tsrFinder_filename_string)
    run_subtract(first_seq_file, rna_blacklist_file, strand_specific=False, output_filename=blacklisted_first_sequencing_file)

    os.system("PolTools tsrFinder " + blacklisted_first_sequencing_file + " " + tsrFinder_param_string + " " +
              hg38_chrom_sizes_file + " -t " + str(max_threads)
    )

    return tsr_file


def build_annotation(tsr_file, search_regions_dict, output_region_filenames, annotations_dict,
                     blacklist_filename, pause_region_radius, percent_for_blacklisting):

    gene_tsr_dict, flow_through_tsrs = map_tsrs_to_search_regions(tsr_file, search_regions_dict)
    max_tsrs_dict, non_max_tsrs_dict = find_max_tsr_in_search_region(gene_tsr_dict)

    truQuant_regions_dict = define_pause_regions_and_gene_bodies(output_region_filenames, max_tsrs_dict,
                                                                 annotations_dict, pause_region_radius)

    mapped_flow_through_tsrs_dict = map_flow_through_tsrs(annotations_dict, flow_through_tsrs)

    mapped_tsrs = max_tsrs_dict, non_max_tsrs_dict, mapped_flow_through_tsrs_dict

    make_blacklisted_regions(blacklist_filename, mapped_tsrs, percent_for_blacklisting)

    return truQuant_regions_dict


def get_counts(sequencing_files, blacklist_filename, output_region_filenames, truQuant_regions_dict, max_threads):
    with multiprocessing.Pool(processes=max_threads) as pool:

        count_information_list = pool.starmap(gather_data,
                                              [(sequencing_file, blacklist_filename, i == 0, output_region_filenames,
                                                truQuant_regions_dict) for i, sequencing_file in
                                               enumerate(sequencing_files)])

    region_data_dict, gene_counts_dict = combine_indv_gene_counts_dict(count_information_list, sequencing_files)

    return region_data_dict, gene_counts_dict


def run_truQuant(annotation_extension, percent_for_blacklisting, sequencing_files, pause_region_radius, tsrFinder_parameters, max_threads):

    output_directory = str(Path(sequencing_files[0]).parent) + "/"

    # Run tsrFinder then define the required files
    tsr_file = run_tsrFinder(sequencing_files[0], max_threads, tsrFinder_parameters)
    tsr_basename = os.path.basename(tsr_file)

    paused_region_filename = output_directory + tsr_basename.replace('-TSR.tab', '-paused_regions.bed')
    gene_body_region_filename = output_directory + tsr_basename.replace('-TSR.tab', '-gene_body_regions.bed')
    blacklist_filename = output_directory + tsr_basename.replace('-TSR.tab', '-blacklisted_regions.bed')
    output_filename = output_directory + os.path.basename(sequencing_files[0].replace('.bed', '-truQuant_output.txt'))

    output_region_filenames = paused_region_filename, gene_body_region_filename

    # Make the regions we are going to be searching for max TSSs in max TSRs
    search_regions_dict, annotations_dict = make_search_regions(annotation_file, annotation_extension)

    # Make the pause regions, gene bodies, and the blacklisted regions
    truQuant_regions_dict = build_annotation(tsr_file, search_regions_dict, output_region_filenames, annotations_dict,
                     blacklist_filename, pause_region_radius, percent_for_blacklisting)

    # Get the number of counts (multithreaded)
    region_data_dict, gene_counts_dict = get_counts(sequencing_files, blacklist_filename, output_region_filenames,
                                                    truQuant_regions_dict, max_threads)

    # Output the data
    output_data(output_filename, region_data_dict, gene_counts_dict, truQuant_regions_dict, sequencing_files)


def parse_input(args):
    def positive_int(num):
        try:
            val = int(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")

        return val

    def between_zero_and_one(number):
        try:
            num = float(number)
            if num < 0 or num > 1:
                raise argparse.ArgumentTypeError(number + " must be positive and less than one")
        except:
            raise argparse.ArgumentTypeError(number + " must be positive and less than one")

        return num

    parser = argparse.ArgumentParser(prog='PolTools truQuant',
                                     description='Annotate and quantify the human genome from PRO-Seq data\n' +
                                     "More information can be found at " +
                                     "https://geoffscollins.github.io/PolTools/truQuant.html")

    parser.add_argument('-a', '--annotation_extension', dest='annotation_extension', metavar='annotation_extension', nargs='?', type=int, default=1000)

    parser.add_argument('-b', '--blacklist_percent', dest='percent_for_blacklisting', metavar='blacklisting_percent', nargs='?', type=between_zero_and_one, default=0.3)

    parser.add_argument('-r', '--pause_radius', dest='pause_region_radius', metavar='pause_region_radius', type=positive_int, nargs='?', default=75)

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?', default=multiprocessing.cpu_count())

    parser.add_argument('seq_file_for_annotation', metavar='sequencing_file_for_annotation', type=str,
                        help='Bed formatted file from a sequencing experiment which will be used to annotate the genome.')

    parser.add_argument('seq_files', metavar='sequencing_files', nargs='*', type=str,
                        help='Bed formatted files from the sequencing experiment')

    parser.add_argument('-d', '--min_seq_depth', dest='min_seq_depth', metavar='min_seq_depth', type=positive_int, nargs='?', default=20)
    parser.add_argument('-m', '--min_avg_transcript_length', dest='min_avg_transcript_length',
                        metavar='min_avg_transcript_length', type=positive_int, nargs='?', default=30)
    parser.add_argument('-l', '--max_fragment_length', dest='max_fragment_length', metavar='max_fragment_length',
                        type=positive_int, nargs='?', default=600)

    args = parser.parse_args(args)

    seq_files = [args.seq_file_for_annotation] + args.seq_files
    annotation_extension = args.annotation_extension
    percent_for_blacklisting = args.percent_for_blacklisting
    pause_region_radius = args.pause_region_radius
    max_threads = args.threads

    # The window size must be 150 bp
    tsrFinder_parameters = 150, args.min_seq_depth, args.min_avg_transcript_length, args.max_fragment_length

    verify_bed_files(seq_files)

    return seq_files, annotation_extension, percent_for_blacklisting, pause_region_radius, tsrFinder_parameters, max_threads


def main(args):
    """
    truQuant <Sequencing Files>
    More information can be found at https://geoffscollins.github.io/PolTools/truQuant.html

    :param args: list of the sequencing files
    :type args: list
    :return:
    """
    check_dependencies("bedtools")

    sequencing_files, annotation_extension, percent_for_blacklisting, pause_region_radius, tsrFinder_parameters, max_threads = parse_input(args)

    run_truQuant(annotation_extension, percent_for_blacklisting, sequencing_files, pause_region_radius, tsrFinder_parameters, max_threads)


if __name__ == "__main__":
    main(sys.argv[1:])
