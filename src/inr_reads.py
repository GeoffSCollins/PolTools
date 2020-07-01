'''
This will take in a a region file and a sequencing file to output the average pausing distance for each gene
'''

# Todo

import sys
import csv

# Has keys of chromosomes and values of a dictionary with keys of strand and values of a list of counts
counts_dict = {}


def build_counts_dict(sequencing_filename):
    with open(sequencing_filename, 'r') as file:
        for i, line in enumerate(file):
            chromosome, left, right, _, _, strand = line.split()

            if strand == "+":
                five_prime_position = int(left)
            else:
                five_prime_position = int(right)

            if chromosome not in counts_dict:
                counts_dict[chromosome] = {}

            if strand not in counts_dict[chromosome]:
                counts_dict[chromosome][strand] = {}

            if five_prime_position not in counts_dict[chromosome][strand]:
                counts_dict[chromosome][strand][five_prime_position] = 0

            counts_dict[chromosome][strand][five_prime_position] += 1


def get_counts(counts_dict, regions_filename):
    # This will return a list of the number of counts

    counts_list = []

    # Loop through each region
    with open(regions_filename) as file:
        for line in file:
            chromosome, left, right, gene_name, _, strand = line.split()

            # Only getting counts from -5 to +5, so we
            loop_left = int(left) + 95
            loop_right = int(right) - 95

            max_counts = 0
            for i in range(loop_left, loop_right + 1):
                if i in counts_dict[chromosome][strand]:
                    if counts_dict[chromosome][strand][i] > max_counts:
                        max_counts = counts_dict[chromosome][strand][i]

            # We have the counts, so just append it
            counts_list.append((gene_name, max_counts))

    return counts_list


if __name__ == '__main__':
    regions_filename = sys.argv[1]

    tfiib_files = [
        "/media/projects/FKBP_DATASETS/2019-03-14-0149_Price_JSantana/DEDUPED/TFIIB_DMSO_2hrs_20190314-dedup_hg38-blacklisted-150.bed",
        "/media/projects/FKBP_DATASETS/2019-03-14-0149_Price_JSantana/DEDUPED/TFIIB_dTAG7_2hrs_20190314-dedup_hg38-blacklisted-150.bed",
        "/media/projects/FKBP_DATASETS/2019-03-14-0149_Price_JSantana/DEDUPED/TFIIB_dTAG7_Flavo_2hrs_20190314-dedup_hg38-blacklisted-150.bed",
        "/media/projects/FKBP_DATASETS/2020-01-17-PriceSantana/DEDUPED/hg38/blacklisted_files/TFIIB_DMSO_20200117-dedup_hg38blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2020-01-17-PriceSantana/DEDUPED/hg38/blacklisted_files/TFIIB_VHL_2hrs_20200117-dedup_hg38blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2020-01-17-PriceSantana/DEDUPED/hg38/blacklisted_files/TFIIB_VHL_Flavo_2hrs_20200117-dedup_hg38blacklisted.bed"
    ]

    taf1_files = [
        "/media/projects/FKBP_DATASETS/2019-12-10-PriceSantana/TAF1_DMSO_2hrs_20191210-dedup-hg38-blacklisted-150.bed",
        "/media/projects/FKBP_DATASETS/2019-12-10-PriceSantana/TAF1_VHL_2hrs_20191210-dedup-hg38-blacklisted-150.bed",
        "/media/projects/FKBP_DATASETS/2019-12-10-PriceSantana/TAF1_VHL_Flavo_2hrs_20191210-dedup-hg38-blacklisted-150.bed",
        "/media/projects/FKBP_DATASETS/2019-10-28-PriceSantana/BED/TAF1/TAF1_DMSO_2hrs_20191027-hg38-dedup-blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2019-10-28-PriceSantana/BED/TAF1/TAF1_VHL_2hrs_20191027-hg38-dedup-blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2019-10-28-PriceSantana/BED/TAF1/TAF1_VHL_Flavo_2hrs_20191027-hg38-dedup-blacklisted.bed"
    ]

    taf4_files = [
        "/media/projects/FKBP_DATASETS/2019-10-28-PriceSantana/BED/TAF4/TAF4_DMSO_2hrs_20191027-hg38-dedup-blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2019-10-28-PriceSantana/BED/TAF4/TAF4_VHL_0.5hrs_20191027-hg38-dedup-blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2019-10-28-PriceSantana/BED/TAF4/TAF4_VHL_2hrs_20191027-hg38-dedup-blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2020-01-17-PriceSantana/DEDUPED/hg38/blacklisted_files/TAF4_DMSO_20200117-dedup_hg38blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2020-01-17-PriceSantana/DEDUPED/hg38/blacklisted_files/TAF4_VHL_0.5hrs-dedup_hg38blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2020-01-17-PriceSantana/DEDUPED/hg38/blacklisted_files/TAF4_VHL_2hrs_20200117-dedup_hg38blacklisted.bed"
    ]

    sequencing_files_list = []

    # Will have keys of gene names and values of another dictionary with keys of dataset and values of counts
    genes_dict = {}

    if "taf1" in regions_filename.lower():
        sequencing_files_list.extend(taf1_files)
    if "taf4" in regions_filename.lower():
        sequencing_files_list.extend(taf4_files)
    if "tfiib" in regions_filename.lower():
        sequencing_files_list.extend(tfiib_files)

    data = []

    for i, sequencing_filename in enumerate(sequencing_files_list):
        counts_dict = {}
        build_counts_dict(sequencing_filename)

        tmp_list = get_counts(counts_dict, regions_filename)

        # Looks like [ (x, x) ]

        # Now we want to get the counts into the genes dict
        for tup in tmp_list:
            gene_name, counts = tup

            if gene_name not in genes_dict:
                genes_dict[gene_name] = {}

            if sequencing_filename not in genes_dict[gene_name]:
                genes_dict[gene_name][sequencing_filename] = 0

            genes_dict[gene_name][sequencing_filename] = counts


    with open(regions_filename.replace(".bed", "_inr_reads.txt"), 'w') as output_file:
        output_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

        for i, gene in enumerate(genes_dict):
            if i == 0:
                # Write all of the headers
                output_writer.writerow(genes_dict[gene].keys())
                output_writer.writerow([gene] + [genes_dict[gene][dataset] for dataset in genes_dict[gene]])

            else:
                # We write the counts
                output_writer.writerow([gene] + [genes_dict[gene][dataset] for dataset in genes_dict[gene]])
