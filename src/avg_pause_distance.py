'''
This will take in a a region file and a sequencing file to output the average pausing distance for each gene
'''

# Todo

import sys

from utils.make_transcripts_dict import build_transcripts_dict

# Has keys of chromosomes and values of a dictionary with keys of strand and values of a list of counts
counts_dict = {}

def get_pausing_distances(region_filename):
    all_pause_distances = [0]*101

    # Go through each gene and get the distances from each one
    with open(region_filename) as file:
        for line in file:
            chromosome, left, right, gene_name, fold_change, strand = line.rstrip().split()

            five_prime_position = int(left) + 100

            # Get all of the transcript lengths at this position
            if five_prime_position in counts_dict[chromosome][strand]:
                # Go through all the transcripts that have this 5' end
                for three_prime_end in counts_dict[chromosome][strand][five_prime_position]:
                    transcript_length = abs(five_prime_position - three_prime_end)

                    if transcript_length <= 100 and transcript_length > 17:
                        all_pause_distances[transcript_length] += 1

    return all_pause_distances


if __name__ == '__main__':
    regions_filename = sys.argv[1]


    tfiib_files = [
        "/media/projects/FKBP_DATASETS/2019-03-14-0149_Price_JSantana/DEDUPED/TFIIB_DMSO_2hrs_20190314-dedup_hg38-blacklisted-150.bed",
        "/media/projects/FKBP_DATASETS/2019-03-14-0149_Price_JSantana/DEDUPED/TFIIB_dTAG7_2hrs_20190314-dedup_hg38-blacklisted-150.bed",
        "/media/projects/FKBP_DATASETS/2020-01-17-PriceSantana/DEDUPED/hg38/blacklisted_files/TFIIB_DMSO_20200117-dedup_hg38blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2020-01-17-PriceSantana/DEDUPED/hg38/blacklisted_files/TFIIB_VHL_2hrs_20200117-dedup_hg38blacklisted.bed",
    ]

    taf1_files = [
        "/media/projects/FKBP_DATASETS/2019-12-10-PriceSantana/TAF1_DMSO_2hrs_20191210-dedup-hg38-blacklisted-150.bed",
        "/media/projects/FKBP_DATASETS/2019-12-10-PriceSantana/TAF1_VHL_2hrs_20191210-dedup-hg38-blacklisted-150.bed",
        "/media/projects/FKBP_DATASETS/2019-10-28-PriceSantana/BED/TAF1/TAF1_DMSO_2hrs_20191027-hg38-dedup-blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2019-10-28-PriceSantana/BED/TAF1/TAF1_VHL_2hrs_20191027-hg38-dedup-blacklisted.bed",
    ]

    taf4_files = [
        "/media/projects/FKBP_DATASETS/2019-10-28-PriceSantana/BED/TAF4/TAF4_DMSO_2hrs_20191027-hg38-dedup-blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2019-10-28-PriceSantana/BED/TAF4/TAF4_VHL_2hrs_20191027-hg38-dedup-blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2020-01-17-PriceSantana/DEDUPED/hg38/blacklisted_files/TAF4_DMSO_20200117-dedup_hg38blacklisted.bed",
        "/media/projects/FKBP_DATASETS/2020-01-17-PriceSantana/DEDUPED/hg38/blacklisted_files/TAF4_VHL_2hrs_20200117-dedup_hg38blacklisted.bed"
    ]

    sequencing_files_list = []

    if "taf1" in regions_filename.lower():
        sequencing_files_list.extend(taf1_files)
    if "taf4" in regions_filename.lower():
        sequencing_files_list.extend(taf4_files)
    if "tfiib" in regions_filename.lower():
        sequencing_files_list.extend(tfiib_files)


    print("Region file: " + regions_filename)

    for i, sequencing_filename in enumerate(sequencing_files_list):
        counts_dict = {}
        build_transcripts_dict(sequencing_filename)

        pausing_pileups = get_pausing_distances(regions_filename)

        avg_len = 0
        total_transcripts = 0
        for i in range(len(pausing_pileups)):
            avg_len += pausing_pileups[i] * i
            total_transcripts += pausing_pileups[i]

        avg_len /= total_transcripts

        print("Average pausing distance for\t" + sequencing_filename.split("/")[-1] + "\t" + str(avg_len))
