import csv

from PolTools.main_programs import truQuant
from PolTools.utils.constants import annotation_file
from PolTools.utils.make_random_filename import generate_random_filename


def _define_pause_regions_and_gene_bodies(max_tsrs_dict, annotations_dict, truQuant_regions_dict, downstream_extension):
    # This function will define the pause region as the max TSR for each gene and gene bodies as the end of
    # the pause region to the TES of the annotation

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
                gene_right = str(int(annotations_dict[gene][2]) + downstream_extension)  # The right position of the annotation
            else:
                pause_left = int(round(float(max_tsrs_dict[gene][-1]))) - 75
                pause_right = int(round(float(max_tsrs_dict[gene][-1]))) + 75

                gene_left = str(int(annotations_dict[gene][1]) - downstream_extension)  # The left position of the annotation
                gene_right = pause_left - 1

            # Fill the truQuant_regions_dict
            if gene not in truQuant_regions_dict:
                truQuant_regions_dict[gene] = {"Pause": [], "Body": []}

            truQuant_regions_dict[gene]["Pause"] = [max_tsrs_dict[gene][0]] + [pause_left, pause_right] + [
                max_tsrs_dict[gene][-2]]
            truQuant_regions_dict[gene]["Body"] = [gene_left, gene_right, int(gene_right) - int(gene_left)]


def _make_blacklisted_regions(blacklist_filename, annotations_dict, max_tsrs_dict, non_max_tsrs_dict, flow_through_tsrs,
                             percent_for_blacklisting):
    # We are blacklisting all non max TSRs that are not inside ANY paused region

    blacklisted_tsrs = []

    # Map all of the TSRs that did not get mapped to search regions
    mapped_flow_through_tsrs_dict = truQuant.map_flow_through_tsrs(annotations_dict, flow_through_tsrs)

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



def blacklist_extended_gene_bodies(tsr_file, downstream_extension):
    """

    :param tsr_file:
    :type tsr_file: str
    :param blacklist_filename:
    :type blacklist_filename: str
    :param downstream_extension:
    :type downstream_extension: int
    :return:
    """
    annotation_extension = 1000
    percent_for_blacklisting = 0.3
    truQuant_regions_dict = {}

    blacklist_filename = generate_random_filename()

    # 1: Make the regions we are going to be searching for max TSSs in max TSRs
    search_regions_dict, annotations_dict = truQuant.make_search_regions(annotation_file, annotation_extension)

    # 2: Make the pause regions and gene bodies
    gene_tsr_dict, flow_through_tsrs = truQuant.map_tsrs_to_search_regions(tsr_file, search_regions_dict)
    max_tsrs_dict, non_max_tsrs_dict = truQuant.find_max_tsr_in_search_region(gene_tsr_dict)
    _define_pause_regions_and_gene_bodies(max_tsrs_dict, annotations_dict, truQuant_regions_dict, downstream_extension)
    _make_blacklisted_regions(blacklist_filename, annotations_dict, max_tsrs_dict, non_max_tsrs_dict, flow_through_tsrs,
                             percent_for_blacklisting)

    return blacklist_filename
