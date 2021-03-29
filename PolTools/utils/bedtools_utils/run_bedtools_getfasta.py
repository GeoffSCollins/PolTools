import os

from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.verify_bed_file import verify_bed_files
from PolTools.utils.constants import hg38_fasta_file


def run_getfasta(regions_file, output_filename=''):
    """
    Runs strand specific bedtools getfasta to get the sequence at the regions in the regions_file.

    :param regions_filename: filename of the regions of the genome to quantify
    :type regions_filename: str
    :param output_filename: optional name of the output file (will be random if not provided)
    :type output_filename: str
    :return: filename of the resultant bedtools fasta output
    :rtype: str
    """

    if output_filename == '':
        output_filename = generate_random_filename(extension=".fa")

    verify_bed_files(regions_file)

    os.system("bedtools getfasta -s -fi " + hg38_fasta_file + " -bed " + regions_file + " > " + output_filename)

    return output_filename
