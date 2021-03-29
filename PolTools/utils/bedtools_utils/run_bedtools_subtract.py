import os

from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.verify_bed_file import verify_bed_files


def run_subtract(sequencing_file, *subtracted_regions, output_filename=''):
    """
    Runs strand specific bedtools subtract to remove regions from a sequencing file.

    :param sequencing_filename: filename of the sequencing data collected
    :type sequencing_filename: str
    :param subtracted_regions: filename of the regions to be subtracted (can be multiple)
    :type subtracted_regions: str
    :param output_filename: optional name of the output file (will be random if not provided)
    :type output_filename: str
    :return: filename of the resultant bedtools subtract output
    :rtype: str
    """
    if output_filename == '':
        output_filename = generate_random_filename()

    verify_bed_files(sequencing_file, list(subtracted_regions))

    if len(subtracted_regions) == 1:
        os.system(
            "bedtools subtract -nonamecheck -s -A -a " + sequencing_file + " -b " + subtracted_regions[0] + " > " + output_filename)

    else:
        command = "bedtools subtract -nonamecheck -s -A -a " + sequencing_file

        for region_filename in subtracted_regions[:-1]:
            command += " -b " + region_filename + " | bedtools subtract -nonamecheck -s -A -a stdin -b "

        command += subtracted_regions[-1] + " > " + output_filename

        os.system(command)

    return output_filename
