from PolTools.utils.verify_bed_file import verify_bed_files

def determine_region_length(regions_filename):
    """
    Determines the length of regions using the first region in the bed file.

    :param regions_filename: filename of the bed formatted regions of which to get the length
    :return: length of the first region in the regions file
    :rtype: int
    """

    verify_bed_files(regions_filename)

    with open(regions_filename) as file:
        first_line = file.readline()

    chromosome, left, right, gene_name, _, strand = first_line.split()

    return int(right) - int(left)
