def determine_region_length(regions_filename):
    """
    Determines the length of regions using the first region in the file.

    :param regions_filename: filename of the regions of which to get the length
    :return: length of the first region in the regions file
    :rtype: int
    """
    with open(regions_filename) as file:
        lines = file.readlines()

    chromosome, left, right, gene_name, _, strand = lines[0]

    return int(right) - int(left)
