from PolTools.utils.make_random_filename import generate_random_filename


def make_read_end_file(sequencing_filename, end):
    """
    Makes a new file containing only the 5' ends of the given file.

    :param sequencing_filename: filename of the sequencing data collected
    :type sequencing_filename: str
    :return: a new filename which contains the 5' ends of the given file
    :rtype: str
    """
    new_filename = generate_random_filename()

    with open(sequencing_filename) as file:
        with open(new_filename, 'w') as output_file:
            for line in file:
                chromosome, left, right, gene_name, score, strand = line.split()

                # Now get the 5' end
                if end == "five":
                    if strand == "+":
                        left = int(left)
                        right = int(left) + 1
                    else:
                        right = int(right)
                        left = int(right) - 1

                    # Output the 5' end to the file
                    output_file.write("\t".join([chromosome, str(left), str(right), gene_name, score, strand]) + "\n")

                elif end == "three":
                    if strand == "+":
                        left = int(right) - 1
                        right = int(right)
                    else:
                        left = int(left)
                        right = int(left) + 1

                    output_file.write("\t".join([chromosome, str(left), str(right), gene_name, score, strand]) + "\n")

                else:
                    raise ValueError("Read end must be either five or three.")

    return new_filename
