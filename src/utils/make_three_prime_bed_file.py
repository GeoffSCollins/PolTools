import csv

def make_three_bed_file(sequencing_filename):
    """
    Makes a new file containing only the 3' ends of the given file.

    :param sequencing_filename: filename of the sequencing data collected
    :type sequencing_filename: str
    :return: a new filename which contains the 3' ends of the given file
    :rtype: str
    """
    with open(sequencing_filename) as file:
        extension = sequencing_filename.split(".")[-1]
        if extension == sequencing_filename:
            new_filename = sequencing_filename + "-3-end.bed"
        else:
            new_filename = sequencing_filename.replace(extension, "3-end." + extension)

        with open(new_filename, 'w') as output_file:
            output_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

            for line in file:
                chromosome, left, right, gene_name, score, strand = line.split()

                # Now get the 3' end
                if strand == "+":
                    left = int(right) - 1
                    right = int(right)
                else:
                    left = int(left)
                    right = int(left) + 1

                # Output the 5' end to the file
                output_writer.writerow([chromosome, left, right, gene_name, score, strand])

    return new_filename