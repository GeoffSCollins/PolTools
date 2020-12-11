import csv

from GC_bioinfo.utils.make_random_filename import generate_random_filename

def make_five_bed_file(sequencing_filename):
    """
    Makes a new file containing only the 5' ends of the given file.

    :param sequencing_filename: filename of the sequencing data collected
    :type sequencing_filename: str
    :return: a new filename which contains the 5' ends of the given file
    :rtype: str
    """
    output_lines = []

    with open(sequencing_filename) as file:
        for line in file:
            chromosome, left, right, gene_name, score, strand = line.split()

            # Now get the 5' end
            if strand == "+":
                left = int(left)
                right = int(left) + 1
            else:
                right = int(right)
                left = int(right) - 1

            # Output the 5' end to the file
            output_lines.append([chromosome, left, right, gene_name, score, strand])


    new_filename = generate_random_filename()

    with open(new_filename, 'w') as output_file:
        output_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

        for line in output_lines:
            output_writer.writerow(line)

    return new_filename