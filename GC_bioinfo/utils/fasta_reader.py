
class NotFastaFileError(Exception):

    def __init__(self, filename):
        self.filename = filename

    def __str__(self):
        return "File " + self.filename + " is not a fasta file"


def read_fasta(filename):
    """
    Returns the sequences from a fasta file.

    :param filename: filename of the fasta file
    :type filename: str
    :return: list
    """
    sequences = []
    with open(filename) as file:
        for i, line in enumerate(file):
            if i % 2 == 0:
                # This line should start with a >
                if line[0] != ">":
                    raise NotFastaFileError(filename)

            else:
                sequences.append(line.rstrip())

    return sequences
