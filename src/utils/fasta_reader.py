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
            if i % 2 == 1:
                # This line has the sequence
                sequences.append(line.rstrip())

    return sequences
