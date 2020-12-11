from collections import defaultdict

def build_counts_dict(sequencing_filename):
    """
    Builds two dictionaries for the 5' ends and 3' ends of transcripts. Dictionaries can be accessed like so:
    five_prime_counts_dict[chromosome][strand][five_prime_position], where five_prime_position is an integer

    :param sequencing_filename: filename of the sequencing data collected
    :type sequencing_filename: str
    :return: two dictionaries, five_prime_counts_dict and three_prime_counts_dict
    :rtype: (dict, dict)
    """

    five_prime_counts_dict = {}
    three_prime_counts_dict = {}

    with open(sequencing_filename) as file:
        for i, line in enumerate(file):
            chromosome, left, right, _, _, strand = line.rstrip().split()

            left = int(left)
            right = int(right)

            if chromosome not in five_prime_counts_dict:
                five_prime_counts_dict[chromosome] = {"+": defaultdict(lambda: 0), "-": defaultdict(lambda: 0)}

            if chromosome not in three_prime_counts_dict:
                three_prime_counts_dict[chromosome] = {"+": defaultdict(lambda: 0), "-": defaultdict(lambda: 0)}

            if strand == "+":
                five_prime_position = left
                three_prime_position = right - 1

            else:
                five_prime_position = right - 1
                three_prime_position = left

            if five_prime_position not in five_prime_counts_dict[chromosome][strand]:
                five_prime_counts_dict[chromosome][strand][five_prime_position] = 1
            else:
                five_prime_counts_dict[chromosome][strand][five_prime_position] += 1

            if three_prime_position not in three_prime_counts_dict[chromosome][strand]:
                three_prime_counts_dict[chromosome][strand][three_prime_position] = 1
            else:
                three_prime_counts_dict[chromosome][strand][three_prime_position] += 1

    return five_prime_counts_dict, three_prime_counts_dict
