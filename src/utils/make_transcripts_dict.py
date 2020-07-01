'''
Defines a function build_transcripts_dict which takes a filename and returns a dictionary for all the transcripts
'''

def build_transcripts_dict(sequencing_filename):
    """
    Builds a dictionary containing all of the transcripts from the given file. Can be accessed like so:
    transcripts_dict[chromosome]["+"][five_prime_position], to get a list of the connected 3' ends

    :param sequencing_filename: filename of the sequencing data collected
    :type sequencing_filename: str
    :return: a dictionary containing the connection between all 5' ends and 3' ends
    :rtype: dict
    """
    transcripts_dict = {}
    with open(sequencing_filename) as file:
        for i, line in enumerate(file):
            chromosome, left, right, _, _, strand = line.rstrip().split()

            left = int(left)
            right = int(right)

            if chromosome not in transcripts_dict:
                transcripts_dict[chromosome] = {"+": {}, "-": {}}

            if strand == "+":
                if left not in transcripts_dict[chromosome]["+"]:
                    transcripts_dict[chromosome]["+"][left] = []

                transcripts_dict[chromosome][strand][left].append(right)
            else:
                if right not in transcripts_dict[chromosome]["-"]:
                    transcripts_dict[chromosome]["-"][right] = []

                transcripts_dict[chromosome][strand][right].append(left)

    return transcripts_dict