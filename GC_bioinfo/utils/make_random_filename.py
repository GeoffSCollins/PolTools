import random
import string
import os.path


def generate_random_filename(extension=".bed"):
    """
    Returns a unique filename in the /tmp directory

    :return: a unique filename in the /tmp directory
    :rtype: str
    """
    random_filename = "/tmp/GC_bioinfo_" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=8)) + extension

    while os.path.isfile(random_filename):
        # Get random filenames until one is unique
        random_filename = "/tmp/GC_bioinfo_" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=8)) + extension

    return random_filename
