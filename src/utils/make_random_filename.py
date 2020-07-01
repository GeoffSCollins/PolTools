import random
import string
import os.path


def generate_random_filename():
    """
    Returns a unique filename in the /tmp directory

    :return: a unique filename in the /tmp directory
    :rtype: str
    """
    random_filename = "/tmp/" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=8)) + ".bed"

    while os.path.isfile(random_filename):
        # Get random filenames until one is unique
        random_filename = "/tmp/" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=8)) + ".bed"

    return random_filename
