import os


def remove_files(*files):
    """
    Deletes the files from the system.

    :param files: files to be deleted
    :type files: str or list or tuple
    """
    for file in files:
        if isinstance(file, list) or isinstance(file, tuple):
            for file_in_list in file:
                remove_files(file_in_list)
        else:
            os.remove(file)
