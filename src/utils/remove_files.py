import os


def remove_files(*files):
    """
    Deletes the files from the system.

    :param str or list files: files to be deleted
    """
    command = "rm "
    for file in files[:-1]:
        if isinstance(file, list):
            for filename in file:
                command += filename + " "
        else:
            command += file + " "

    if isinstance(files[-1], list):
        for filename in files[-1]:
            command += filename + " "
    else:
        command += files[-1]

    os.system(command)
