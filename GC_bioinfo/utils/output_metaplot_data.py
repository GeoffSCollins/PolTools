from itertools import chain

from GC_bioinfo.utils.print_tab_delimited import print_tab_delimited


def output_metaplot_data(averages, region_length, prime_name):
    """

    :param averages: averages list from the metaplots programs
    :type averages: list
    :param region_length: length of the region
    :type region_length: int
    :param prime_name: either "five prime" or "three prime"
    :type prime_name: str
    :return:
    """
    avgs_data, files = [x for x in zip(*averages)]

    # Merge all of the lists together
    merged_list = [list(chain.from_iterable(x)) for x in zip(*avgs_data)]

    # 5. Put the data into a file
    header = ["Position"]
    # Write the header first
    for file in files:
        if prime_name:
            # Include a space before print the prime name
            header.append(file.split("/")[-1] + " " + prime_name + " sense strand")
            header.append(file.split("/")[-1] + " " + prime_name + " divergent strand")
        else:
            header.append(file.split("/")[-1] + " sense strand")
            header.append(file.split("/")[-1] + " divergent strand")

    print_tab_delimited(header)

    for i, base_list in enumerate(merged_list):
        position = i - (region_length / 2)
        if position >= 0:
            position += 1

        print_tab_delimited([position] + base_list)