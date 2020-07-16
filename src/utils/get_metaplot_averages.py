def get_averages(input_2d_list, region_length):
    """

    :param input_2d_list:
    :param region_length:
    :return:
    """
    # This will take in a 2d list containing the regions and counts at that specific base
    # and output a list of the averages at each base
    averages_list = [0] * region_length

    # We loop through each base in the region
    for current_base_position, counts in enumerate(input_2d_list[0]):
        current_sum = 0

        # The total is the number of regions
        total = len(input_2d_list)

        # Loop through all of the regions
        for region in input_2d_list:
            current_sum += region[current_base_position]

        avg = current_sum / total
        averages_list[current_base_position] = avg

    return averages_list