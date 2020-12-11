def average_vertically(input_2d_list):
    """
    Takes in a 2d list and gets the average of each x position

    :param input_2d_list:
    :return: a list of the averaged values
    :rtype: list
    """
    height = len(input_2d_list)

    if height == 0:
        return []

    width = len(input_2d_list[0])

    averages_list = [0] * width

    # We loop through each base in the region
    for curr_position in range(width):
        current_sum = 0

        # Loop through all of the regions
        for region in input_2d_list:
            current_sum += region[curr_position]

        avg = current_sum / height
        averages_list[curr_position] = avg

    return averages_list