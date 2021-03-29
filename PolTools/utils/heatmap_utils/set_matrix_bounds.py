from PolTools.utils.make_random_filename import generate_random_filename


def set_matrix_bounds(matrix_filename, minimum_value, maximum_value):
    """
    Sets minimum and maximum values for the matrix by replacing values smaller than the minimum value with the minimum
    value and replacing values larger than the maximum value with the maximum value

    :param matrix_filename:
    :type matrix_filename: str
    :param minimum_value:
    :type minimum_value: float
    :param maximum_value:
    :type maximum_value: float
    :return: new matrix filename
    :rtype: str

    """

    original_matrix = []
    with open(matrix_filename) as file:
        for line in file:
            original_matrix.append([float(val) for val in line.rstrip().split()])

    width = len(original_matrix[0])
    height = len(original_matrix)

    new_matrix = []
    for row in range(height):
        new_matrix.append([0 for _ in range(width)])

    for row in range(height):
        for col in range(width):
            curr_value = original_matrix[row][col]

            if curr_value > maximum_value:
                new_matrix[row][col] = maximum_value
            elif curr_value < minimum_value:
                new_matrix[row][col] = minimum_value
            else:
                new_matrix[row][col] = curr_value

    new_matrix_filename = generate_random_filename()

    with open(new_matrix_filename, 'w') as file:
        for row in new_matrix:
            str_row = [str(val) for val in row]
            file.write("\t".join(str_row) + "\n")

    return new_matrix_filename
