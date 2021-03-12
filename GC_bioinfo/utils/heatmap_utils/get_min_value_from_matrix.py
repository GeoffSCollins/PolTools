def get_min_value_from_matrix(matrix_filename):
    """
    Returns the minimum value of a matrix file

    :param matrix_filename: str
    :rtype: float
    """
    matrix = []
    with open(matrix_filename) as file:
        for line in file:
            matrix.extend([float(val) for val in line.rstrip().split()])

    return min(matrix)
