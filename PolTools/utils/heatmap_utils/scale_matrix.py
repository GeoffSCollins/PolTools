
from PolTools.utils.make_random_filename import generate_random_filename

def scale_matrix(original_matrix_filename, scale_factor):
    """
    Multiplies each value in the matrix by the scale factor

    :param original_matrix_filename:
    :type original_matrix_filename: str
    :param scale_factor:
    :type scale_factor: float
    :return:
    """
    matrix = []
    with open(original_matrix_filename) as file:
        for line in file:
            if line != "":
                matrix.append([float(val) * scale_factor for val in line.split()])

    normalized_matrix_filename = generate_random_filename('.matrix')
    with open(normalized_matrix_filename, 'w') as file:
        for row in matrix:
            output_row = [str(val) for val in row]
            file.write("\t".join(output_row) + "\n")

    return normalized_matrix_filename