import sys

from PolTools.utils.make_random_filename import generate_random_filename


def add_matrices(list_of_matrices):
    # Temporarily set the width and the height of the output matrix to 0
    width = 0
    height = 0

    if not list_of_matrices:
        sys.stderr.write("At least one matrix must be provided to add_matrices.\n")
        sys.exit(1)

    for i, matrix_filename in enumerate(list_of_matrices):
        current_matrix = []
        with open(matrix_filename) as file:
            for line in file:
                if line:
                    current_matrix.append([float(val) for val in line.split()])

        if i == 0:
            # Get the width and height
            width = len(current_matrix[0])
            height = len(current_matrix)

            # Make a matrix with all 0's
            matrix = []

            for i in range(height):
                matrix.append([0 for _ in range(width)])

        # If the current size is not equal to the first matrix size, exit
        if len(current_matrix[0]) != width or len(current_matrix) != height:
            sys.stderr.write("To add matrices together, they must all be the same size")
            sys.exit(1)

        # Now add the values
        for row in range(height):
            for col in range(width):
                matrix[row][col] += current_matrix[row][col]

    added_matrix_filename = generate_random_filename('.matrix')
    with open(added_matrix_filename, 'w') as file:
        for row in matrix:
            string_row = [str(val) for val in row]

            file.write("\t".join(string_row) + "\n")

    return added_matrix_filename
