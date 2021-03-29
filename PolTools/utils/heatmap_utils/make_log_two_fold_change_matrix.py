import math
import sys

from PolTools.utils.make_random_filename import generate_random_filename


def make_log_two_fold_change_matrix(numerator_filename, denominator_filename):
    # Going to divide the first by the second
    first_matrix = []
    with open(numerator_filename) as file:
        for line in file:
            first_matrix.append([float(val) for val in line.rstrip().split()])

    second_matrix = []
    with open(denominator_filename) as file:
        for line in file:
            second_matrix.append([float(val) for val in line.rstrip().split()])

    # Verify the matricies are the same size
    first_matrix_width = len(first_matrix[0])
    first_matrix_height = len(first_matrix)

    second_matrix_width = len(second_matrix[0])
    second_matrix_height = len(second_matrix)

    if not (first_matrix_width == second_matrix_width) and not (first_matrix_height == second_matrix_height):
        sys.stderr.write("The matricies are not the same size. Exiting ...")
        sys.exit(1)

    # Now divide the first by the second
    fold_change_matrix = []
    for i in range(first_matrix_height):
        fold_change_matrix.append([0 for _ in range(first_matrix_width)])

    for row in range(first_matrix_height):
        for col in range(first_matrix_width):
            numerator = first_matrix[row][col]
            denominator = second_matrix[row][col]

            if numerator == 0:
                numerator = 1

            if denominator == 0:
                denominator = 1

            # Take the log2 of this value
            log_two_fold_change = math.log2(numerator / denominator)

            fold_change_matrix[row][col] = log_two_fold_change

    # Write to a file
    fold_change_matrix_filename = generate_random_filename(".matrix")
    with open(fold_change_matrix_filename, 'w') as file:
        for row in fold_change_matrix:
            str_row = [str(val) for val in row]
            file.write("\t".join(str_row) + "\n")

    return fold_change_matrix_filename
