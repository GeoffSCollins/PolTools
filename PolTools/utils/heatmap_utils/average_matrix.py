from PolTools.utils.make_random_filename import generate_random_filename

def average_matrix(filename, num_lines_to_average):
    """
    Vertically averages a matrix by num_lines_to_average

    :param filename:
    :type filename str
    :param num_lines_to_average:
    :type num_lines_to_average: int
    :return: filename of averaged matrix
    :rtype: str
    """

    with open(filename) as file:
        lines = [line.rstrip() for line in file if line.rstrip()]

    averaged_matrix_filename = generate_random_filename()

    with open(averaged_matrix_filename, 'w') as file:
        curr_sums = [0] * len(lines[0].split())
        for i, line in enumerate(lines):
            vals = line.split()

            for j, val in enumerate(vals):
                curr_sums[j] += float(val)

            if (i+1) % num_lines_to_average == 0:
                # Average and print out
                avgs = [str(curr_sum / num_lines_to_average) for curr_sum in curr_sums]

                file.write("\t".join(avgs) + "\n")
                curr_sums = [0] * len(curr_sums)

    return averaged_matrix_filename
