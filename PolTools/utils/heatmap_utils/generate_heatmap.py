import os
import sys

from PolTools.utils.constants import generate_heatmap_location
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.remove_files import remove_files
from PolTools.utils.heatmap_utils.set_matrix_bounds import set_matrix_bounds
from PolTools.utils.heatmap_utils.get_max_value_from_matrix import get_max_value_from_matrix
from PolTools.utils.heatmap_utils.get_min_value_from_matrix import get_min_value_from_matrix


class Ticks:
    def __init__(self, minor_tick_mark_interval_size, major_tick_mark_interval_size, offset=None, width=1):
        """
        :param minor_tick_mark_interval_size: Number of pixels per each minor tick mark
        :type minor_tick_mark_interval_size: float
        :param major_tick_mark_interval_size: Number of pixels per each major tick mark
        :type major_tick_mark_interval_size: float
        """

        if minor_tick_mark_interval_size <= 0:
            sys.stderr.write("The minor tick mark interval size must be greater than 0. It was set to " + str(
                minor_tick_mark_interval_size))
            sys.exit(1)

        if major_tick_mark_interval_size <= 0:
            sys.stderr.write("The major tick mark interval size must be greater than 0. It was set to " + str(
                major_tick_mark_interval_size))
            sys.exit(1)

        self.minor_tick_mark_interval_size = minor_tick_mark_interval_size
        self.major_tick_mark_interval_size = major_tick_mark_interval_size
        self.offset = offset
        self.width = 1

    def get_locations(self, matrix_width):
        """

        :param matrix_width:
        :return: minor_tick_locations, major_tick_locations
        """

        minor_tick_locations = []
        major_tick_locations = []
        minor_seen = set()
        major_seen = set()

        # I don't want tick marks on the first column so I add 1 to each of the seen sets
        minor_seen.add(0)
        major_seen.add(0)

        for i in range(1, matrix_width):
            if (i // self.minor_tick_mark_interval_size) not in minor_seen:
                minor_tick_locations.append(i-1)
                minor_seen.add(i // self.minor_tick_mark_interval_size)

            if (i // self.major_tick_mark_interval_size) not in major_seen:
                major_tick_locations.append(i-1)
                major_seen.add(i // self.major_tick_mark_interval_size)

        # Add the offset if necessary
        if self.offset:
            offset_minor_tick_locations = [val + self.offset for val in minor_tick_locations]
            offset_major_tick_locations = [val + self.offset for val in major_tick_locations]

            minor_tick_locations = offset_minor_tick_locations
            major_tick_locations = offset_major_tick_locations

        if self.width != 1:
            wider_minor_ticks = []
            wider_major_ticks = []

            for location in minor_tick_locations:
                for i in range(self.width):
                    wider_minor_ticks.append(location + i)

            for location in major_tick_locations:
                for i in range(self.width):
                    major_tick_locations.append(location + i)

        return minor_tick_locations, major_tick_locations


def make_ticks_matrix(matrix_width, matrix_height, max_black_value, ticks):
    # Now let's make the tick marks matrix
    row_list = [0 for _ in range(matrix_width)]
    ticks_matrix = [row_list.copy() for _ in range(matrix_height)]

    # Make the intervals of the major and minor tick marks
    minor_tick_locations, major_tick_locations = ticks.get_locations(matrix_width)

    for row in range(int(matrix_height / 4)):
        for col in minor_tick_locations:
            ticks_matrix[row][col] = max_black_value

    for row in range(int(3 * (matrix_height / 4))):
        for col in major_tick_locations:
            ticks_matrix[row][col] = max_black_value

    return ticks_matrix


def add_ticks_matrix(old_matrix_filename, ticks):
    matrix = []
    all_values = []
    with open(old_matrix_filename) as file:
        for line in file:
            if line != "":
                matrix.append([float(val) for val in line.split()])
                all_values.extend([float(val) for val in line.split()])

    max_value = max(all_values)
    matrix_width = len(matrix[0])
    matrix_height = len(matrix)

    if ticks:
        ticks_matrix = make_ticks_matrix(matrix_width, 50, max_value, ticks)

        # Add the ticks matrix to the bottom of the matrix
        matrix.extend(ticks_matrix)

    # Now write the new matrix to a file
    new_matrix_filename = generate_random_filename('.matrix')
    with open(new_matrix_filename, 'w') as file:
        for row in matrix:
            output_row = [str(val) for val in row]
            file.write("\t".join(output_row) + "\n")

    return new_matrix_filename


def generate_heatmap(matrix_filename, color_scheme, output_filename, gamma, min_value, max_value, center_value="default", ticks=None):
    """
    Will create a heatmap of dimensions width x height at the output_filename location

    :param center_value:
    :param min_value:
    :param max_value:
    :param ticks: Ticks object
    :param matrix_filename:
    :type matrix_filename: str
    :param color_scheme: either "gray" or "red/blue"
    :type color_scheme: str
    :param output_filename:
    :type output_filename: str
    :param gamma:
    :type gamma: float
    :return:
    """

    # Make sure the averaged_matrix_filename exists
    if not os.path.isfile(matrix_filename):
        sys.stderr.write("The filename " + matrix_filename + " was not found. Exiting ...\n")
        sys.exit(1)

    if color_scheme not in ["gray", "red/blue", "cyan/magenta"]:
        sys.stderr.write("The chosen color scheme (" + color_scheme + ") is not available. Exiting ...\n")
        sys.exit(1)

    if gamma < 0:
        sys.stderr.write("The gamma must be positive. It was set to " + str(gamma) + "\n")
        sys.exit(1)

    # Set the max and min values
    if min_value == None:
        min_value = get_min_value_from_matrix(matrix_filename)
    if max_value == None:
        max_value = get_max_value_from_matrix(matrix_filename)

    finalized_matrix = set_matrix_bounds(matrix_filename, min_value, max_value)

    if min_value < 0 and max_value > 0 and center_value == 'default':
        center_value = 0

    if center_value == "default":
        center_value = (max_value + min_value) / 2

    if ticks:
        new_matrix_filename = add_ticks_matrix(finalized_matrix, ticks)

        os.system("/usr/bin/Rscript " + generate_heatmap_location + " " +
                  " ".join([new_matrix_filename, color_scheme, output_filename, str(gamma), str(min_value),
                            str(max_value), str(center_value)]))

        remove_files(new_matrix_filename)

    else:
        os.system("/usr/bin/Rscript " + generate_heatmap_location + " " +
                  " ".join([finalized_matrix, color_scheme, output_filename, str(gamma), str(min_value),
                            str(max_value), str(center_value)]))

    remove_files(finalized_matrix)
