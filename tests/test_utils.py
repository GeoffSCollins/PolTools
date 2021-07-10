import unittest.mock

from PolTools.utils.heatmap_utils.add_matrices import add_matrices
from PolTools.utils.heatmap_utils.average_matrix import average_matrix
from PolTools.utils.get_region_length import determine_region_length
from PolTools.utils.build_counts_dict import build_counts_dict
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.make_read_end_file import make_read_end_file
from PolTools.utils.make_transcripts_dict import build_transcripts_dict
from PolTools.utils.bedtools_utils.run_bedtools_coverage import run_coverage
from PolTools.utils.bedtools_utils.run_bedtools_getfasta import run_getfasta
from PolTools.utils.bedtools_utils.run_bedtools_subtract import run_subtract
from PolTools.utils.heatmap_utils.scale_matrix import scale_matrix
from PolTools.utils.heatmap_utils.set_matrix_bounds import set_matrix_bounds
from PolTools.utils.verify_bed_file import verify_bed_files
from PolTools.utils.verify_region_length_is_even import verify_region_length_is_even

from PolTools.utils.remove_files import remove_files

from quieter import Quieter


class TestAddMatrices(unittest.TestCase):

    def test_different_sized_matrices(self):
        # Make example 2x2 matrix
        two = generate_random_filename()
        with open(two, 'w') as file:
            file.write("1 1\n1 1")

        # Make an example 3x3 matrix
        three = generate_random_filename()
        with open(three, 'w') as file:
            file.write("2 2 2\n2 2 2")

        with self.assertRaises(SystemExit):
            with Quieter():
                add_matrices([two, three])

        remove_files(two, three)

    def test_file_not_matrix(self):
        not_matrix = generate_random_filename()
        with open(not_matrix, 'w') as file:
            file.write("asfd asd fasdf ")

        with self.assertRaises(ValueError):
            add_matrices([not_matrix])

        remove_files(not_matrix)

    def test_add_two(self):
        mat_one = generate_random_filename()
        with open(mat_one, 'w') as file:
            file.write("1 1\n1 1")

        mat_two = generate_random_filename()
        with open(mat_two, 'w') as file:
            file.write("1 1\n1 1")

        result_mat = add_matrices([mat_one, mat_two])

        with open(result_mat) as file:
            r = []
            for line in file:
                r.append([float(val) for val in line.split()])

        self.assertEqual(r, [[2, 2], [2, 2]])

        remove_files(mat_one, mat_two, result_mat)


class TestAverageMatrix(unittest.TestCase):

    def get_test_mat(self):
        mat = generate_random_filename()
        with open(mat, 'w') as file:
            file.write("""1 1
                    2 1
                    3 1
                    4 1
                    5 1
                    6 2
                    7 2
                    8 2
                    9 2
                    10 2
                    """)
        return mat

    def test_file_not_matrix(self):
        not_matrix = generate_random_filename()
        with open(not_matrix, 'w') as file:
            file.write("asfd asd fasdf ")

        with self.assertRaises(ValueError):
            add_matrices([not_matrix])

        remove_files(not_matrix)

    def test_even_average(self):
        test_mat = self.get_test_mat()
        result_mat = average_matrix(test_mat, 5)

        with open(result_mat) as file:
            r = []
            for line in file:
                r.append([float(val) for val in line.split()])

        self.assertEqual(r, [[3, 1], [8, 2]])

        remove_files(test_mat, result_mat)

    def test_odd_average(self):
        test_mat = self.get_test_mat()
        result_mat = average_matrix(test_mat, 4)

        with open(result_mat) as file:
            r = []
            for line in file:
                r.append([float(val) for val in line.split()])

        # The last line is discarded because it cannot be evenly divided
        self.assertEqual(r, [[2.5, 1], [6.5, 1.75]])

        remove_files(test_mat, result_mat)


class TestBlacklistExtendedGeneBodies(unittest.TestCase):
    # Todo
    pass


class TestDetermineRegionLength(unittest.TestCase):

    def test_invalid_file(self):
        invalid_filename = generate_random_filename()

        with open(invalid_filename, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            determine_region_length(invalid_filename)

        remove_files(invalid_filename)

    def test_one_region(self):
        one_region_filename = generate_random_filename()

        with open(one_region_filename, 'w') as file:
            file.write("chr1    1   11  name    score   +\n")

        self.assertEqual(determine_region_length(one_region_filename), 10)

        remove_files(one_region_filename)

    def test_multiple_regions(self):
        multiple_regions_filename = generate_random_filename()

        with open(multiple_regions_filename, 'w') as file:
            file.write("chr1    1   11  name    score   +\n")
            file.write("chr1    11   21  name    score   +\n")
            file.write("chr1    21   31  name    score   +\n")
            file.write("chr1    31   41  name    score   +\n")

        self.assertEqual(determine_region_length(multiple_regions_filename), 10)

        remove_files(multiple_regions_filename)



class TestMakeFiveAndThreeDict(unittest.TestCase):

    def test_invalid_file(self):
        invalid_filename = generate_random_filename()

        with open(invalid_filename, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            build_counts_dict(invalid_filename, "five")

        remove_files(invalid_filename)

    def test_positive_strand_one_read(self):
        pos_strand_filename = generate_random_filename()

        with open(pos_strand_filename, 'w') as file:
            file.write("chr1    1   11  name    score   +\n")

        five_dict = build_counts_dict(pos_strand_filename, "five")

        self.assertEqual(five_dict["chr1"]["+"][1], 1)

        three_dict = build_counts_dict(pos_strand_filename, "three")
        self.assertEqual(three_dict["chr1"]["+"][10], 1)

        remove_files(pos_strand_filename)

    def test_positive_strand_multiple_reads(self):
        pos_strand_filename = generate_random_filename()

        with open(pos_strand_filename, 'w') as file:
            file.write("chr1    1   11  name    score   +\n")
            file.write("chr1    1   7  name    score   +\n")
            file.write("chr1    1   99  name    score   +\n")
            file.write("chr1    100   789  name    score   +\n")
            file.write("chr1    51   789  name    score   +\n")

        five_dict = build_counts_dict(pos_strand_filename, "five")

        self.assertEqual(five_dict["chr1"]["+"][1], 3)
        self.assertEqual(five_dict["chr1"]["+"][100], 1)
        self.assertEqual(five_dict["chr1"]["+"][51], 1)

        three_dict = build_counts_dict(pos_strand_filename, "three")

        self.assertEqual(three_dict["chr1"]["+"][10], 1)
        self.assertEqual(three_dict["chr1"]["+"][6], 1)
        self.assertEqual(three_dict["chr1"]["+"][98], 1)
        self.assertEqual(three_dict["chr1"]["+"][788], 2)

        remove_files(pos_strand_filename)


    def test_negative_strand_one_read(self):
        neg_strand_filename = generate_random_filename()

        with open(neg_strand_filename, 'w') as file:
            file.write("chr1    1   11  name    score   -\n")

        five_dict = build_counts_dict(neg_strand_filename, "five")

        self.assertEqual(five_dict["chr1"]["-"][10], 1)

        three_dict = build_counts_dict(neg_strand_filename, "three")
        self.assertEqual(three_dict["chr1"]["-"][1], 1)

        remove_files(neg_strand_filename)

    def test_negative_strand_multiple_reads(self):
        neg_strand_filename = generate_random_filename()

        with open(neg_strand_filename, 'w') as file:
            file.write("chr1    1   11  name    score   -\n")
            file.write("chr1    1   7  name    score   -\n")
            file.write("chr1    1   99  name    score   -\n")
            file.write("chr1    100   789  name    score   -\n")
            file.write("chr1    51   789  name    score   -\n")

        five_dict = build_counts_dict(neg_strand_filename, "five")

        self.assertEqual(five_dict["chr1"]["-"][10], 1)
        self.assertEqual(five_dict["chr1"]["-"][6], 1)
        self.assertEqual(five_dict["chr1"]["-"][98], 1)
        self.assertEqual(five_dict["chr1"]["-"][788], 2)

        three_dict = build_counts_dict(neg_strand_filename, "three")

        self.assertEqual(three_dict["chr1"]["-"][1], 3)
        self.assertEqual(three_dict["chr1"]["-"][100], 1)
        self.assertEqual(three_dict["chr1"]["-"][51], 1)

        remove_files(neg_strand_filename)


class TestMakeFivePrimeBedFile(unittest.TestCase):

    def test_invalid_file(self):
        invalid_filename = generate_random_filename()

        with open(invalid_filename, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            make_read_end_file(invalid_filename, 'five')

        remove_files(invalid_filename)

    def test_one_read_positive_strand(self):
        one_read_filename = generate_random_filename()

        with open(one_read_filename, 'w') as file:
            file.write("chr1    1   11  name    score   +\n")

        five_prime_file = make_read_end_file(one_read_filename, 'five')

        with open(five_prime_file) as file:
            results = file.readline().split()

        self.assertEqual(results, ["chr1", "1", "2", "name", "score", "+"])

        remove_files(one_read_filename, five_prime_file)

    def test_positive_already_five_prime(self):
        one_read_filename = generate_random_filename()

        with open(one_read_filename, 'w') as file:
            file.write("chr1    1   2  name    score   +\n")

        five_prime_file = make_read_end_file(one_read_filename, 'five')

        with open(five_prime_file) as file:
            results = file.readline().split()

        self.assertEqual(results, ["chr1", "1", "2", "name", "score", "+"])

        remove_files(one_read_filename, five_prime_file)

    def test_one_read_negative_strand(self):
        one_read_filename = generate_random_filename()

        with open(one_read_filename, 'w') as file:
            file.write("chr1    1   11  name    score   -\n")

        five_prime_file = make_read_end_file(one_read_filename, 'five')

        with open(five_prime_file) as file:
            results = file.readline().split()

        self.assertEqual(results, ["chr1", "10", "11", "name", "score", "-"])

        remove_files(one_read_filename, five_prime_file)

    def test_negative_alread_five_prime(self):
        one_read_filename = generate_random_filename()

        with open(one_read_filename, 'w') as file:
            file.write("chr1    1   2  name    score   -\n")

        five_prime_file = make_read_end_file(one_read_filename, 'five')

        with open(five_prime_file) as file:
            results = file.readline().split()

        self.assertEqual(results, ["chr1", "1", "2", "name", "score", "-"])

        remove_files(one_read_filename, five_prime_file)

    def test_two_reads_opposite_strands(self):
        filename = generate_random_filename()

        with open(filename, 'w') as file:
            file.write("chr1    1   11  name    score   -\n")
            file.write("chr1    21   50  name    score   +\n")

        five_prime_file = make_read_end_file(filename, 'five')

        with open(five_prime_file) as file:
            results = [line.rstrip().split() for line in file.readlines()]

        expected = [
            ["chr1", "10", "11", "name", "score", "-"],
            ["chr1", "21", "22", "name", "score", "+"]
        ]

        self.assertEqual(results, expected)

        remove_files(filename, five_prime_file)


class TestMakeThreePrimeBedFile(unittest.TestCase):

    def test_invalid_file(self):
        invalid_filename = generate_random_filename()

        with open(invalid_filename, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            make_read_end_file(invalid_filename, "three")

        remove_files(invalid_filename)

    def test_one_read_positive_strand(self):
        one_read_filename = generate_random_filename()

        with open(one_read_filename, 'w') as file:
            file.write("chr1    1   11  name    score   +\n")

        three_prime_file = make_read_end_file(one_read_filename, "three")

        with open(three_prime_file) as file:
            results = file.readline().split()

        self.assertEqual(results, ["chr1", "10", "11", "name", "score", "+"])

        remove_files(one_read_filename, three_prime_file)

    def test_one_read_negative_strand(self):
        one_read_filename = generate_random_filename()

        with open(one_read_filename, 'w') as file:
            file.write("chr1    1   11  name    score   -\n")

        three_prime_file = make_read_end_file(one_read_filename, "three")

        with open(three_prime_file) as file:
            results = file.readline().split()

        self.assertEqual(results, ["chr1", "1", "2", "name", "score", "-"])

        remove_files(one_read_filename, three_prime_file)

    def test_two_reads_opposite_strands(self):
        filename = generate_random_filename()

        with open(filename, 'w') as file:
            file.write("chr1    1   11  name    score   -\n")
            file.write("chr1    21   50  name    score   +\n")

        three_prime_file = make_read_end_file(filename, "three")

        with open(three_prime_file) as file:
            results = [line.rstrip().split() for line in file.readlines()]

        expected = [
            ["chr1", "1", "2", "name", "score", "-"],
            ["chr1", "49", "50", "name", "score", "+"]
        ]

        self.assertEqual(results, expected)

        remove_files(filename, three_prime_file)


class TestMakeTranscriptsDict(unittest.TestCase):

    def test_invalid_file(self):
        invalid_filename = generate_random_filename()

        with open(invalid_filename, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            build_transcripts_dict(invalid_filename)

        remove_files(invalid_filename)

    def test_one_read_positive_strand(self):
        one_read_filename = generate_random_filename()

        with open(one_read_filename, 'w') as file:
            file.write("chr1    1   11  name    score   +\n")

        result = build_transcripts_dict(one_read_filename)

        self.assertEqual(result["chr1"]["+"][1][10], 1)

        remove_files(one_read_filename)

    def test_one_read_negative_strand(self):
        one_read_filename = generate_random_filename()

        with open(one_read_filename, 'w') as file:
            file.write("chr1    1   11  name    score   -\n")

        result = build_transcripts_dict(one_read_filename)

        self.assertEqual(result["chr1"]["-"][10][1], 1)

        remove_files(one_read_filename)

    def test_two_reads_opposite_strands(self):
        filename = generate_random_filename()

        with open(filename, 'w') as file:
            file.write("chr1    1   11  name    score   -\n")
            file.write("chr1    21   50  name    score   +\n")

        result = build_transcripts_dict(filename)

        self.assertEqual(result["chr1"]["-"][10][1], 1)
        self.assertEqual(result["chr1"]["+"][21][49], 1)

        remove_files(filename)


class TestRunBedtoolsCoverage(unittest.TestCase):

    def test_invalid_file(self):
        invalid_regions_filename = generate_random_filename()
        invalid_sequencing_filename = generate_random_filename()

        with open(invalid_regions_filename, 'w') as file:
            file.write("Invalid file")

        with open(invalid_sequencing_filename, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            run_coverage(invalid_regions_filename, invalid_sequencing_filename)

        remove_files(invalid_regions_filename, invalid_sequencing_filename)

    def test_correct_usage(self):
        regions_filename = generate_random_filename()
        sequencing_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write("chr1\t1\t10\tregion\t0\t+")

        with open(sequencing_filename, 'w') as file:
            file.write("chr1\t3\t7\tread\t0\t+")

        output_filename = run_coverage(regions_filename, sequencing_filename)

        with open(output_filename) as file:
            result = file.readline().split()

        expected = ["chr1", "1", "10", "region", "0", "+", "1", "4", "9", "0.4444444"]

        self.assertEqual(result, expected)

        remove_files(regions_filename, sequencing_filename, output_filename)



    def test_correct_usage_with_flag(self):
        regions_filename = generate_random_filename()
        sequencing_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write("chr1\t1\t10\tregion\t0\t+")

        with open(sequencing_filename, 'w') as file:
            file.write("chr1\t3\t7\tread\t0\t+")

        output_filename = run_coverage(regions_filename, sequencing_filename, flags=["-d"])

        with open(output_filename) as file:
            result = [line.split() for line in file.readlines() if line]

        expected = [
            ["chr1", "1", "10", "region", "0", "+", "1", "0"],
            ["chr1", "1", "10", "region", "0", "+", "2", "0"],
            ["chr1", "1", "10", "region", "0", "+", "3", "1"],
            ["chr1", "1", "10", "region", "0", "+", "4", "1"],
            ["chr1", "1", "10", "region", "0", "+", "5", "1"],
            ["chr1", "1", "10", "region", "0", "+", "6", "1"],
            ["chr1", "1", "10", "region", "0", "+", "7", "0"],
            ["chr1", "1", "10", "region", "0", "+", "8", "0"],
            ["chr1", "1", "10", "region", "0", "+", "9", "0"]
            ]

        self.assertEqual(result, expected)

        remove_files(regions_filename, sequencing_filename, output_filename)


class TestRunBedtoolsSubtract(unittest.TestCase):

    def test_invalid_file(self):
        invalid_subtracted_region = generate_random_filename()
        invalid_sequencing_filename = generate_random_filename()

        with open(invalid_subtracted_region, 'w') as file:
            file.write("Invalid file")

        with open(invalid_sequencing_filename, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            run_subtract(invalid_subtracted_region, invalid_sequencing_filename)

        remove_files(invalid_subtracted_region, invalid_sequencing_filename)

    def test_proper_removal(self):
        regions_filename = generate_random_filename()
        sequencing_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write("chr1\t1\t10\tregion\t0\t+")

        with open(sequencing_filename, 'w') as file:
            file.write("chr1\t3\t7\tread\t0\t+")

        output_filename = run_subtract(regions_filename, sequencing_filename)

        with open(output_filename) as file:
            result = file.readline()

        self.assertEqual(result, "")

        remove_files(regions_filename, sequencing_filename, output_filename)

    def test_correct_usage_no_removal(self):
        subtracted_region_one = generate_random_filename()
        sequencing_filename = generate_random_filename()

        with open(subtracted_region_one, 'w') as file:
            file.write("chr1\t100\t110\tregion\t0\t+")

        with open(sequencing_filename, 'w') as file:
            file.write("chr1\t3\t7\tread\t0\t+")

        output_filename = run_subtract(sequencing_filename, subtracted_region_one)

        with open(output_filename) as file:
            result = file.readline().split()

        expected = ["chr1", "3", "7", "read", "0", "+"]

        self.assertEqual(result, expected)

        remove_files(subtracted_region_one, sequencing_filename, output_filename)


class TestRunBedtoolsGetFasta(unittest.TestCase):

    def test_invalid_file(self):
        invalid_region_filename = generate_random_filename()

        with open(invalid_region_filename, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            run_getfasta(invalid_region_filename)

        remove_files(invalid_region_filename)

    def test_valid_request(self):
        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write("chr1\t1\t10\tregion\t0\t+")

        output_filename = run_getfasta(regions_filename)

        with open(output_filename) as file:
            _ = file.readline()
            result = file.readline().rstrip()

        self.assertEqual(result, "N"*9)

        remove_files(regions_filename, output_filename)


class TestScaleMatrix(unittest.TestCase):

    def test_invalid_file(self):
        invalid_matrix_filename = generate_random_filename()

        with open(invalid_matrix_filename, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            scale_matrix(invalid_matrix_filename, 1)

        remove_files(invalid_matrix_filename)

    def test_proper_scale(self):
        mat_filename = generate_random_filename()

        with open(mat_filename, 'w') as file:
            file.write("1 2\n3 4")

        scaled_matrix = scale_matrix(mat_filename, 2)

        with open(scaled_matrix) as file:
            result = result = [[float(val) for val in line.split()] for line in file if line]

        expected = [
            [2, 4],
            [6, 8]
        ]

        self.assertEqual(result, expected)

        remove_files(mat_filename, scaled_matrix)

    def test_negative_scale(self):
        mat_filename = generate_random_filename()

        with open(mat_filename, 'w') as file:
            file.write("1 2\n3 4")

        scaled_matrix = scale_matrix(mat_filename, -2)

        with open(scaled_matrix) as file:
            result = [[float(val) for val in line.split()] for line in file if line]

        expected = [
            [-2, -4],
            [-6, -8]
        ]

        self.assertEqual(result, expected)

        remove_files(mat_filename, scaled_matrix)

    def test_invalid_scale_factor(self):
        mat_filename = generate_random_filename()

        with open(mat_filename, 'w') as file:
            file.write("1 2\n3 4")

        with self.assertRaises(Exception):
            scale_matrix(mat_filename, "invalid")

        remove_files(mat_filename)


class TestSetMatrixBounds(unittest.TestCase):

    def test_invalid_file(self):
        invalid_matrix_filename = generate_random_filename()

        with open(invalid_matrix_filename, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            set_matrix_bounds(invalid_matrix_filename, 0, 1)

        remove_files(invalid_matrix_filename)

    def test_basic_bounds(self):
        mat_filename = generate_random_filename()

        with open(mat_filename, 'w') as file:
            file.write("-8 -1.25\n0.5 1.8")

        bounded_matrix = set_matrix_bounds(mat_filename, 0, 1)

        with open(bounded_matrix) as file:
            result = [[float(val) for val in line.split()] for line in file if line]

        expected = [
            [0, 0],
            [0.5, 1]
        ]

        self.assertEqual(result, expected)

        remove_files(mat_filename, bounded_matrix)

    def test_negative_bounds(self):
        mat_filename = generate_random_filename()

        with open(mat_filename, 'w') as file:
            file.write("-8 -1.25\n0.5 1.8")

        bounded_matrix = set_matrix_bounds(mat_filename, -5, -1)

        with open(bounded_matrix) as file:
            result = [[float(val) for val in line.split()] for line in file if line]

        expected = [
            [-5, -1.25],
            [-1, -1]
        ]

        self.assertEqual(result, expected)

        remove_files(mat_filename, bounded_matrix)


class TestVerifyBedFile(unittest.TestCase):

    def test_invalid_file(self):
        invalid_bed_file = generate_random_filename()

        with open(invalid_bed_file, 'w') as file:
            file.write("Invalid file")

        with self.assertRaises(Exception):
            verify_bed_files(invalid_bed_file)

        remove_files(invalid_bed_file)

    def test_valid_file(self):
        correct_bed_file = generate_random_filename()

        with open(correct_bed_file, 'w') as file:
            file.write("chr1\t1\t10\tregion\t0\t+")

        verify_bed_files(correct_bed_file)

        remove_files(correct_bed_file)


class TestVerifyRegionLengthIsEven(unittest.TestCase):
    def test_not_even(self):
        region_length = 11
        with self.assertRaises(Exception):
            verify_region_length_is_even(region_length)

    def test_even(self):
        region_length = 12
        verify_region_length_is_even(region_length)


if __name__ == '__main__':
    unittest.main()
