import unittest.mock

from PolTools.main_programs import nucleotide_heatmap
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.remove_files import remove_files

from quieter import Quieter

class TestNucleotideHeatmap(unittest.TestCase):

    def test_incorrect_number_of_arguments(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                nucleotide_heatmap.parse_args([])

        with self.assertRaises(SystemExit):
            with Quieter():
                nucleotide_heatmap.parse_args(["max_tss_file"])

        with self.assertRaises(SystemExit):
            with Quieter():
                nucleotide_heatmap.parse_args(["max_tss_file", "region_width"])

        with self.assertRaises(SystemExit):
            with Quieter():
                nucleotide_heatmap.parse_args(["max_tss_file", "region_width", "vertical_average", "extra"])

        max_tss_file = generate_random_filename()
        with open(max_tss_file, 'w') as file:
            file.write(
                "\t".join(['chr1', '1', '2', 'name', '0', '+'])
            )

        result = nucleotide_heatmap.parse_args([max_tss_file, '50', '2000', '2'])
        self.assertEqual(result, (max_tss_file, 50, 2000, 2))

        remove_files(max_tss_file)

    def test_expand_region(self):
        max_tss_file = generate_random_filename()

        with open(max_tss_file, 'w') as file:
            file.write(
                "\t".join(["chr1", "925739", "925740", "SAMD11", "0", "+"]) + "\n" +
                "\t".join(["chr1", "959255", "959256", "NOC2L", "0", "-"]) + "\n"
            )

        region_width = 20
        expanded_region = nucleotide_heatmap.expand_region(max_tss_file, region_width)

        result = []
        with open(expanded_region) as file:
            for line in file:
                result.append(line.split())

        expected = [
            ["chr1", "925729", "925749", "SAMD11", "0", "+"],
            ["chr1", "959246", "959266", "NOC2L", "0", "-"]
        ]

        self.assertEqual(result, expected)

        remove_files(expanded_region, max_tss_file)

    def test_get_sequences(self):
        expanded_region_file = generate_random_filename()

        with open(expanded_region_file, 'w') as file:
            file.write(
                "\t".join(["chr1", "925729", "925749", "SAMD11", "0", "+"]) + "\n" +
                "\t".join(["chr1", "959246", "959266", "NOC2L", "0", "-"]) + "\n"
            )

        result = nucleotide_heatmap.get_sequences(expanded_region_file)

        expected = [
            "TGCAGAGCCCAGCAGATCCC",
            "TGGGGTGCACGCTTCGGGTT"
        ]

        self.assertEqual(result, expected)

        remove_files(expanded_region_file)

    def test_convert_sequences_to_matrix_file(self):
        sequence = ["ATGC"]

        test_nts = ["A", "T", "G", "C"]

        for i, nt in enumerate(test_nts):
            matrix_filename = nucleotide_heatmap.convert_sequences_to_matrix_file(sequence, nt, 1)

            with open(matrix_filename) as file:
                result = [int(val) for val in file.readline().split()]

            expected = [0, 0, 0, 0]

            # Set the current nt to 1 because it was found
            expected[i] = 1

            self.assertEqual(result, expected)

            remove_files(matrix_filename)

if __name__ == '__main__':
    unittest.main()
