import unittest.mock

from PolTools.main_programs import gene_body_heatmap, gene_body_combined_heatmap, gene_body_fold_change_heatmap

from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.remove_files import remove_files
from quieter import Quieter

class TestGeneBodyHeatmap(unittest.TestCase):

    def test_make_incremented_regions(self):
        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write(
                "\t".join(
                    ["gene_name", "chromosome", "pause_left", "pause_right", "strand", "total_reads", "max_tss",
                     "max_tss_five_prime_reads", "avg_tss", "std_tss", "gene_body_left", "gene_body_right"]
                ) + "\n" +
                "\t".join(
                    # gene_name, chromosome, pause_left, pause_right, strand, total_reads, max_tss, max_tss_five_prime_reads, avg_tss
                    # std_tss, gene_body_left, gene_body_right
                    ["positive_gene", "chr1", "100", "250", "+", "20", "200", "20", "200", "0", "250", "300"]
                ) + "\n" +
                "\t".join(
                    ["negative_gene", "chr1", "4850", "5000", "-", "20", "4900", "20", "4900", "0", "4700", "4850"]
                )
            )

        downstream_distance = 100
        interval_size = 50
        upstream_distance = 50

        region_intervals_filename = gene_body_heatmap.make_incremented_regions(regions_filename,
                                                                               downstream_distance,
                                                                               interval_size,
                                                                               upstream_distance)

        result = []

        with open(region_intervals_filename) as file:
            for line in file:
                result.append(line.split())

        # Not finding the positive genes for some reason
        expected = [
            ["chr1", "150", "200", "positive_gene", "20", "+"],
            ["chr1", "200", "250", "positive_gene", "20", "+"],
            ["chr1", "250", "300", "positive_gene", "20", "+"],
            ["chr1", "300", "350", "positive_gene", "20", "+"],
            ["chr1", "350", "400", "positive_gene", "20", "+"],

            ["chr1", "4900", "4950", "negative_gene", "20", "-"],
            ["chr1", "4850", "4900", "negative_gene", "20", "-"],
            ["chr1", "4800", "4850", "negative_gene", "20", "-"],
            ["chr1", "4750", "4800", "negative_gene", "20", "-"],
            ["chr1", "4700", "4750", "negative_gene", "20", "-"],
            ["chr1", "4650", "4700", "negative_gene", "20", "-"],
            ["chr1", "4600", "4650", "negative_gene", "20", "-"],
        ]

        self.assertEqual(result, expected)

        remove_files(regions_filename, region_intervals_filename)

    def test_read_coverage_file(self):

        coverage_file = generate_random_filename()

        with open(coverage_file, 'w') as file:
            file.write(
                "\t".join(["chr1", "left", "right", "positive_gene", "0", "+", "1", "0", "0", "0"]) + "\n" +
                "\t".join(["chr1", "left", "right", "positive_gene", "0", "+", "2", "0", "0", "0"]) + "\n" +
                "\t".join(["chr1", "left", "right", "positive_gene", "0", "+", "4", "0", "0", "0"]) + "\n" +
                "\t".join(["chr1", "left", "right", "positive_gene", "0", "+", "10", "0", "0", "0"]) + "\n" +
                "\t".join(["chr1", "left", "right", "positive_gene", "0", "+", "10", "0", "0", "0"]) + "\n" +

                "\t".join(["chr1", "left", "right", "negative_gene", "0", "-", "1", "0", "0", "0"]) + "\n" +
                "\t".join(["chr1", "left", "right", "negative_gene", "0", "-", "1", "0", "0", "0"]) + "\n" +
                "\t".join(["chr1", "left", "right", "negative_gene", "0", "-", "1", "0", "0", "0"]) + "\n" +
                "\t".join(["chr1", "left", "right", "negative_gene", "0", "-", "1", "0", "0", "0"]) + "\n" +
                "\t".join(["chr1", "left", "right", "negative_gene", "0", "-", "1", "0", "0", "0"]) + "\n" +

                "\t".join(["chr1", "left", "right", "third_gene", "0", "-", "5", "0", "0", "0"]) + "\n" +
                "\t".join(["chr1", "left", "right", "third_gene", "0", "-", "5", "0", "0", "0"]) + "\n" +
                "\t".join(["chr1", "left", "right", "third_gene", "0", "-", "5", "0", "0", "0"]) + "\n"
            )

        # Make the matrix width 7, meaning that two more columns of 0s need to be added
        matrix_filename, matrix_height = gene_body_heatmap.read_coverage_file(coverage_file, 7)

        result = []
        with open(matrix_filename) as file:
            for line in file:
                result.append(line.split())

        expected = [
            ["5", "5", "5", "0", "0", "0", "0"],
            ["1", "2", "4", "10", "10", "0", "0"],
            ["1", "1", "1", "1", "1", "0", "0"]
        ]

        self.assertEqual(result, expected)
        self.assertEqual(matrix_height, 3)

        remove_files(coverage_file, matrix_filename)


if __name__ == '__main__':
    unittest.main()
