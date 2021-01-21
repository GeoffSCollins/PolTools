import unittest.mock
import io

from GC_bioinfo.main_programs import TES_heatmap, TES_combined_heatmap, TES_fold_change_heatmap

from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files
from quiter import Quieter


# TODO
class TestTESHeatmap(unittest.TestCase):
    def test_invalid_number_of_arguments(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            with Quieter():
                TES_heatmap.main([])

    def test_make_incremented_regions(self):

        self.maxDiff = None

        truQuant_output_file = generate_random_filename('-truQuant_output.txt')

        tQ_text = """  Gene    Chromosome      Pause Region Left       Pause Region Right      Strand  Total 5' Reads  MaxTSS  MaxTSS 5' Reads Weighted Pause Region Center    STDEV of TSSs   Gene Body Left  Gene Body Right Gene Body Distance      seq_file.bed Pause Region   seq_file.bed Gene Body
  negative_gene   chr1    5000  5150  -       194     5100  46      5100  13.306459171023036      4000  5000  600   194     18
  positive_gene  chr1    1000  1150  +       234     1100  27      1100  25.417791063821863      1150  2000  850    234     17"""

        tQ_text = [line.split() for line in tQ_text.split("\n") if line]

        with open(truQuant_output_file, 'w') as file:
            for line in tQ_text:
                file.write("\t".join(line) + "\n")

        downstream_distance = 0
        upstream_distance = 200
        bp_width = 1000
        interval_size = 200

        incremented_regions_file = TES_heatmap.make_incremented_regions(truQuant_output_file,
                                                                        downstream_distance,
                                                                        upstream_distance,
                                                                        bp_width,
                                                                        interval_size)


        result = []
        with open(incremented_regions_file) as file:
            for line in file:
                result.append(line.split())

        expected = [
            ["chr1", "5100", "5300", "negative_gene", "46", "-"],
            ["chr1", "4900", "5100", "negative_gene", "46", "-"],
            ["chr1", "4700", "4900", "negative_gene", "46", "-"],
            ["chr1", "4500", "4700", "negative_gene", "46", "-"],
            ["chr1", "4300", "4500", "negative_gene", "46", "-"],

            ["chr1", "900", "1100", "positive_gene", "27", "+"],
            ["chr1", "1100", "1300", "positive_gene", "27", "+"],
            ["chr1", "1300", "1500", "positive_gene", "27", "+"],
            ["chr1", "1500", "1700", "positive_gene", "27", "+"],
            ["chr1", "1700", "1900", "positive_gene", "27", "+"],
        ]

        remove_files(truQuant_output_file)
        self.assertEqual(result, expected)

    def test_quantify_intervals(self):
        pass

    def test_read_coverage_file(self):
        pass

    def test_build_matrix(self):
        pass


class TestCombineTESHeatmap(unittest.TestCase):
    def test_invalid_number_of_arguments(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            with Quieter():
                TES_combined_heatmap.main([])


class TestTESFoldChangeHeatmap(unittest.TestCase):
    def test_invalid_number_of_arguments(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            with Quieter():
                TES_fold_change_heatmap.main([])


if __name__ == '__main__':
    unittest.main()
