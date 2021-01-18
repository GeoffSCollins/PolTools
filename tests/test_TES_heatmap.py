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
        pass

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
