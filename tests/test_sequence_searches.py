import unittest.mock
from pathlib import Path

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import sequence_searches
from sequence_searches import InvalidSearchException

import io


class TestSequenceSearches(unittest.TestCase):
    polr2a_region_file = str(Path(__file__).parent) + "/test_files/POLR2A_inr.bed"
    # Has the sequence "CGAGTTCGCT GCTCAGAAGC" (+1 on right side of space)

    ccnt1_region_file = str(Path(__file__).parent) + "/test_files/CCNT1_inr.bed"
    # Has the sequence "AGTGCCTGCA GCCTTCGCCG" (+1 on right side of space)

    def _parse_stdoutput(self, output):
        # Returns the True or False value of the test
        return output.split("\n")[1].split()[-1]

    def test_no_arguments(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            sequence_searches.main([])

    def test_only_regions_file(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            sequence_searches.main([self.polr2a_region_file])

    def test_incorrect_search(self):
        # Should raise an InvalidSearch error because of the underscore
        with self.assertRaises(InvalidSearchException):
            sequence_searches.main([self.polr2a_region_file, 'TATA,-5_7'])

    def test_search_not_in_region(self):
        # Should raise an InvalidSearch error because polr2a file is -10 to +10
        with self.assertRaises(InvalidSearchException):
            sequence_searches.main([self.polr2a_region_file, 'TATA,-20:-10'])

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_strand_one_upstream_sequence_not_found(self, mock_stdout):
        sequence_searches.main([self.polr2a_region_file, 'TTT,-10:-5'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "False")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_strand_one_upstream_sequence_found(self, mock_stdout):
        sequence_searches.main([self.polr2a_region_file, 'CGAGT,-10:-5'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "True")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_strand_one_upstream_sequence_not_found(self, mock_stdout):
        sequence_searches.main([self.ccnt1_region_file, 'TTT,-10:-5'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "False")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_strand_one_upstream_sequence_found(self, mock_stdout):
        sequence_searches.main([self.ccnt1_region_file, 'AGTGC,-10:-5'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "True")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_strand_one_downstream_sequence_not_found(self, mock_stdout):
        sequence_searches.main([self.polr2a_region_file, 'TTT,5:10'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "False")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_strand_one_downstream_sequence_found(self, mock_stdout):
        sequence_searches.main([self.polr2a_region_file, 'GAAGC,5:10'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "True")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_strand_one_downstream_sequence_not_found(self, mock_stdout):
        sequence_searches.main([self.ccnt1_region_file, 'TTT,5:10'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "False")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_strand_one_downstream_sequence_found(self, mock_stdout):
        sequence_searches.main([self.ccnt1_region_file, 'CGCCG,5:10'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "True")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_strand_inr_sequence_not_found(self, mock_stdout):
        sequence_searches.main([self.polr2a_region_file, 'TTT,-3:3'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "False")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_strand_inr_sequence_found(self, mock_stdout):
        sequence_searches.main([self.polr2a_region_file, 'GCTGCT,-3:3'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "True")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_strand_inr_sequence_found(self, mock_stdout):
        sequence_searches.main([self.ccnt1_region_file, 'TTT,-3:3'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "False")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_strand_inr_sequence_not_found(self, mock_stdout):
        sequence_searches.main([self.ccnt1_region_file, 'GCAGCC,-3:3'])
        output = self._parse_stdoutput(mock_stdout.getvalue())
        self.assertEqual(output, "True")


if __name__ == '__main__':
    unittest.main()
