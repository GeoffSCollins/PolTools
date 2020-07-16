import unittest.mock
from pathlib import Path

import sys
import os
import io

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import base_distribution


class TestBaseDistribution(unittest.TestCase):

    polr2a_region_file = str(Path(__file__).parent) + "/test_files/POLR2A_inr.bed"
    # Has the sequence "CGAGTTCGCT GCTCAGAAGC" (+1 on right side of space)

    ccnt1_region_file = str(Path(__file__).parent) + "/test_files/CCNT1_inr.bed"
    # Has the sequence "AAGTGCCTGC AGCCTTCGCC" (+1 on right side of space)

    def get_sequence(self, std_output):
        sequence = ""

        for i, line in enumerate(std_output.split("\n")):
            if i != 0 and len(line.split()) != 0:
                position, a, t, g, c = line.split()

                if a == '1.0':
                    sequence += 'A'
                    continue
                if t == '1.0':
                    sequence += 'T'
                    continue
                if g == '1.0':
                    sequence += 'G'
                    continue
                if c == '1.0':
                    sequence += 'C'
                    continue

                # Should never get to this location as that would mean no base was found
                return False

        return sequence


    def test_no_arguments(self):
        with self.assertRaises(SystemExit):
            base_distribution.main([])

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_polr2a(self, mock_stdout):
        base_distribution.main([self.polr2a_region_file])

        result = self.get_sequence(mock_stdout.getvalue())
        self.assertEqual(result, "CGAGTTCGCTGCTCAGAAGC")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_ccnt1(self, mock_stdout):
        base_distribution.main([self.ccnt1_region_file])

        result = self.get_sequence(mock_stdout.getvalue())
        self.assertEqual(result, "AAGTGCCTGCAGCCTTCGCC")

    def test_two_files(self):
        with self.assertRaises(SystemExit):
            base_distribution.main(['', ''])


if __name__ == '__main__':
    unittest.main()
