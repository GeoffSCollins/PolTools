import unittest.mock
from pathlib import Path

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import five_prime_metaplot, three_prime_metaplot

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
            five_prime_metaplot.main([])

        with self.assertRaises(SystemExit):
            three_prime_metaplot.main([])

    def test_only_regions_file(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            five_prime_metaplot.main([self.polr2a_region_file])

        with self.assertRaises(SystemExit):
            three_prime_metaplot.main([self.polr2a_region_file])

if __name__ == '__main__':
    unittest.main()
