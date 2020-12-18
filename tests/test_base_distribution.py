import unittest.mock
from pathlib import Path

import io

from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files

from GC_bioinfo.main_programs import base_distribution

from quiter import Quieter


class TestBaseDistribution(unittest.TestCase):

    polr2a_region_file = str(Path(__file__).parent) + "/test_files/POLR2A_inr.bed"
    # Has the sequence "CGAGTTCGCT GCTCAGAAGC" (+1 on right side of space)

    ccnt1_region_file = str(Path(__file__).parent) + "/test_files/CCNT1_inr.bed"
    # Has the sequence "AAGTGCCTGC AGCCTTCGCC" (+1 on right side of space)

    two_regions_file = str(Path(__file__).parent) + "/test_files/two_regions.bed"
    # Has both POLR2A and CCNT1

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
            with Quieter():
                base_distribution.main([])

    def test_two_arguments(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                base_distribution.main([self.polr2a_region_file, self.ccnt1_region_file])

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_polr2a(self, stdout):
        base_distribution.main([self.polr2a_region_file])

        result = self.get_sequence(stdout.getvalue())
        self.assertEqual(result, "CGAGTTCGCTGCTCAGAAGC")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_ccnt1(self, stdout):
        base_distribution.main([self.ccnt1_region_file])

        result = self.get_sequence(stdout.getvalue())
        self.assertEqual(result, "AAGTGCCTGCAGCCTTCGCC")

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_with_n_in_sequence(self, stdout):
        region_file = generate_random_filename()

        with open(region_file, 'w') as file:
            file.write("chr1\t1\t11\tname\t0\t+")

        self.assertFalse(self.get_sequence(stdout.getvalue()))

        base_distribution.main([region_file])
        remove_files(region_file)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_two_regions(self, stdout):
        base_distribution.main([self.two_regions_file])

        result = []

        for i, line in enumerate(stdout.getvalue().split("\n")):
            if i != 0 and len(line.split()) != 0:
                position, a, t, g, c = line.split()

                curr_val = {
                    "A": float(a),
                    "C": float(c),
                    "G": float(g),
                    "T": float(t)
                }

                result.append(curr_val)

        # Just do the expected for the first 4
        expected = [
            {
                "A": 0.5,
                "C": 0.5,
                "G": 0,
                "T": 0
            },
            {
                "A": 0.5,
                "C": 0,
                "G": 0.5,
                "T": 0
            },
            {
                "A": 0.5,
                "C": 0,
                "G": 0.5,
                "T": 0
            },
            {
                "A": 0,
                "C": 0,
                "G": 0.5,
                "T": 0.5
            },

        ]

        self.assertEqual(result[0], expected[0])
        self.assertEqual(result[1], expected[1])
        self.assertEqual(result[2], expected[2])
        self.assertEqual(result[3], expected[3])


if __name__ == '__main__':
    unittest.main()
