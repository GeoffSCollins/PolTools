import unittest.mock
from pathlib import Path

import sys

sys.path.append("../GC_bioinfo")

from utils.remove_files import remove_files

from other_programs import rolling_average_bigwigs

from quiet_stderr import Quieter

"""
python3 rolling_average_bigwigs.py 101 /media/genomes/hg38/hg38.chrom.sizes /media/genomes/hg38/hg38.fa
python3 rolling_average_bigwigs.py 31 /media/genomes/hg38/hg38.chrom.sizes /media/genomes/hg38/hg38.fa


python3 GC_bioinfo/GC_bioinfo/rolling_average_bigwigs.py 1 /media/genomes/KF297339.1/KF297339.1.chrom.size /media/genomes/KF297339.1/KF297339.1.fa
python3 GC_bioinfo/GC_bioinfo/rolling_average_bigwigs.py 5 /media/genomes/KF297339.1/KF297339.1.chrom.size /media/genomes/KF297339.1/KF297339.1.fa
python3 GC_bioinfo/GC_bioinfo/rolling_average_bigwigs.py 11 /media/genomes/KF297339.1/KF297339.1.chrom.size /media/genomes/KF297339.1/KF297339.1.fa
python3 GC_bioinfo/GC_bioinfo/rolling_average_bigwigs.py 31 /media/genomes/KF297339.1/KF297339.1.chrom.size /media/genomes/KF297339.1/KF297339.1.fa
python3 GC_bioinfo/GC_bioinfo/rolling_average_bigwigs.py 101 /media/genomes/KF297339.1/KF297339.1.chrom.size /media/genomes/KF297339.1/KF297339.1.fa

"""


class TestRollingAverageBigWigs(unittest.TestCase):

    chrom_sizes_file = str(Path(__file__).parent) + "/test_files/hg38_chrom_sizes.txt"
    sample_fasta_file = str(Path(__file__).parent) + "/test_files/hg38_sample.fa"
    small_fasta_file = str(Path(__file__).parent) + "/test_files/small_hg38_sample.fa"


    def test_no_arguments(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                rolling_average_bigwigs.main([])

    def test_one_argument(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                rolling_average_bigwigs.main([1])

    def test_two_arguments(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                rolling_average_bigwigs.main([1, self.chrom_sizes_file])

    def test_make_bedgraphs(self):
        sys.stderr.write("HERE: " + self.small_fasta_file)
        with Quieter():
            bedgraphs = rolling_average_bigwigs.make_bedgraphs(1, self.small_fasta_file)

        a_bedgraph, t_bedgraph, g_bedgraph, c_bedgraph = bedgraphs

        with open(a_bedgraph) as file:
            a_result = [line.split() for line in file.readlines()]

        a_expected = [
            ["chr1", "0", "1", "1"],
            ["chr1", "1", "2", "1"],
            ["chr1", "2", "3", "1"],
            ["chr1", "3", "4", "1"],
            ["chr1", "9", "10", "1"]
        ]
        self.assertEqual(a_expected, a_result)

        with open(t_bedgraph) as file:
            t_result = [line.split() for line in file.readlines()]

        t_expected = [
            ["chr1", "4", "5", "1"],
            ["chr1", "7", "8", "1"]
        ]
        self.assertEqual(t_expected, t_result)

        with open(g_bedgraph) as file:
            g_result = [line.split() for line in file.readlines()]

        g_expected = [
            ["chr1", "5", "6", "1"]
        ]
        self.assertEqual(g_expected, g_result)

        with open(c_bedgraph) as file:
            c_result = [line.split() for line in file.readlines()]

        c_expected = [
            ["chr1", "6", "7", "1"],
            ["chr1", "8", "9", "1"]

        ]
        self.assertEqual(c_expected, c_result)

        remove_files(bedgraphs)



if __name__ == '__main__':
    unittest.main()
