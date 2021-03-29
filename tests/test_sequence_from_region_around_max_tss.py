import unittest.mock

import io

from PolTools.main_programs import sequence_from_region_around_max_tss
from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.remove_files import remove_files
from PolTools.utils.constants import hg38_chrom_sizes_random_file

from quieter import Quieter

class TestSequenceFromRegionAroundMaxTSS(unittest.TestCase):

    def test_arguments(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                sequence_from_region_around_max_tss.parse_args([])

        with self.assertRaises(SystemExit):
            with Quieter():
                sequence_from_region_around_max_tss.parse_args(["max_tss_file"])

        with self.assertRaises(SystemExit):
            with Quieter():
                sequence_from_region_around_max_tss.parse_args(["max_tss_file", 'left'])

        max_tss_file = generate_random_filename()

        with open(max_tss_file, 'w') as file:
            file.write(
                "\t".join(["chr16", "53607", "53608", "POLR3K", "0", "-"]) + "\n" +
                "\t".join(["chr16", "53872", "53873", "SNRNP25", "0", "+"]) + "\n"
            )

        result = sequence_from_region_around_max_tss.parse_args([max_tss_file, '-5', '10'])
        search = [
            [
                "-",
                5
            ],
            [
                '+',
                10
            ]
        ]
        self.assertEqual(result, (False, max_tss_file, search))

        remove_files(max_tss_file)

    def test_get_regions_file(self):

        max_tss_file = generate_random_filename()

        with open(max_tss_file, 'w') as file:
            file.write(
                "\t".join(["chr16", "53607", "53608", "POLR3K", "0", "-"]) + "\n" +
                "\t".join(["chr16", "53872", "53873", "SNRNP25", "0", "+"]) + "\n"
            )

        search = [
            ["-", 5],
            ["+", 5]
        ]

        chrom_sizes = {}
        with open(hg38_chrom_sizes_random_file) as file:
            for line in file:
                chrom, size = line.split()
                chrom_sizes[chrom] = int(size)

        region_file, gene_names = sequence_from_region_around_max_tss.get_regions_file(max_tss_file, search, chrom_sizes)

        result = []

        with open(region_file) as file:
            for line in file:
                result.append(line.split())

        expected = [
            ["chr16", "53603", "53613", "POLR3K", "0", "-"],
            ["chr16", "53867", "53877", "SNRNP25", "0", "+"]
        ]

        self.assertEqual(result, expected)

        remove_files(max_tss_file, region_file)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_complete_run(self, stdout):
        max_tss_file = generate_random_filename()

        with open(max_tss_file, 'w') as file:
            file.write(
                "\t".join(["chr16", "53607", "53608", "POLR3K", "0", "-"]) + "\n" +
                "\t".join(["chr16", "53872", "53873", "SNRNP25", "0", "+"]) + "\n"
            )

        sequence_from_region_around_max_tss.main([max_tss_file, "-5", "+5"])

        result = [line for line in stdout.getvalue().split("\n") if line]

        expected = [
            ">POLR3K",
            "GAGTTGGAGC",
            ">SNRNP25",
            "CTGGCAGTGC"
        ]

        self.assertEqual(result, expected)
        remove_files(max_tss_file)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_upstream(self, stdout):
        # Try -36 to -19
        max_tss_file = generate_random_filename()

        with open(max_tss_file, 'w') as file:
            file.write(
                "\t".join(["chr16", "53607", "53608", "POLR3K", "0", "-"]) + "\n" +
                "\t".join(["chr16", "53872", "53873", "SNRNP25", "0", "+"]) + "\n"
            )

        sequence_from_region_around_max_tss.main([max_tss_file, "-36", "-19"])

        result = [line for line in stdout.getvalue().split("\n") if line]

        expected = [
            ">POLR3K",
            "ATCGGCGCTGAGCGGCAG",
            ">SNRNP25",
            "GCGCGTGCGCGCTTGGCC"
        ]

        self.assertEqual(result, expected)
        remove_files()

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_downstream(self, stdout):
        # Try +5 to +35
        max_tss_file = generate_random_filename()

        with open(max_tss_file, 'w') as file:
            file.write(
                "\t".join(["chr16", "53607", "53608", "POLR3K", "0", "-"]) + "\n" +
                "\t".join(["chr16", "53872", "53873", "SNRNP25", "0", "+"]) + "\n"
            )

        # sequence_from_region_around_max_tss.main([max_tss_file, "+5:+35"])
        sequence_from_region_around_max_tss.main([max_tss_file, "+5", "+35"])

        result = [line for line in stdout.getvalue().split("\n") if line]

        expected = [
            ">POLR3K",
            "CCTGCGGAGTTCGAGACCATGCTGCTGTTCT",
            ">SNRNP25",
            "CGGGCAGAGCCCGGCTGAGAGGGGCGGCCCT"
        ]

        self.assertEqual(result, expected)

        remove_files(max_tss_file)



if __name__ == '__main__':
    unittest.main()
