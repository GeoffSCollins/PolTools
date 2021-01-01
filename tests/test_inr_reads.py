import unittest.mock
import io

from quiter import Quieter

import GC_bioinfo.main_programs.inr_reads as inr_reads
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files

class TestInrReads(unittest.TestCase):

    def test_no_arguments(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            with Quieter():
                inr_reads.main([])

    def test_one_argument(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            with Quieter():
                inr_reads.main(['placeholder'])

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_strand_region_size_one(self, stdout):
        # First make a bed file
        bed_file = generate_random_filename()

        with open(bed_file, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "50", "100", "name", "0", "+"]
                ) + "\n"
            )

        # Make a region file on the 5' end of the read in the bed file
        region_file = generate_random_filename()

        with open(region_file, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "50", "51", "region", "0", "+"]
                ) + '\n' +
                "\t".join(
                    ["chr1", "40", "41", "second_region", "0", "+"]
                ) + "\n"
            )

        # Now run inr_reads
        inr_reads.main([region_file, bed_file])

        result = stdout.getvalue()
        # Remove the headers from the result
        result = result.split("\n")[1:]

        # Put result into a list
        result = [line.split("\t") for line in result if line]

        expected = [
            ["region", "1"],
            ["second_region", "0"]
        ]

        self.assertEqual(result, expected)

        remove_files(bed_file, region_file)


    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_strand_region_size_one(self, stdout):
        # First make a bed file
        bed_file = generate_random_filename()

        with open(bed_file, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "50", "100", "name", "0", "-"]
                ) + "\n"
            )

        # Make a region file on the 5' end of the read in the bed file
        region_file = generate_random_filename()

        with open(region_file, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "99", "100", "region", "0", "-"]
                ) + '\n' +
                "\t".join(
                    ["chr1", "90", "91", "second_region", "0", "-"]
                ) + "\n"
            )

        # Now run inr_reads
        inr_reads.main([region_file, bed_file])

        result = stdout.getvalue()
        # Remove the headers from the result
        result = result.split("\n")[1:]

        # Put result into a list
        result = [line.split("\t") for line in result if line]

        expected = [
            ["region", "1"],
            ["second_region", "0"]
        ]

        self.assertEqual(result, expected)

        remove_files(bed_file, region_file)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_strand_region_size_twenty(self, stdout):
        # First make a bed file
        bed_file = generate_random_filename()

        with open(bed_file, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "50", "100", "name", "0", "+"]
                ) + "\n"
            )

        # Make a region file on the 5' end of the read in the bed file
        region_file = generate_random_filename()

        with open(region_file, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "40", "60", "region", "0", "+"]
                ) + '\n' +
                "\t".join(
                    ["chr1", "30", "50", "second_region", "0", "+"]
                ) + "\n"
            )

        # Now run inr_reads
        inr_reads.main([region_file, bed_file])

        result = stdout.getvalue()
        # Remove the headers from the result
        result = result.split("\n")[1:]

        # Put result into a list
        result = [line.split("\t") for line in result if line]

        expected = [
            ["region", "1"],
            ["second_region", "0"]
        ]

        self.assertEqual(result, expected)

        remove_files(bed_file, region_file)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_strand_region_size_twenty(self, stdout):
        # First make a bed file
        bed_file = generate_random_filename()

        with open(bed_file, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "50", "100", "name", "0", "-"]
                ) + "\n"
            )

        # Make a region file on the 5' end of the read in the bed file
        region_file = generate_random_filename()

        with open(region_file, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "90", "110", "region", "0", "-"]
                ) + '\n' +
                "\t".join(
                    ["chr1", "80", "100", "second_region", "0", "-"]
                ) + "\n"
            )

        # Now run inr_reads
        inr_reads.main([region_file, bed_file])

        result = stdout.getvalue()
        # Remove the headers from the result
        result = result.split("\n")[1:]

        # Put result into a list
        result = [line.split("\t") for line in result if line]

        expected = [
            ["region", "1"],
            ["second_region", "0"]
        ]

        self.assertEqual(result, expected)

        remove_files(bed_file, region_file)


if __name__ == '__main__':
    unittest.main()
