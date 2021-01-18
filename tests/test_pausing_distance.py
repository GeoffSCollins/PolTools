import io
import multiprocessing

import unittest.mock

import GC_bioinfo.main_programs.pausing_distance_distribution_from_maxTSS as pausing
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files

from quiter import Quieter

class TestPausingDistance(unittest.TestCase):

    def test_arguments(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                pausing.parse_input(['regions_filename'])

        regions_file = generate_random_filename()
        with open(regions_file, 'w') as file:
            file.write(
                "\t".join(['chr1', '1', '3', 'name', '0', '+'])
            )

        max_threads = multiprocessing.cpu_count()

        result = pausing.parse_input([regions_file, 'seq_file'])
        self.assertEqual(result, (regions_file, ['seq_file'], max_threads))

        # Test the threading option
        result = pausing.parse_input([regions_file, 'seq_file', '-t', '2'])
        self.assertEqual(result, (regions_file, ['seq_file'], 2))
        result = pausing.parse_input([regions_file, 'seq_file', '--threads', '2'])
        self.assertEqual(result, (regions_file, ['seq_file'], 2))

        remove_files(regions_file)

    def make_seq_file(self):
        seq_filename = generate_random_filename()

        with open(seq_filename, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "10", "20", "name", "0", "+"]
                ) + "\n" +
                "\t".join(
                    ["chr1", "10", "20", "name", "0", "+"]
                ) + "\n" +
                "\t".join(
                    ["chr1", "50", "80", "name", "0", "-"]
                ) + "\n" +
                "\t".join(
                    ["chr1", "50", "80", "name", "0", "-"]
                ) + "\n"
            )

        return seq_filename

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_strand_region_size_one(self, stdout):
        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "10", "11", "region", "0", "+"]
                ) + "\n"
            )

        seq_filename = self.make_seq_file()

        pausing.main([regions_filename, seq_filename])

        remove_files(regions_filename, seq_filename)

        result = stdout.getvalue()

        # Remove headers from result and split each line into a list
        result = [line.split() for line in result.split("\n")[1:]]

        # The result should only have the transcript length (10) and two reads at that length
        result = result[10]

        expected = ["10", "2"]

        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_strand_region_size_twenty(self, stdout):
        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "0", "20", "region", "0", "+"]
                ) + "\n"
            )

        seq_filename = self.make_seq_file()

        pausing.main([regions_filename, seq_filename])

        remove_files(regions_filename, seq_filename)

        result = stdout.getvalue()

        # Remove headers from result and split each line into a list
        result = [line.split() for line in result.split("\n")[1:]]

        # The result should only have the transcript length (10) and two reads at that length
        result = result[10]

        expected = ["10", "2"]

        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_strand_region_size_one(self, stdout):
        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "79", "80", "region", "0", "-"]
                ) + "\n"
            )

        seq_filename = self.make_seq_file()

        pausing.main([regions_filename, seq_filename])

        remove_files(regions_filename, seq_filename)

        result = stdout.getvalue()

        # Remove headers from result and split each line into a list
        result = [line.split() for line in result.split("\n")[1:]]

        # The result should only have the transcript length (10) and two reads at that length
        result = result[30]

        expected = ["30", "2"]

        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_strand_region_size_twenty(self, stdout):
        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write(
                "\t".join(
                    ["chr1", "70", "90", "region", "0", "-"]
                ) + "\n"
            )

        seq_filename = self.make_seq_file()

        pausing.main([regions_filename, seq_filename])

        remove_files(regions_filename, seq_filename)

        result = stdout.getvalue()

        # Remove headers from result and split each line into a list
        result = [line.split() for line in result.split("\n")[1:]]

        # The result should only have the transcript length (10) and two reads at that length
        result = result[30]

        expected = ["30", "2"]

        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()
