import unittest.mock

import multiprocessing
import io

from GC_bioinfo.main_programs import divergent_pileup_metaplot
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files

from quiter import Quieter


class TestDivergentPileup(unittest.TestCase):

    def make_example_file(self):
        filename = generate_random_filename()

        with open(filename, 'w') as file:
            file.write("chr1\t1\t10\tname\t0\t+\n")
            file.write("chr1\t80\t90\tname\t0\t-\n")

        return filename

    def test_arguments(self):
        # No arguments throws error
        with self.assertRaises(SystemExit):
            with Quieter():
                divergent_pileup_metaplot.parse_args([])

        # One argument throws an error
        with self.assertRaises(SystemExit):
            with Quieter():
                divergent_pileup_metaplot.parse_args(['one'])

        regions_file = generate_random_filename()
        with open(regions_file, 'w') as file:
            file.write(
                "\t".join(['chr1', '1', '3', 'name', '0', '+'])
            )

        region_length = 2
        max_cpu_count = multiprocessing.cpu_count()

        # Two args and three args are fine
        result = divergent_pileup_metaplot.parse_args([regions_file, 'seq_file_one'])
        self.assertEqual(result, (regions_file, ['seq_file_one'], region_length, max_cpu_count))
        result = divergent_pileup_metaplot.parse_args([regions_file, ['seq_file_one', 'seq_file_two']])
        # I don't know why, but in the test arguments with a list value must be inside a list and inside a string :O
        self.assertEqual(result, (regions_file, ["['seq_file_one', 'seq_file_two']"], region_length, max_cpu_count))

        # Test number of threads is working
        result = divergent_pileup_metaplot.parse_args([regions_file, 'seq_file_one', '-t', '4'])
        self.assertEqual(result, (regions_file, ['seq_file_one'], region_length, 4))
        result = divergent_pileup_metaplot.parse_args([regions_file, 'seq_file_one', '--threads', '4'])
        self.assertEqual(result, (regions_file, ['seq_file_one'], region_length, 4))

        # Threads can't be 0 (negative values aren't recorded)
        with self.assertRaises(SystemExit):
            with Quieter():
                divergent_pileup_metaplot.parse_args([regions_file, 'seq_file_one', '-t', '0'])

        remove_files(regions_file)

    def test_split_bed_file(self):
        filename = self.make_example_file()

        fw_filename, rv_filename = divergent_pileup_metaplot.split_bed_file(filename)

        with open(fw_filename) as file:
            fw_result = file.readline().split()

        with open(rv_filename) as file:
            rv_result = file.readline().split()

        fw_expected = ["chr1", "1", "10", "name", "0", "+"]
        rv_expected = ["chr1", "80", "90", "name", "0", "-"]

        self.assertEqual(fw_expected, fw_result)
        self.assertEqual(rv_expected, rv_result)

        remove_files(filename, fw_filename, rv_filename)

    def test_make_rev_region_file(self):
        filename = self.make_example_file()

        rev_regions_filename = divergent_pileup_metaplot.make_rev_region_file(filename)

        with open(rev_regions_filename) as file:
            result = [line.split() for line in file if line]

        expected = [
            ["chr1", "1", "10", "name", "0", "-"],
            ["chr1", "80", "90", "name", "0", "+"]
        ]

        self.assertEqual(expected, result)

        remove_files(filename, rev_regions_filename)

    def test_get_pileups_helper(self):
        reads_filename = self.make_example_file()

        seq_files = divergent_pileup_metaplot.split_bed_file(reads_filename)

        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write("chr1\t0\t20\tregion\t0\t+")

        from GC_bioinfo.utils.get_region_length import determine_region_length
        region_length = determine_region_length(regions_filename)

        region_files = divergent_pileup_metaplot.split_bed_file(regions_filename)

        result = divergent_pileup_metaplot.get_pileups_helper(region_files, seq_files, region_length, opposite_strand=False)
        expected = [0] + [1]*9 + [0]*10

        remove_files(reads_filename, seq_files, regions_filename, region_files)

        self.assertEqual(result, expected)

    def test_get_pileups_helper_opposite_strand(self):
        reads_filename = self.make_example_file()

        seq_files = divergent_pileup_metaplot.split_bed_file(reads_filename)

        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write("chr1\t0\t20\tregion\t0\t+")

        from GC_bioinfo.utils.get_region_length import determine_region_length
        region_length = determine_region_length(regions_filename)

        region_files = divergent_pileup_metaplot.split_bed_file(regions_filename)

        result = divergent_pileup_metaplot.get_pileups_helper(region_files, seq_files, region_length, opposite_strand=True)
        expected = [0]*10 + [-1] * 9 + [0]

        remove_files(reads_filename, seq_files, regions_filename, region_files)

        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_complete_run(self, stdout):
        reads_filename = generate_random_filename()

        with open(reads_filename, 'w') as file:
            file.write("chr1\t1\t10\tname\t0\t+\n")
            file.write("chr1\t10\t20\tname\t0\t-\n")

        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write("chr1\t0\t20\tname\t0\t+\n")

        divergent_pileup_metaplot.main([regions_filename, reads_filename])

        output = stdout.getvalue().split("\n")[1:]

        result = []

        for line in output:
            if line:
                result.append(
                    tuple([float(val) for val in line.split()])
                )

        remove_files(reads_filename, regions_filename)

        position = list(range(-10, 0)) + list(range(1, 11))
        fw_expected = [0] + [1]*9 + [0]*10
        rv_expected = [0]*10 + [-1]*10

        expected = list(zip(position, fw_expected, rv_expected))

        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_two_sequencing_files(self, stdout):
        reads_filename = generate_random_filename()
        reads_filename_two = generate_random_filename()

        with open(reads_filename, 'w') as file:
            file.write("chr1\t1\t10\tname\t0\t+\n")
            file.write("chr1\t10\t20\tname\t0\t-\n")

        with open(reads_filename_two, 'w') as file:
            file.write("chr1\t1\t10\tname\t0\t-\n")
            file.write("chr1\t10\t20\tname\t0\t+\n")

        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write("chr1\t0\t20\tname\t0\t+\n")

        divergent_pileup_metaplot.main([regions_filename, reads_filename, reads_filename_two])

        output = stdout.getvalue().split("\n")
        header = output[0].split("\t")

        reads_basename = reads_filename.split("/")[-1]
        reads_basename_two = reads_filename_two.split("/")[-1]

        expected_header = ["Position", reads_basename + " sense strand", reads_basename + " divergent strand",
                           reads_basename_two + " sense strand", reads_basename_two + " divergent strand"]

        self.assertEqual(header, expected_header)

        result = []

        for line in output[1:]:
            if line:
                result.append(
                    tuple([float(val) for val in line.split()])
                )

        remove_files(reads_filename, regions_filename, reads_filename_two)

        position = list(range(-10, 0)) + list(range(1, 11))
        fw_expected = [0] + [1] * 9 + [0] * 10
        rv_expected = [0] * 10 + [-1] * 10

        fw_expected_two = [0] * 10 + [1] * 10
        rv_expected_two = [0] + [-1] * 9 + [0] * 10

        expected = list(zip(position, fw_expected, rv_expected, fw_expected_two, rv_expected_two))

        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
