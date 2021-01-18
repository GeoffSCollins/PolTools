import unittest.mock
import io
import multiprocessing

from GC_bioinfo.main_programs import five_prime_metaplot, three_prime_metaplot

from GC_bioinfo.utils.run_metaplot import parse_input
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files
from quiter import Quieter

class TestFiveAndThreeMetaplots(unittest.TestCase):
    """
    This test file also tests the run_metaplot util
    """

    def test_parse_input(self):
        # No arguments throws error
        with self.assertRaises(SystemExit):
            with Quieter():
                parse_input([], 'five')

        with self.assertRaises(SystemExit):
            with Quieter():
                parse_input([], 'three')

        # Needs a region file and at least one seq file
        regions_file = generate_random_filename()
        with open(regions_file, 'w') as file:
            file.write(
                "\t".join(['chr1', '1', '3', 'name', '0', '+'])
            )

        region_length = 2
        max_threads = multiprocessing.cpu_count()

        with self.assertRaises(SystemExit):
            with Quieter():
                parse_input([], 'five')

        with self.assertRaises(SystemExit):
            with Quieter():
                parse_input([], 'three')

        # These will work!!
        result = parse_input([regions_file, 'seq_file'], 'five')
        self.assertEqual(result, (regions_file, ['seq_file'], 2, max_threads))
        result = parse_input([regions_file, 'seq_file'], 'three')
        self.assertEqual(result, (regions_file, ['seq_file'], 2, max_threads))

        # Test the threading is working
        result = parse_input([regions_file, 'seq_file', '-t', '4'], 'three')
        self.assertEqual(result, (regions_file, ['seq_file'], 2, 4))

        result = parse_input([regions_file, 'seq_file', '--threads', '4'], 'three')
        self.assertEqual(result, (regions_file, ['seq_file'], 2, 4))

        remove_files(regions_file)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_read_five(self, stdout):
        region_filename = generate_random_filename()
        seq_filename = generate_random_filename()

        with open(region_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "0", "10", "name", "0", "+"])
            )

        with open(seq_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "2", "9", "name", "0", "+"])
            )

        five_prime_metaplot.main([region_filename, seq_filename])

        # Get the result from stdout by splitting into a list and making the output floats where possible
        result = [line for line in stdout.getvalue().split("\n") if line]
        result[0] = result[0].split("\t")

        for i, line in enumerate(result[1:]):
            result[i+1] = [float(val) for val in line.split()]

        seq_file_basename = seq_filename.split("/")[-1]

        expected = [
            ["Position", seq_file_basename + " 5' sense strand", seq_file_basename + " 5' divergent strand"],
            [-5, 0, 0],
            [-4, 0, 0],
            [-3, 1, 0],
            [-2, 0, 0],
            [-1, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [3, 0, 0],
            [4, 0, 0],
            [5, 0, 0]
        ]

        remove_files(region_filename, seq_filename)
        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_positive_read_three(self, stdout):
        region_filename = generate_random_filename()
        seq_filename = generate_random_filename()

        with open(region_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "0", "10", "name", "0", "+"])
            )

        with open(seq_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "2", "9", "name", "0", "+"])
            )

        three_prime_metaplot.main([region_filename, seq_filename])

        # Get the result from stdout by splitting into a list and making the output floats where possible
        result = [line for line in stdout.getvalue().split("\n") if line]
        result[0] = result[0].split("\t")

        for i, line in enumerate(result[1:]):
            result[i + 1] = [float(val) for val in line.split()]

        seq_file_basename = seq_filename.split("/")[-1]

        expected = [
            ["Position", seq_file_basename + " 3' sense strand", seq_file_basename + " 3' divergent strand"],
            [-5, 0, 0],
            [-4, 0, 0],
            [-3, 0, 0],
            [-2, 0, 0],
            [-1, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [3, 0, 0],
            [4, 1, 0],
            [5, 0, 0]
        ]

        remove_files(region_filename, seq_filename)
        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_read_five(self, stdout):
        region_filename = generate_random_filename()
        seq_filename = generate_random_filename()

        with open(region_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "0", "10", "name", "0", "+"])
            )

        with open(seq_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "2", "9", "name", "0", "-"])
            )

        five_prime_metaplot.main([region_filename, seq_filename])

        # Get the result from stdout by splitting into a list and making the output floats where possible
        result = [line for line in stdout.getvalue().split("\n") if line]
        result[0] = result[0].split("\t")

        for i, line in enumerate(result[1:]):
            result[i + 1] = [float(val) for val in line.split()]

        seq_file_basename = seq_filename.split("/")[-1]

        expected = [
            ["Position", seq_file_basename + " 5' sense strand", seq_file_basename + " 5' divergent strand"],
            [-5, 0, 0],
            [-4, 0, 0],
            [-3, 0, 0],
            [-2, 0, 0],
            [-1, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [3, 0, 0],
            [4, 0, -1],
            [5, 0, 0]
        ]

        remove_files(region_filename, seq_filename)
        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_negative_read_three(self, stdout):
        region_filename = generate_random_filename()
        seq_filename = generate_random_filename()

        with open(region_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "0", "10", "name", "0", "+"])
            )

        with open(seq_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "2", "9", "name", "0", "-"])
            )

        three_prime_metaplot.main([region_filename, seq_filename])

        # Get the result from stdout by splitting into a list and making the output floats where possible
        result = [line for line in stdout.getvalue().split("\n") if line]
        result[0] = result[0].split("\t")

        for i, line in enumerate(result[1:]):
            result[i + 1] = [float(val) for val in line.split()]

        seq_file_basename = seq_filename.split("/")[-1]

        expected = [
            ["Position", seq_file_basename + " 3' sense strand", seq_file_basename + " 3' divergent strand"],
            [-5, 0, 0],
            [-4, 0, 0],
            [-3, 0, -1],
            [-2, 0, 0],
            [-1, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [3, 0, 0],
            [4, 0, 0],
            [5, 0, 0]
        ]

        remove_files(region_filename, seq_filename)
        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_multiple_seq_files_five(self, stdout):
        region_filename = generate_random_filename()
        seq_filename = generate_random_filename()
        seq_filename_two = generate_random_filename()

        with open(region_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "0", "10", "name", "0", "+"])
            )

        with open(seq_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "2", "9", "name", "0", "-"])
            )

        with open(seq_filename_two, 'w') as file:
            file.write(
                "\t".join(["chr1", "1", "8", "name", "0", "-"])
            )

        five_prime_metaplot.main([region_filename, seq_filename, seq_filename_two])

        # Get the result from stdout by splitting into a list and making the output floats where possible
        result = [line for line in stdout.getvalue().split("\n") if line]

        result[0] = result[0].split("\t")

        for i, line in enumerate(result[1:]):
            result[i + 1] = [float(val) for val in line.split()]

        seq_file_basename = seq_filename.split("/")[-1]
        seq_file_two_basename = seq_filename_two.split("/")[-1]

        expected = [
            ["Position", seq_file_basename + " 5' sense strand", seq_file_basename + " 5' divergent strand",
             seq_file_two_basename + " 5' sense strand", seq_file_two_basename + " 5' divergent strand"],
            [-5, 0, 0, 0, 0],
            [-4, 0, 0, 0, 0],
            [-3, 0, 0, 0, 0],
            [-2, 0, 0, 0, 0],
            [-1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [2, 0, 0, 0, 0],
            [3, 0, 0, 0, -1],
            [4, 0, -1, 0, 0],
            [5, 0, 0, 0, 0]
        ]

        remove_files(region_filename, seq_filename, seq_filename_two)
        self.assertEqual(result[0][3], expected[0][3])
        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_multiple_seq_files_three(self, stdout):
        region_filename = generate_random_filename()
        seq_filename = generate_random_filename()
        seq_filename_two = generate_random_filename()

        with open(region_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "0", "10", "name", "0", "+"])
            )

        with open(seq_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "2", "9", "name", "0", "-"])
            )

        with open(seq_filename_two, 'w') as file:
            file.write(
                "\t".join(["chr1", "1", "8", "name", "0", "-"])
            )

        three_prime_metaplot.main([region_filename, seq_filename, seq_filename_two])

        # Get the result from stdout by splitting into a list and making the output floats where possible
        result = [line for line in stdout.getvalue().split("\n") if line]
        result[0] = result[0].split("\t")

        for i, line in enumerate(result[1:]):
            result[i + 1] = [float(val) for val in line.split()]

        seq_file_basename = seq_filename.split("/")[-1]
        seq_file_two_basename = seq_filename_two.split("/")[-1]

        expected = [
            ["Position", seq_file_basename + " 3' sense strand", seq_file_basename + " 3' divergent strand",
             seq_file_two_basename + " 3' sense strand", seq_file_two_basename + " 3' divergent strand"],
            [-5, 0, 0, 0, 0],
            [-4, 0, 0, 0, -1],
            [-3, 0, -1, 0, 0],
            [-2, 0, 0, 0, 0],
            [-1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [2, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [4, 0, 0, 0, 0],
            [5, 0, 0, 0, 0]
        ]

        remove_files(region_filename, seq_filename, seq_filename_two)
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
