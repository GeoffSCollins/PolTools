import unittest.mock
import io
import multiprocessing

from GC_bioinfo.main_programs import metaplot

from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files
from quieter import Quieter

class TestFiveAndThreeMetaplots(unittest.TestCase):
    """
    This test file also tests the run_metaplot util
    """

    def test_parse_input(self):
        # No arguments throws error
        with self.assertRaises(SystemExit):
            with Quieter():
                metaplot.parse_input([])

        # Needs a region file and at least one seq file
        regions_file = generate_random_filename()
        with open(regions_file, 'w') as file:
            file.write(
                "\t".join(['chr1', '1', '3', 'name', '0', '+'])
            )

        region_length = 2
        max_threads = multiprocessing.cpu_count()

        # These will work!!
        result = metaplot.parse_input(['five', regions_file, 'seq_file'])
        self.assertEqual(result, ('five', regions_file, ['seq_file'], 2, max_threads))
        result = metaplot.parse_input(['three', regions_file, 'seq_file'])
        self.assertEqual(result, ('three', regions_file, ['seq_file'], 2, max_threads))

        # Test the threading is working
        result = metaplot.parse_input(['five', regions_file, 'seq_file', '-t', '4'])
        self.assertEqual(result, ('five', regions_file, ['seq_file'], 2, 4))

        result = metaplot.parse_input(['three', regions_file, 'seq_file', '--threads', '4'])
        self.assertEqual(result, ('three', regions_file, ['seq_file'], 2, 4))

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

        metaplot.main(['five', region_filename, seq_filename])

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

        metaplot.main(['three', region_filename, seq_filename])

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

        metaplot.main(['five', region_filename, seq_filename])

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

        metaplot.main(['three', region_filename, seq_filename])

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

        metaplot.main(['five', region_filename, seq_filename, seq_filename_two])

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

        metaplot.main(['three', region_filename, seq_filename, seq_filename_two])

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

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_complete_run(self, stdout):
        reads_filename = generate_random_filename()

        with open(reads_filename, 'w') as file:
            file.write("chr1\t1\t10\tname\t0\t+\n")
            file.write("chr1\t10\t20\tname\t0\t-\n")

        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write("chr1\t0\t20\tname\t0\t+\n")

        metaplot.main(['whole', regions_filename, reads_filename])

        output = stdout.getvalue().split("\n")[1:]

        result = []

        for line in output:
            if line:
                result.append(
                    tuple([float(val) for val in line.split()])
                )

        remove_files(reads_filename, regions_filename)

        position = list(range(-10, 0)) + list(range(1, 11))
        fw_expected = [0] + [1] * 9 + [0] * 10
        rv_expected = [0] * 10 + [-1] * 10

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

        metaplot.main(['whole', regions_filename, reads_filename, reads_filename_two])

        output = stdout.getvalue().split("\n")
        header = output[0].split("\t")

        reads_basename = reads_filename.split("/")[-1]
        reads_basename_two = reads_filename_two.split("/")[-1]

        expected_header = ["Position", reads_basename + " whole sense strand", reads_basename + " whole divergent strand",
                           reads_basename_two + " whole sense strand", reads_basename_two + " whole divergent strand"]

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
