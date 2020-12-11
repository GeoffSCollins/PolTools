import unittest.mock
import sys
import io

from main_programs import five_prime_metaplot, three_prime_metaplot

from utils.make_random_filename import generate_random_filename
from utils.remove_files import remove_files
from quiet_stderr import Quieter

sys.path.append("../GC_bioinfo")

class TestFiveAndThreeMetaplots(unittest.TestCase):
    """
    This test file also tests the run_metaplot util
    """

    def test_no_arguments(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            with Quieter():
                five_prime_metaplot.main([])

        with self.assertRaises(SystemExit):
            with Quieter():
                three_prime_metaplot.main([])

    def test_only_regions_file(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            with Quieter():
                five_prime_metaplot.main(["placeholder"])

        with self.assertRaises(SystemExit):
            with Quieter():
                three_prime_metaplot.main(["placeholder"])

    def test_three_arguemnts(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            with Quieter():
                five_prime_metaplot.main(["placeholder"])

        with self.assertRaises(SystemExit):
            with Quieter():
                three_prime_metaplot.main(["placeholder"])

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
