import io
import multiprocessing

import unittest.mock

from GC_bioinfo.main_programs import tps_distance_per_gene

from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.make_random_filename import generate_random_filename

from quieter import Quieter

class TestTPSDistancePerGene(unittest.TestCase):
    def test_arguments(self):
        # Should print the usage
        with self.assertRaises(SystemExit):
            with Quieter():
                tps_distance_per_gene.parse_args([])

        with self.assertRaises(SystemExit):
            with Quieter():
                tps_distance_per_gene.parse_args(["regions filename"])

        regions_file = generate_random_filename()
        with open(regions_file, 'w') as file:
            file.write(
                "\t".join(['chr1', '1', '3', 'name', '0', '+'])
            )

        max_threads = multiprocessing.cpu_count()
        result = tps_distance_per_gene.parse_args([regions_file, 'seq_file'])
        self.assertEqual(result, (regions_file, ['seq_file'], 2, max_threads))

        result = tps_distance_per_gene.parse_args([regions_file, 'seq_file', '-t', '2'])
        self.assertEqual(result, (regions_file, ['seq_file'], 2, 2))

        result = tps_distance_per_gene.parse_args([regions_file, 'seq_file', '--threads', '2'])
        self.assertEqual(result, (regions_file, ['seq_file'], 2, 2))



    def test_get_pausing_distances_helper(self):
        region_filename = generate_random_filename()

        with open(region_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "100", "101", "positive_gene", "0", "+"]) + "\n" +
                "\t".join(["chr1", "9999", "10000", "negative_gene", "0", "-"]) + "\n"
            )

        transcripts_dict = {
            "chr1": {
                "+": {
                    100: [200, 300, 800, 200, 200, 283]
                },
                "-": {
                    9999: [9000, 9000, 9000, 8000, 7050, 6542]
                }
            }
        }

        result = tps_distance_per_gene.get_pausing_distances_helper(region_filename, transcripts_dict, 1)

        expected = {
            "positive_gene": 101,
            "negative_gene": 999
        }

        self.assertDictEqual(result, expected)

        remove_files(region_filename)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_complete_run(self, stdout):
        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "100", "101", "positive_gene", "0", "+"]) + "\n" +
                "\t".join(["chr1", "9999", "10000", "negative_gene", "0", "-"]) + "\n"
            )

        sequencing_file = generate_random_filename()

        with open(sequencing_file, 'w') as file:
            file.write(
                "\t".join(["chr1", "100", "201", "name", "0", "+"]) + "\n" +
                "\t".join(["chr1", "100", "301", "name", "0", "+"]) + "\n" +
                "\t".join(["chr1", "100", "801", "name", "0", "+"]) + "\n" +
                "\t".join(["chr1", "100", "201", "name", "0", "+"]) + "\n" +
                "\t".join(["chr1", "100", "201", "name", "0", "+"]) + "\n" +
                "\t".join(["chr1", "100", "284", "name", "0", "+"]) + "\n" +

                "\t".join(["chr1", "9000", "10000", "name", "0", "-"]) + "\n" +
                "\t".join(["chr1", "9000", "10000", "name", "0", "-"]) + "\n" +
                "\t".join(["chr1", "9000", "10000", "name", "0", "-"]) + "\n" +
                "\t".join(["chr1", "8000", "10000", "name", "0", "-"]) + "\n" +
                "\t".join(["chr1", "7050", "10000", "name", "0", "-"]) + "\n" +
                "\t".join(["chr1", "6542", "10000", "name", "0", "-"]) + "\n"
            )

        tps_distance_per_gene.main([regions_filename, sequencing_file])

        result = stdout.getvalue()

        # Eliminate the headers
        result = result.split("\n")[1:]

        # Put the result into a list
        result = [line.split() for line in result if line]

        expected = [
            ["positive_gene", "101"],
            ["negative_gene", "999"]
        ]

        self.assertEqual(result, expected)

        remove_files(regions_filename, sequencing_file)


if __name__ == '__main__':
    unittest.main()
