import unittest.mock

from GC_bioinfo.main_programs import read_through_transcription
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files

from quiter import Quieter

class TestReadThroughTranscription(unittest.TestCase):
    def test_make_incremented_regions(self):

        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "10", "100", "positive", "0", "+"]) + "\n" +
                "\t".join(["chr1", "2000", "9815", "negative", "0", "-"]) + "\n"
            )


        upstream_distance = 10
        downstream_distance = 40
        interval_size = 10

        region_intervals_filename = read_through_transcription.make_incremented_regions(regions_filename,
                                                                                        upstream_distance,
                                                                                        downstream_distance,
                                                                                        interval_size)
        with open(region_intervals_filename) as file:
            result = []
            for line in file:
                result.append(line.split())

        expected = [
            ["chr1", "90", "100", "positive", "0", "+"],
            ["chr1", "100", "110", "positive", "0", "+"],
            ["chr1", "110", "120", "positive", "0", "+"],
            ["chr1", "120", "130", "positive", "0", "+"],
            ["chr1", "130", "140", "positive", "0", "+"],

            ["chr1", "2000", "2010", "negative", "0", "-"],
            ["chr1", "1990", "2000", "negative", "0", "-"],
            ["chr1", "1980", "1990", "negative", "0", "-"],
            ["chr1", "1970", "1980", "negative", "0", "-"],
            ["chr1", "1960", "1970", "negative", "0", "-"],
        ]

        self.assertEqual(result, expected)

        remove_files(regions_filename, region_intervals_filename)

    def test_no_arguments(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                read_through_transcription.main([])

    def test_incorrect_number_of_arguemnts(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                read_through_transcription.main(['Regions Filename', 'TSR Filename', 'Upstream Distance', 'Downstream Distance', 'Sequencing Files'])

        with self.assertRaises(SystemExit):
            with Quieter():
                read_through_transcription.main(['Regions Filename', 'TSR Filename', 'Upstream Distance', 'Downstream Distance'])

        with self.assertRaises(SystemExit):
            with Quieter():
                read_through_transcription.main(['Regions Filename', 'TSR Filename', 'Upstream Distance'])

        with self.assertRaises(SystemExit):
            with Quieter():
                read_through_transcription.main(['Regions Filename', 'TSR Filename'])

        with self.assertRaises(SystemExit):
            with Quieter():
                read_through_transcription.main(['Regions Filename'])


if __name__ == '__main__':
    unittest.main()
