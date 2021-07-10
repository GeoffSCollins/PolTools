import unittest.mock
import io
from pathlib import Path

from PolTools.main_programs import make_regions_file_centered_on_max_tss

from quieter import Quieter

class TestMakeRegionsFileCenteredOnMaxTSS(unittest.TestCase):
    truQuant_file = str(Path(__file__).parent) + "/test_files/sample-truQuant_output.txt"

    def test_no_input(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                make_regions_file_centered_on_max_tss.parse_args([])

        with self.assertRaises(SystemExit):
            with Quieter():
                make_regions_file_centered_on_max_tss.parse_args(['truQuant_file'])

        result = make_regions_file_centered_on_max_tss.parse_args(['truQuant_file', '10'])
        self.assertEqual(result, ('truQuant_file', 10, False))

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_get_max_tss(self, stdout):
        make_regions_file_centered_on_max_tss.main([self.truQuant_file, "1"])

        result_lines = stdout.getvalue().split("\n")
        result = [line.split() for line in result_lines if line]

        expected = [
            ["chr1", "925738", "925739", "SAMD11", "86", "+"],
            ["chr1", "959254", "959255", "NOC2L", "241", "-"]
        ]

        self.assertEqual(result, expected)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_twenty_bp_region(self, stdout):
        make_regions_file_centered_on_max_tss.main([self.truQuant_file, "20"])

        result_lines = stdout.getvalue().split("\n")
        result = [line.split() for line in result_lines if line]

        expected = [
            ["chr1", "925728", "925748", "SAMD11", "86", "+"],
            ["chr1", "959245", "959265", "NOC2L", "241", "-"]
        ]

        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()