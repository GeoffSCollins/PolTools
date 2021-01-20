import multiprocessing

import unittest.mock

import GC_bioinfo.main_programs.truQuant as truQuant

from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files

from quiter import Quieter

class TestTruQuant(unittest.TestCase):

    def test_arguments(self):
        with self.assertRaises(SystemExit):
            with Quieter():
                truQuant.parse_input([])

        seq_file_for_annotation = generate_random_filename()

        with open(seq_file_for_annotation, 'w') as file:
            file.write(
                "\t".join(['chr1', '1', '3', 'name', '0', '+'])
            )

        max_threads = multiprocessing.cpu_count()

        result = truQuant.parse_input([seq_file_for_annotation])
        self.assertEqual(result, ([seq_file_for_annotation], 1000, 0.3, 75, max_threads))

        result = truQuant.parse_input([seq_file_for_annotation, '-a', '500', '-b' '0.6', '-r', '10', '-t', '20'])
        self.assertEqual(result, ([seq_file_for_annotation], 500, 0.6, 10, 20))

        result = truQuant.parse_input([seq_file_for_annotation, '-a', '500', '--blacklisting_percent' '0.6', '-r', '10', '-t', '20'])
        self.assertEqual(result, ([seq_file_for_annotation], 500, 0.6, 10, 20))

        # Try making the blacklist percent greater than one
        with self.assertRaises(SystemExit):
            with Quieter():
                truQuant.parse_input([seq_file_for_annotation, '-b', '1.5'])

        remove_files(seq_file_for_annotation)


    def test_make_search_regions(self):
        # Tests both positive and negative strands with differing extensions

        regions_filename = generate_random_filename()

        with open(regions_filename, 'w') as file:
            file.write("\t".join(["chr1", "100", "200", "+", "positive_strand_test", "108", "111"]) + "\n")
            file.write("\t".join(["chr1", "2", "505", "-", "negative_strand_test", "498", "501"]) + "\n")

        search_regions_dict, annotations_dict = truQuant.make_search_regions(regions_filename, 10)

        expected_search_regions_dict = {
            "chr1": [["chr1", "90", "108", "positive_strand_test", "0", "+"],
                     ["chr1", "501", "515", "negative_strand_test", "0", "-"]
                     ]
        }

        expected_annotations_dict = {
            "positive_strand_test": ["chr1", "100", "200", "positive_strand_test", "0", "+"],
            "negative_strand_test": ["chr1", "2", "505", "negative_strand_test", "0", "-"]
        }

        self.assertDictEqual(search_regions_dict, expected_search_regions_dict)
        self.assertDictEqual(annotations_dict, expected_annotations_dict)

        search_regions_dict, annotations_dict = truQuant.make_search_regions(regions_filename, 100)

        expected_search_regions_dict = {
            "chr1": [["chr1", "0", "108", "positive_strand_test", "0", "+"],
                     ["chr1", "501", "605", "negative_strand_test", "0", "-"]
                     ]
        }

        expected_annotations_dict = {
            "positive_strand_test": ["chr1", "100", "200", "positive_strand_test", "0", "+"],
            "negative_strand_test": ["chr1", "2", "505", "negative_strand_test", "0", "-"]
        }

        self.assertDictEqual(search_regions_dict, expected_search_regions_dict)
        self.assertDictEqual(annotations_dict, expected_annotations_dict)

        remove_files(regions_filename)

    def test_map_tsrs_to_search_regions(self):
        # Test # TSRs
        # One that is contained in the search region, one with partial overlap in the 5' end,
        # one with partial overlap in the 3' end. One with no overlap before the TSR. One with no overlap after the TSR
        # One on the opposite strand

        # Need to define a TSR file and a search regions dict
        search_regions_dict = {
            "chr1": [["chr1", "100", "200", "positive_strand_test", "0", "+"],
                     ["chr1", "500", "600", "negative_strand_test", "0", "-"]
                     ]
        }

        tsr_filename = generate_random_filename(".tab")

        additional_columns = ["tss_left", "tss_right", "tss_strength", "avg_tss"]

        with open(tsr_filename, 'w') as file:
            file.write("\t".join(["chr1", "40", "60", "no_overlap", "0", "+"] + additional_columns) + "\n")
            file.write("\t".join(["chr1", "300", "320", "no_overlap2", "0", "+"] + additional_columns) + "\n")
            file.write("\t".join(["chr1", "90", "110", "partial_overlap_5'", "0", "+"] + additional_columns) + "\n")
            file.write("\t".join(["chr1", "190", "210", "partial_overlap_3'", "0", "+"] + additional_columns) + "\n")
            file.write("\t".join(["chr1", "140", "160", "complete_overlap", "0", "+"] + additional_columns) + "\n")
            file.write("\t".join(["chr1", "140", "160", "opposite_strand", "0", "-"] + additional_columns) + "\n")

        gene_tsr_dict, flow_through_tsrs = truQuant.map_tsrs_to_search_regions(tsr_filename, search_regions_dict)

        expected_gene_tsr_dict = {
            "positive_strand_test": [
                ["chr1", "90", "110", "partial_overlap_5'", "0", "+", "avg_tss"],
                ["chr1", "190", "210", "partial_overlap_3'", "0", "+", "avg_tss"],
                ["chr1", "140", "160", "complete_overlap", "0", "+", "avg_tss"]
            ]
        }

        expected_flow_through_tsrs = [
            ["chr1", "40", "60", "no_overlap", "0", "+", "avg_tss"],
            ["chr1", "300", "320", "no_overlap2", "0", "+", "avg_tss"],
            ["chr1", "140", "160", "opposite_strand", "0", "-", "avg_tss"]
        ]

        self.assertDictEqual(gene_tsr_dict, expected_gene_tsr_dict)
        self.assertEqual(flow_through_tsrs, expected_flow_through_tsrs)

        remove_files(tsr_filename)

    def test_find_max_tsr_in_search_region(self):
        # Need to define a gene_tsr_dict first
        gene_tsr_dict = {
            "gene_name": [
                ["chr1", "20", "40", "read_sum", "30", "+", "avg_tss"],
                ["chr1", "40", "60", "read_sum", "90", "+", "avg_tss"],
                ["chr1", "60", "80", "read_sum", "18", "+", "avg_tss"],
                ["chr1", "120", "140", "read_sum", "90", "+", "avg_tss"]
            ]
        }

        max_tsrs_dict, non_max_tsrs_dict = truQuant.find_max_tsr_in_search_region(gene_tsr_dict)

        expected_max_tsrs_dict = {
            "gene_name": ["chr1", "40", "60", "read_sum", "90", "+", "avg_tss"]
        }

        expected_non_max_tsrs_dict = {
            "gene_name": [
                ["chr1", "20", "40", "read_sum", "30", "+", "avg_tss"],
                ["chr1", "60", "80", "read_sum", "18", "+", "avg_tss"],
                ["chr1", "120", "140", "read_sum", "90", "+", "avg_tss"]
            ]
        }

        self.assertDictEqual(max_tsrs_dict, expected_max_tsrs_dict)
        self.assertDictEqual(non_max_tsrs_dict, expected_non_max_tsrs_dict)

    def test_define_pause_regions_and_gene_bodies(self):
        # Need to define output filenames, max_tsrs_dict, annotations_dict, and pause_region_radius

        pause_region_filename = generate_random_filename()
        gene_body_filename = generate_random_filename()

        region_filenames = pause_region_filename, gene_body_filename

        max_tsrs_dict = {
            "positive_gene": ["chr1", "140", "160", "read_sum", "90", "+", "155.2134"],
            "negative_gene": ["chr1", "200", "220", "read_sum", "532", "-", "208.3854"]
        }

        annotations_dict = {
            "positive_gene": ["chr1", "20", "900", "positive_gene", "0", "+"],
            "negative_gene": ["chr1", "85", "750", "negative_gene", "0", "-"]
        }

        pause_region_radius = 75

        truQuant_regions_dict = truQuant.define_pause_regions_and_gene_bodies(region_filenames,
                                                                              max_tsrs_dict,
                                                                              annotations_dict,
                                                                              pause_region_radius)

        pause_regions = []
        with open(pause_region_filename) as file:
            for line in file:
                pause_regions.append(line.split())

        gene_bodies = []
        with open(gene_body_filename) as file:
            for line in file:
                gene_bodies.append(line.split())

        expected_pause_regions = [
            ["chr1", "80", "230", "positive_gene", "90", "+"],
            ["chr1", "134", "284", "negative_gene", "532", "-"]
        ]

        expected_gene_bodies = [
            ["chr1", "230", "900", "positive_gene", "90", "+"],
            ["chr1", "85", "134", "negative_gene", "532", "-"]
        ]

        expected_truQuant_regions_dict = {
            "positive_gene": {
                "Pause": ["chr1", 80, 230, "+"],
                "Body": [230, 900, 670]
            },
            "negative_gene": {
                "Pause": ["chr1", 134, 284, "-"],
                "Body": [85, 134, 49]
            }
        }

        self.assertEqual(pause_regions, expected_pause_regions)
        self.assertEqual(gene_bodies, expected_gene_bodies)
        self.assertDictEqual(truQuant_regions_dict, expected_truQuant_regions_dict)

        # Verify the pause region and gene body files
        remove_files(region_filenames)

    def test_map_flow_through_tsrs(self):

        annotations_dict = {
            "positive_gene": ["chr1", "200", "1000", "positive_gene", "0", "+"]
        }

        flow_through_tsrs = [
            # [tsr_chromosome, tsr_left, tsr_right, tsr_read_sum, tsr_strength, tsr_strand, avg_tss]
            ["chr1", "150", "250", "partial_overlap_5'", "0", "+", "153.532"],
            ["chr1", "950", "1050", "partial_overlap_3'", "0", "+", "1000.213"],
            ["chr1", "80", "180", "no_overlap_5'", "0", "+", "175.924"],
            ["chr1", "1200", "1300", "no_overlap_3'", "0", "+", "1243.438"],
            ["chr1", "400", "500", "complete_overlap", "0", "+", "450.000"],
            ["chr1", "400", "500", "opposite_strand", "0", "-", "450.000"],
        ]

        mapped_flow_through_tsrs_dict = truQuant.map_flow_through_tsrs(annotations_dict, flow_through_tsrs)

        expected_mapped_flow_through_tsrs_dict = {
            # [tsr_chromosome, tsr_left, tsr_right, gene_name, tsr_counts, tsr_strand]
            "positive_gene": [
                ["chr1", "150", "250", "positive_gene", "0", "+"],
                ["chr1", "950", "1050", "positive_gene", "0", "+"],
                ["chr1", "400", "500", "positive_gene", "0", "+"]
            ]

        }

        self.assertDictEqual(mapped_flow_through_tsrs_dict, expected_mapped_flow_through_tsrs_dict)

    def test_make_blacklisted_regions(self):

        blacklist_filename = generate_random_filename()

        max_tsrs_dict = {
            "gene_name": ["chr1", "40", "60", "read_sum", "90", "+", "avg_tss"]
        }

        non_max_tsrs_dict = {
            "gene_name": [
                ["chr1", "90", "110", "read_sum", "45", "+", "avg_tss"],
                ["chr1", "70", "90", "read_sum", "20", "+", "avg_tss"]
            ]
        }

        mapped_flow_through_tsrs_dict = {
            "gene_name": [
                ["chr1", "450", "470", "read_sum", "30", "+"],
                ["chr1", "900", "920", "read_sum", "1", "+"]
            ]
        }

        mapped_tsrs = max_tsrs_dict, non_max_tsrs_dict, mapped_flow_through_tsrs_dict

        percent_for_blacklisting = 0.3

        truQuant.make_blacklisted_regions(blacklist_filename, mapped_tsrs, percent_for_blacklisting)

        blacklist = []
        with open(blacklist_filename) as file:
            for line in file:
                blacklist.append(line.split())

        expected_blacklist = [
            ["chr1", "450", "470", "gene_name", "30", "+"],
            ["chr1", "90", "110", "gene_name", "45", "+"]
        ]

        self.assertEqual(blacklist, expected_blacklist)

        remove_files(blacklist_filename)

    def test_get_counts(self):

        pause_regions_filename = generate_random_filename()

        with open(pause_regions_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "100", "250", "positive_gene", "90", "+"]) + "\n" +
                "\t".join(["chr1", "700", "850", "negative_gene", "1523", "-"]) + "\n"
            )

        gene_body_filename = generate_random_filename()

        with open(gene_body_filename, 'w') as file:
            file.write(
                "\t".join(["chr1", "251", "750", "positive_gene", "90", "+"]) + "\n" +
                "\t".join(["chr1", "200", "700", "negative_gene", "1523", "-"]) + "\n"
            )

        blacklisted_sequencing_file = generate_random_filename()

        with open(blacklisted_sequencing_file, 'w') as file:
            file.write(
                "\t".join(["chr1", "80", "220", "5'not_counted", "0", "+"]) + "\n" +
                "\t".join(["chr1", "259", "285", "5'not_counted", "0", "+"]) + "\n" +
                "\t".join(["chr1", "132", "220", "5'count", "0", "+"]) + "\n" +

                "\t".join(["chr1", "132", "800", "5'count", "0", "-"]) + "\n" +
                "\t".join(["chr1", "750", "783", "5'count", "0", "-"]) + "\n" +
                "\t".join(["chr1", "750", "900", "5'not_counted", "0", "-"]) + "\n" +
                "\t".join(["chr1", "500", "600", "5'not_counted", "0", "-"]) + "\n"
            )

        indv_gene_counts_dict = truQuant.get_counts_in_paused_region(pause_regions_filename,
                                                                     blacklisted_sequencing_file)

        indv_gene_counts_dict = truQuant.get_counts_in_gene_bodies(gene_body_filename,
                                                                   blacklisted_sequencing_file,
                                                                   indv_gene_counts_dict)

        expected_indv_gene_counts_dict = {
            "positive_gene": {
                "Pause": 1,
                "Body": 1
            },
            "negative_gene": {
                "Pause": 2,
                "Body": 1
            }
        }

        self.assertDictEqual(indv_gene_counts_dict, expected_indv_gene_counts_dict)

        remove_files(pause_regions_filename, gene_body_filename, blacklisted_sequencing_file)


if __name__ == '__main__':
    unittest.main()
