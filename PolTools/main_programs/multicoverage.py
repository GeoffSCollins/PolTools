import multiprocessing
import sys
import argparse

from PolTools.utils.make_read_end_file import make_read_end_file
from PolTools.utils.print_tab_delimited import print_tab_delimited
from PolTools.utils.remove_files import remove_files
from PolTools.utils.bedtools_utils.run_bedtools_coverage import run_coverage


def organize_counts(coverage_file, seq_file):
    ret_dict = {}
    with open(coverage_file) as file:
        for line in file:
            chrom, left, right, name, score, strand, counts, *_ = line.split()

            tsr = (chrom, left, right, name, strand)
            ret_dict[tsr] = counts

    return seq_file, ret_dict


def output_data(data):
    # Print the header first
    filenames = [filename.split("/")[-1] for filename, _ in data]
    print_tab_delimited(["Chromosome", "Left", "Right", "Name", "Strand"] + filenames)

    for _, tsr in data[0][1]:
        print_tab_delimited(list(tsr) + [file_data[tsr] for _, file_data in data])


def gather_data(read_type, seq_file, regions_file):

    if read_type in ['five', 'three']:
        modified_seq_filename = make_read_end_file(seq_file, read_type)
        need_to_remove_modified_seq_filename = True
    else:
        modified_seq_filename = seq_file
        need_to_remove_modified_seq_filename = False

    # Then quantify the 5' ends
    coverage_file = run_coverage(regions_file, modified_seq_filename)

    data = organize_counts(coverage_file, seq_file)

    remove_files(coverage_file)

    if need_to_remove_modified_seq_filename:
        remove_files(modified_seq_filename)

    return data


def parse_args(args):
    def positive_int(num):
        try:
            val = int(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")

        return val

    parser = argparse.ArgumentParser(prog='PolTools multicoverage',
        description='Quantify multiple sequencing files using the same regions\n' +
                     "More information can be found at " +
                     "https://geoffscollins.github.io/PolTools/multicoverage.html")

    parser.add_argument('read_type', metavar='read type', type=str, choices=["five", "three", "whole"],
                        help='either five, three, or whole')

    parser.add_argument('regions_filename', metavar='regions_filename', type=str,
                        help='Bed formatted regions file with an even region length or a region length of one.')

    parser.add_argument('seq_files', metavar='sequencing_files', nargs='+', type=str,
                        help='Bed formatted files from the sequencing experiment')

    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?',
                        default=multiprocessing.cpu_count())

    args = parser.parse_args(args)
    read_type = args.read_type
    regions_file = args.regions_filename
    sequencing_files = args.seq_files
    max_threads = args.threads

    return read_type, regions_file, sequencing_files, max_threads


def main(args):
    read_type, regions_file, sequencing_files, max_threads = parse_args(args)

    with multiprocessing.Pool(processes=max_threads) as pool:
        data = pool.starmap(gather_data, [(read_type, seq_file, regions_file) for seq_file in sequencing_files])

    output_data(data)


if __name__ == '__main__':
    main(sys.argv[1:])

