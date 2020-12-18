import multiprocessing
import sys

from GC_bioinfo.utils.make_five_prime_bed_file import make_five_bed_file
from GC_bioinfo.utils.make_three_prime_bed_file import make_three_bed_file
from GC_bioinfo.utils.print_tab_delimited import print_tab_delimited
from GC_bioinfo.utils.remove_files import remove_files
from GC_bioinfo.utils.run_bedtools_coverage import run_coverage

tsrs_dict = {}


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
    print_tab_delimited(["Chromosome", "Left", "Right", "Name", "Strand"] + [tup[0].split("/")[-1] for tup in data])

    for tsr in data[0][1]:
        print_tab_delimited(list(tsr) + [tup[1][tsr] for tup in data])


def gather_data(read_type, seq_file, tsr_file):
    if read_type == "five":
        modified_seq_filename = make_five_bed_file(seq_file)
        need_to_remove_modified_seq_filename = True

    if read_type == "three":
        modified_seq_filename = make_three_bed_file(seq_file)
        need_to_remove_modified_seq_filename = True

    if read_type == "whole":
        modified_seq_filename = seq_file
        need_to_remove_modified_seq_filename = False


    # Then quantify the 5' ends
    coverage_file = run_coverage(tsr_file, modified_seq_filename)

    data = organize_counts(coverage_file, seq_file)

    remove_files(coverage_file)

    if need_to_remove_modified_seq_filename:
        remove_files(modified_seq_filename)

    return data


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo multicoverage <five/three/whole> <regions file> <sequencing files>\n")
    sys.stderr.write(
        "More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/multicoverage.rst\n")


def parse_args(args):
    if len(args) < 3:
        print_usage()
        sys.exit(1)

    read_type = sys.argv[1].lower()

    if read_type not in ["five", "three", "whole"]:
        print_usage()
        sys.exit(1)

    tsr_file = sys.argv[2]
    sequencing_files = sys.argv[3:]

    return read_type, tsr_file, sequencing_files


def main(args):
    read_type, tsr_file, sequencing_files = parse_args(args)

    pool = multiprocessing.Pool(processes=len(sequencing_files))
    data = pool.starmap(gather_data, [(read_type, seq_file, tsr_file) for seq_file in sequencing_files])
    output_data(data)


if __name__ == '__main__':
    main(sys.argv[1:])

