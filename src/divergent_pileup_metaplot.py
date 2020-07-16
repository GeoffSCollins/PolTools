import sys
import os
import multiprocessing

from utils.run_bedtools_coverage import run_coverage
from utils.remove_files import remove_files
from utils.make_random_filename import generate_random_filename
from utils.verify_bed_file import verify_bed_files
from utils.get_region_length import determine_region_length
from utils.output_metaplot_data import output_metaplot_data
from utils.check_dependencies import check_dependencies


def split_bed_file(bed_file):
    with open(bed_file.replace(".bed", "-FW.bed"), 'w') as fw_file:
        with open(bed_file.replace(".bed", "-RV.bed"), 'w') as rv_file:
            with open(bed_file) as file:
                for line in file:
                    if "+" in line:
                        fw_file.write(line)
                    else:
                        rv_file.write(line)


def get_pileups_helper(region_filename, sequencing_file, opposite_strand):
    split_bed_file(sequencing_file)

    fw_seq_file = sequencing_file.replace(".bed", "-FW.bed")
    rv_seq_file = sequencing_file.replace(".bed", "-RV.bed")

    fw_region_file = region_filename.replace(".bed", "-FW.bed")
    rv_region_file = region_filename.replace(".bed", "-RV.bed")

    # Run bedtools coverage on all of them
    fw_output = run_coverage(fw_region_file, fw_seq_file, flags=["-d"])
    rv_output = run_coverage(rv_region_file, rv_seq_file, flags=["-d"])

    # Combine the two files together
    combined_file = generate_random_filename()
    # os.system("cat " + fw_output + " " + rv_output + " > " + combined_file)
    os.system("cat " + rv_output + " > " + combined_file)

    # Counts list will start with the -500 nt and go to +500 nt
    counts_list = [0] * region_length

    with open(combined_file) as file:
        for line in file:
            if not opposite_strand:
                if "+" in line:
                    position = int(line.split()[-2]) - 1  # Subtract 1 because position starts at 1
                else:
                    # position = (region_length - 1) - int(line.split()[-2])  # Subtract 1 because position starts at 1
                    position = region_length - int(line.split()[-2])

                counts_list[position] += int(line.split()[-1])
            else:
                # If it is the opposite strand, we reverse the position and make the counts negative
                if "-" in line:
                    position = int(line.split()[-2]) - 1  # Subtract 1 because position starts at 1
                else:
                    # position = (region_length - 1) - int(line.split()[-2])  # Subtract 1 because position starts at 1
                    position = region_length - int(line.split()[-2])

                counts_list[position] -= int(line.split()[-1])

    remove_files(fw_output, fw_seq_file, rv_output, rv_seq_file, combined_file)

    return counts_list


def make_rev_region_file(region_filename):
    # We are wanting to get the opposite strand as well
    rev_region_filename = region_filename.replace(".bed", "rev_strands.bed")

    # Make a new version of the regions file
    with open(region_filename) as file:
        with open(rev_region_filename, 'w') as outfile:
            for line in file:
                if "+" in line:
                    outfile.write(line.replace("+", "-"))
                else:
                    outfile.write(line.replace("-", "+"))

    return rev_region_filename


def get_averages(input_2d_list):
    # This will take in a 2d list containing the regions and counts at that specific base
    # and output a list of the averages at each base
    averages_list = [0] * region_length

    # We loop through each base in the region
    for current_base_position, counts in enumerate(input_2d_list[0]):
        current_sum = 0

        # The total is the number of regions
        total = len(input_2d_list)

        # Loop through all of the regions
        for region in input_2d_list:
            current_sum += region[current_base_position]

        avg = current_sum / total
        averages_list[current_base_position] = avg

    return averages_list


def get_pileups(region_filename, sequencing_file, rev_region_filename):
    pileups = get_pileups_helper(region_filename, sequencing_file, False)
    rev_pileups = get_pileups_helper(rev_region_filename, sequencing_file, True)

    same_averages = []
    rev_averages = []

    return list(zip(pileups, rev_pileups)), sequencing_file


def parse_args(args):
    if len(args) < 2:
        print_usage()
        sys.exit(1)

    regions_filename = args[0]
    sequencing_files_list = args[1:]
    verify_bed_files(regions_filename, sequencing_files_list)

    global region_length
    region_length = determine_region_length(regions_filename)

    if region_length % 2 != 0:
        print("The region length is not even, so the +1 nt position cannot be determined. Exiting ...")
        sys.exit(1)

    return regions_filename, sequencing_files_list


def print_usage():
    print("Usage: ")
    print("python3 divergent_pileup_metaplot.py <Regions Filename> <Sequencing Files>")
    print("More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/divergent_pileup_metaplot.rst")


def run_divergent_pileup_plots(regions_filename, sequencing_files_list):
    # Make the opposite stranded regions file
    rev_region_filename = make_rev_region_file(regions_filename)

    # Split the files by strand
    split_bed_file(regions_filename)
    split_bed_file(rev_region_filename)

    pool = multiprocessing.Pool(processes=len(sequencing_files_list))
    averages = pool.starmap(get_pileups, [(regions_filename, sequencing_file, rev_region_filename) for sequencing_file in sequencing_files_list])

    # Remove all the stranded region files
    remove_files(regions_filename.replace(".bed", "-FW.bed"), regions_filename.replace(".bed", "-RV.bed"))
    remove_files(rev_region_filename.replace(".bed", "-FW.bed"), rev_region_filename.replace(".bed", "-RV.bed"), rev_region_filename)

    output_metaplot_data(averages, region_length)


def main(args):
    """
    python3 divergent_pileup_metaplot.py <Regions Filename> <Sequencing Files>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/divergent_pileup_metaplot.rst

    :param args: arguments provided to the program
    :type args: list
    :return:
    """
    check_dependencies("bedtools")
    regions_filename, sequencing_files_list = parse_args(args)
    run_divergent_pileup_plots(regions_filename, sequencing_files_list)


if __name__ == '__main__':
    main(sys.argv[1:])
