import os
import sys
import datetime
from pathlib import Path


from PolTools.utils.remove_files import remove_files


def print_usage():
    sys.stderr.write("PolTools <program> <args>\n")


def parse_args(cli_args):
    # They need to provide the program and the arguments for the program

    if len(cli_args) == 0:
        # This means that the user needs prompting on what programs can be run
        print_usage()
        sys.exit(1)

    program = cli_args[0]
    command = program + ".py " + ' '.join(cli_args[1:])

    return program, command


def verify_program_exists(program):
    if program not in programs_list:
        sys.stderr.write(program + " is not in the program list. Exiting ...\n")
        sys.exit(1)


def run_program(command):
    _file_path = str(Path(__file__).parent.absolute())
    os.system("python3 " + _file_path + "/main_programs/" + command)


def _clean():
    temp_files = ["/tmp/" + file for file in os.listdir("/tmp") if file.startswith("PolTools")]

    for file in temp_files:
        curr_file = Path(file)

        last_modified_time = datetime.datetime.fromtimestamp(
            curr_file.stat().st_mtime
        )

        current_time = datetime.datetime.now()

        minutes_old = (current_time - last_modified_time).total_seconds() / 60.0

        # If the file is older than 45 minutes, then delete the file
        if minutes_old > 45:
            remove_files(file)


def _build():
    # First get the hg38 fasta file
    dir_path = Path(__file__).parent.absolute()
    static_path = os.path.join(dir_path, 'static/')

    if os.geteuid() != 0:
        sys.stderr.write("Install should be as a root user\n")
        sys.exit(1)

    if not os.path.exists(static_path + "hg38.fa"):
        os.system("wget -P " + static_path + " http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz")
        os.system("gunzip " + static_path + "hg38.fa.gz -c > " + static_path + "hg38.fa")

    # Add a cron function to clean regularly?
    response = input("Do you want to schedule the automatic temporary file cleaner to run hourly? (y/n)\n")
    while response not in ['y', 'n']:
        response = input("Do you want to schedule the automatic temporary file cleaner to run hourly? (y/n)\n")

    if response == 'y':
        with open("/etc/cron.hourly/PolTools_clean", 'w') as file:
            file.write("#!/bin/bash\n")
            file.write("sudo PolTools clean\n")

        # Add the exectuable flag to the file
        os.system("chmod +x /etc/cron.hourly/PolTools_clean")

    sys.stderr.write("Finished!\n")


def main(cli_args):
    program, command = parse_args(cli_args)

    if program == 'clean':
        _clean()
    elif program == 'build':
        _build()
    else:
        verify_program_exists(program)
        run_program(command)


if __name__ == '__main__':
    programs_list = ["base_distribution", "inr_reads", "make_regions_file_centered_on_max_tss",
                     "pausing_distance_distribution_from_maxTSS", "read_through_transcription", "sequence_searches",
                     "tps_distance_per_gene", "truQuant", "multicoverage", "tsrFinder", "gene_body_fold_change_heatmap",
                     "gene_body_heatmap", "quantify_gene_body_fold_change_heatmap", "TES_heatmap",
                     "TES_fold_change_heatmap", "nucleotide_heatmap", "track_links_from_bw",
                     "sequence_from_region_around_max_tss", "region_heatmap", "region_fold_change_heatmap", "metaplot"]

    main(sys.argv[1:])
