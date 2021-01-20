import os
import sys
from pathlib import Path


def print_usage():
    sys.stderr.write("GC_bioinfo <program> <args>\n")


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


def main(cli_args):
    program, command = parse_args(cli_args)
    verify_program_exists(program)
    run_program(command)


if __name__ == '__main__':
    programs_list = ["base_distribution", "divergent_pileup_metaplot", "five_prime_metaplot", "inr_reads",
                     "make_regions_file_centered_on_max_tss", "pausing_distance_distribution_from_maxTSS",
                     "read_through_transcription", "sequence_searches", "three_prime_metaplot", "tps_distance_per_gene",
                     "truQuant", "multicoverage", "tsrFinder", "gene_body_combined_heatmap",
                     "gene_body_fold_change_heatmap", "gene_body_heatmap", "quantify_gene_body_fold_change_heatmap",
                     "TES_heatmap", "TES_combined_heatmap", "TES_fold_change_heatmap",
                     "nucleotide_heatmap", "track_links_from_bw", "sequence_from_region_around_max_tss"]

    main(sys.argv[1:])
