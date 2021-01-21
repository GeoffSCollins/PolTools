import sys

from GC_bioinfo.utils.run_metaplot import run_metaplot, parse_input


def main(args):
    """
    GC_bioinfo three_prime_metaplot <Regions Filename> <Sequencing Files>
    More information can be found at https://geoffscollins.github.io/GC_bioinfo/three_prime_metaplot.html

    :param args: arguments provided to the program
    :type args: list
    :return:
    """

    regions_filename, sequencing_files_list, region_length, max_threads = parse_input(args, "three")
    run_metaplot(regions_filename, sequencing_files_list, region_length, "three", max_threads)

if __name__ == '__main__':
    main(sys.argv[1:])