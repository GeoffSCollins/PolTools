import sys

from GC_bioinfo.utils.run_metaplot import run_metaplot, parse_input


def main(args):
    """
    GC_bioinfo three_prime_metaplot <Regions Filename> <Sequencing Files>
    More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/three_prime_metaplot.rst

    :param args: arguments provided to the program
    :type args: list
    :return:
    """

    regions_filename, sequencing_files_list, region_length = parse_input(args, "three")
    run_metaplot(regions_filename, sequencing_files_list, region_length, "three")

if __name__ == '__main__':
    main(sys.argv[1:])