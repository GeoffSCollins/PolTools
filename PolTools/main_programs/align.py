import os
import sys
import subprocess
import argparse

from PolTools.utils.make_random_filename import generate_random_filename
from PolTools.utils.remove_files import remove_files

from multiprocessing import Pool
from collections import defaultdict
from statistics import mean
from shutil import rmtree

# TODO: Move the files to pricenas

global MAIN_DIRECTORY


# Get input
def parse_args(args):
    def positive_int(num):
        try:
            val = int(num)
            if val <= 0:
                raise Exception("Go to the except")
        except:
            raise argparse.ArgumentTypeError(num + " must be positive")

        return val

    parser = argparse.ArgumentParser(prog='PolTools align',
                                     description='Align raw data to generate bed and bigwig files.\n' +
                                                 "More information can be found at " +
                                                 "https://geoffscollins.github.io/PolTools/align.html")

    parser.add_argument('-u', '--url', dest='url', metavar='url', type=str)
    parser.add_argument('-f', '--folder-name', dest='folder_name', metavar='folder_name', type=str)
    parser.add_argument('-l', '--umi-length', dest='umi_length', metavar='umi_length', type=positive_int, default=8)
    parser.add_argument('--min-insert', dest='min_insert', metavar='min_insert', type=positive_int, default=32)
    parser.add_argument('--max-insert', dest='max_insert', metavar='max_insert', type=positive_int, default=600)
    parser.add_argument('-t', '--threads', dest='threads', metavar='threads', type=positive_int, nargs='?', default=10)

    parser.add_argument('-i', '--index', dest='index', metavar='index', type=str)
    parser.add_argument('-c', '--chrom-sizes-file', dest='chrom_sizes_file', metavar='chrom_sizes_file', type=str)
    parser.add_argument('-s', '--server-name', dest='server_name', metavar='server_name', type=str, default=None)
    parser.add_argument('--spike-in-genome', dest='spike_in_genome', metavar='spike_in_genome', type=str, default=None)
    parser.add_argument('-g', '--genome', dest='genomes', metavar='genome', type=str, nargs="+", action="extend")
    parser.add_argument('-x', '--conversion-file', dest='sample_name_conversion_filename', metavar='conversion_file', type=str, default=None)

    return parser.parse_args(args)


def read_sample_name_conversion_table(filename):
    if not filename or not os.path.isfile(filename):
        return None, ['.fastq.gz']

    downloading_prefixes = []

    sample_name_conversion = {}
    with open(filename) as file:
        for i, line in enumerate(file):
            if i != 0:
                download_prefix, prefix, new_name, group_name = line.split()
                downloading_prefixes.append(download_prefix)

                sample_name_conversion[prefix] = [new_name, group_name]

    return sample_name_conversion, downloading_prefixes


def build_directory_structure():
    os.mkdir(MAIN_DIRECTORY)
    os.mkdir(MAIN_DIRECTORY + "original_files")
    os.mkdir(MAIN_DIRECTORY + "trimmed")
    os.mkdir(MAIN_DIRECTORY + "sam")
    os.mkdir(MAIN_DIRECTORY + "bed")
    os.mkdir(MAIN_DIRECTORY + "bigwig")
    os.mkdir(MAIN_DIRECTORY + "logs")


def download_data(url, prefixes):
    # Download the data in the original_files folder
    # -q quiet, -nH is no host directories, -P is the directory prefix
    original_files_directory = MAIN_DIRECTORY + "original_files"
    # Comma separated list to -A
    prefixes = [prefix + "*.fastq.gz" for prefix in prefixes]
    accepted_file_string = ",".join(prefixes)
    os.system("wget -q --cut-dirs=100 -nH -r -A " + accepted_file_string + " -P " + original_files_directory + " " + url)


def _trim_helper(prefix):
    # Do the trimming and write the stdout and stderr to the logging file
    result = subprocess.check_output(
        "trim_galore --paired --small_rna --dont_gzip --quality 0 --length 34 --output_dir " +
        MAIN_DIRECTORY + "trimmed" + " " + MAIN_DIRECTORY + "original_files/" + prefix + "*",
        shell=True,
        stderr=subprocess.STDOUT,
        text=True
    )

    logs_directory = MAIN_DIRECTORY + "logs"

    with open(logs_directory + "/trimming.log", 'a') as file:
        file.write(result)


def trim_adapters(threads):
    # Group the samples by the prefix -- string until the first underscore
    prefixes = set()
    for filename in os.listdir(MAIN_DIRECTORY + "original_files"):
        prefixes.add(filename.split("_lane")[0])

    # Trim adapters for each group with the same prefix
    with Pool(threads) as pool:
        pool.map(_trim_helper, prefixes)


def _align_and_dedup(read_one_file, read_two_file, sample_name, align_params, threads):
    trim5 = align_params['trim5']
    trim3 = align_params['trim3']
    minins = align_params['minins']
    maxins = align_params['maxins']
    index = align_params['index']

    sam_directory = MAIN_DIRECTORY + "sam"
    bed_directory = MAIN_DIRECTORY + "bed"

    sam_file = sam_directory + "/" + sample_name + ".sam"

    result = subprocess.check_output(
        "bowtie --trim5 " + trim5 + " --trim3 " + trim3 + " --minins " + minins + " --maxins " + maxins +
        " --fr --best --allow-contain --sam --fullref --chunkmbs 5000 --threads " + str(threads) + " " + index +
        " -1 " + read_one_file +
        " -2 " + read_two_file +
        " > " + sam_file,
        shell=True,
        stderr=subprocess.STDOUT,
        text=True
    )

    logs_directory = MAIN_DIRECTORY + "logs"

    with open(logs_directory + "/align.log", 'a') as file:
        file.write(sample_name + "\n")
        file.write(result + "\n\n")

    # Also dedup and log
    result = subprocess.check_output(
        "dedup -n " + trim5 + " -f -p -t " + str(threads) +
        " -1 " + read_one_file +
        " -2 " + read_two_file +
        " -S " + sam_file,
        shell=True,
        stderr=subprocess.STDOUT,
        text=True
    )

    with open(logs_directory + "/dedup.log", 'a') as file:
        file.write(result + "\n\n")

    # Move the bed file to the bed directory
    os.system('mv ' + sam_file.replace(".sam", "-dedup.bed") + ' ' + bed_directory)


def align(align_params, sample_name_conversion, threads):
    trimmed_directory = MAIN_DIRECTORY + "trimmed"
    bed_directory = MAIN_DIRECTORY + "bed"

    # Group the trimmed files by prefix
    prefixes = defaultdict(list)
    for filename in os.listdir(trimmed_directory):
        if 'trimming_report' not in filename:
            prefix = '_'.join(filename.split('_')[:-6])
            # Remove the lane string if there is one
            if 'lane' in prefix:
                location = prefix.find('lane')
                prefix = prefix[:location] + prefix[location + 6:]

            prefixes[prefix].append(filename)

    # If there is more than one lane with the prefix, then we align and dedup separately but combine the bed files at the end
    for prefix, filenames in prefixes.items():
        if len(filenames) == 2:
            prefixes[prefix] = [sorted(filenames)]
        else:
            # Group by lane
            lane_strings = defaultdict(list)

            for filename in filenames:
                lane_string = filename.split('_')[-3]
                lane_strings[lane_string].append(filename)

            new_list = []

            # Verify that there are two files for each lane string and make sure we add them to the prefixes dict
            for lane_filenames in lane_strings.values():
                if len(lane_filenames) != 2:
                    print("Something odd happened with lanes. Here are the samples grouped: " + str(lane_filenames))
                else:
                    new_list.append(lane_filenames)

            prefixes[prefix] = sorted(new_list)

    for prefix, groups in prefixes.items():
        need_to_combine_bed_files = False if len(groups) == 1 else True

        # If the sample name is being converted, do it
        if sample_name_conversion and prefix in sample_name_conversion:
            sample_name = sample_name_conversion[prefix][0]
        else:
            sample_name = prefix

        lane_specific_bed_files = []

        for i, group in enumerate(groups):
            read_one_file, read_two_file = group

            read_one_file = trimmed_directory + "/" + read_one_file
            read_two_file = trimmed_directory + "/" + read_two_file

            # If the lanes need to be combined, then we will add a number to the end of the sample name
            if need_to_combine_bed_files:
                lane_specific_bed_files.append(bed_directory + "/" + sample_name + str(i) + ".bed")
                _align_and_dedup(read_one_file, read_two_file, sample_name + str(i), align_params, threads)
            else:
                _align_and_dedup(read_one_file, read_two_file, sample_name, align_params, threads)

        if need_to_combine_bed_files:
            # Combine the bed files that begin with the sample name and write a new file without the numbers at the end
            os.system(
                'cat ' + bed_directory + "/" + sample_name + "*.bed" +
                ' > ' + bed_directory + "/" + sample_name + ".bed"
            )

            # Remove all the lane specific files
            remove_files(lane_specific_bed_files)


def _calculate_spike_in_normalization_table(sample_name_conversion, spike_in_genome, genomes):
    def get_group_name(filename):
        sample_name = filename.replace('-dedup.bed', '')

        for prefix in sample_name_conversion:
            if sample_name_conversion[prefix][0] == sample_name:
                return sample_name_conversion[prefix][1]

        return "unknown"


    bed_directory = MAIN_DIRECTORY + "bed"

    normalization_group_data = defaultdict(dict)

    # Loop through each bed file
    for filename in os.listdir(bed_directory):
        # Get the total number of reads and split it for each genome
        group_name = get_group_name(filename)

        normalization_group_data[group_name][filename] = {
            'total reads': 0,
            'library size normalization factor': 0.0,
            'spike in normalization factor': 0.0
        }

        for genome in genomes:
            normalization_group_data[group_name][filename][genome] = 0

        with open(bed_directory + "/" + filename) as file:
            for line in file:
                chrom, left, right, name, score, strand = line.split()

                # Add the read to total counts and the specific genome
                normalization_group_data[group_name][filename]['total reads'] += 1

                if 'chr' in chrom:
                    normalization_group_data[group_name][filename]['hg38'] += 1

                else:
                    normalization_group_data[group_name][filename][chrom] += 1


    # For each group, calculate the normalization factors and write it to the file
    with open(MAIN_DIRECTORY + "normalization_table.tsv", 'w') as file:
        file.write("\t".join(
            ["Sample", "Total Reads"] + genomes + ["Library Size Correction Factor", "Spike-in Correction Factor",
                                                   "Final Correction Factor", "Group"]
        ) + "\n")

        for group_name, group_data in normalization_group_data.items():

            # Now calculate the normalization factors
            average_counts = mean([sample_data['total reads'] for sample_data in group_data.values()])

            for sample_data in group_data.values():
                sample_data['library size normalization factor'] = average_counts / sample_data['total reads']

            average_normalized_spike_in_reads = mean(
                [sample_data['library size normalization factor'] * sample_data[spike_in_genome]
                 for sample_data in group_data.values()]
            )

            for sample_data in group_data.values():
                sample_data['spike in normalization factor'] = average_normalized_spike_in_reads / sample_data[spike_in_genome]

            # Write the data for each sample
            for sample_filename, sample_data in group_data.items():
                final_correction_factor = sample_data['library size normalization factor'] * \
                                          sample_data['spike in normalization factor']
                file.write(
                    "\t".join(
                        [sample_filename.replace("-dedup.bed", ""), str(sample_data['total reads'])] +
                        [str(sample_data[genome]) for genome in genomes] +
                        [str(sample_data['library size normalization factor'])] +
                        [str(sample_data['spike in normalization factor'])] +
                        [str(final_correction_factor), group_name]
                    ) + "\n"
                )


def _calculate_library_size_normalization_table(sample_name_conversion):
    def get_group_name(filename):
        sample_name = filename.replace('-dedup.bed', '')

        for prefix in sample_name_conversion:
            if sample_name_conversion[prefix][0] == sample_name:
                return sample_name_conversion[prefix][1]

        return "unknown"


    bed_directory = MAIN_DIRECTORY + "bed"

    normalization_group_data = defaultdict(dict)

    for filename in os.listdir(bed_directory):
        group_name = get_group_name(filename)

        with open(bed_directory + "/" + filename) as file:
            for i, _ in enumerate(file):
                pass

        normalization_group_data[group_name][filename] = i + 1

    # Write the data to a file
    with open(MAIN_DIRECTORY + "normalization_table.tsv", 'w') as file:
        for group_name, group_data in normalization_group_data.items():

            average_library_size = mean([count for count in group_data.values()])

            # Write the headers
            file.write(
                "\t".join(
                    ["Sample", "Total Reads", "Correction Factor"]
                ) + "\n"
            )

            for filename, counts in group_data.items():
                file.write(
                    "\t".join(
                        [filename.replace("-dedup.bed", ""), str(counts), str(average_library_size / counts), group_name]
                    ) + "\n"
                )


def calculate_normalization_table(sample_name_conversion, spike_in_genome=None, genomes=None):
    # If they provide arguments, then it must be spike_in normalization
    if spike_in_genome and genomes:
        _calculate_spike_in_normalization_table(sample_name_conversion, spike_in_genome, genomes)
    else:
        # Calculate based on total library size
        _calculate_library_size_normalization_table(sample_name_conversion)


def _generate_bigwig(filename, normalization_factor, chrom_sizes_file, hosting_location):
    bigwig_directory = MAIN_DIRECTORY + "bigwig"

    bedgraph_files = []
    links = []

    link_template = "track type=bigWig visibility=full name='<name>' autoScale=on alwaysZero=on " + \
                    "windowingFunction=maximum negateValues=off color=0,0,0 altColor=0,0,0 " + \
                    "bigDataUrl=<hostLocation>/<file>"

    tmp_file = generate_random_filename(".tmp")

    os.system('cat ' + filename + """ | awk '$1 != "chrEBV" {print $0}' > """ + tmp_file)
    os.system('bedSort ' + tmp_file + ' ' + filename)

    for strand, strand_string in [("+", "FW"), ("-", "RV")]:
        os.system(
            'bedtools genomecov -scale ' + normalization_factor +
            ' -i ' + filename + ' -g ' + chrom_sizes_file +
            ' -bg -strand ' + strand +
            ' > ' + filename.replace('.bed', '-' + strand_string + '.bedGraph')
        )
        os.system(
            'bedtools genomecov -scale ' + normalization_factor +
            ' -i ' + filename + ' -g ' + chrom_sizes_file +
            ' -bg -strand ' + strand + ' -5' +
            ' > ' + filename.replace('.bed', '-' + strand_string + '-5.bedGraph')
        )
        os.system(
            'bedtools genomecov -scale ' + normalization_factor +
            ' -i ' + filename + ' -g ' + chrom_sizes_file +
            ' -bg -strand ' + strand + ' -3' +
            ' > ' + filename.replace('.bed', '-' + strand_string + '-3.bedGraph')
        )

        bedgraph_files.append(filename.replace('.bed', '-' + strand_string + '.bedGraph'))
        bedgraph_files.append(filename.replace('.bed', '-' + strand_string + '-5.bedGraph'))
        bedgraph_files.append(filename.replace('.bed', '-' + strand_string + '-3.bedGraph'))

    # Convert bedGraphs to bigwigs and add track links
    for file in bedgraph_files:
        os.system('bedSort ' + file + ' ' + file)
        bigwig_file = bigwig_directory + "/" + file.split('/')[-1].replace('.bedGraph', '.bw')
        os.system('bedGraphToBigWig ' + file + ' ' + chrom_sizes_file + ' ' + bigwig_file)

        curr_link = link_template

        if "-RV.bw" in bigwig_file or "-RV-5.bw" in bigwig_file or "-RV-3.bw" in bigwig_file:
            curr_link = curr_link.replace("negateValues=off", "negateValues=on")

        file_basename = file.split("/")[-1]

        # Get the sample names
        sample_name = file_basename
        sample_name = sample_name.replace(".bw", "").replace("_", " ").replace("-", " ")

        # Now put the sample name in
        curr_link = curr_link.replace("<name>", sample_name).replace("<file>", file_basename)

        if hosting_location:
            curr_link = curr_link.replace("<hostLocation>", hosting_location)

        links.append(curr_link)

    # Remove all the bedgraph and tmp files
    remove_files(bedgraph_files, tmp_file)

    # Sort the links so they are in the order of FW, RV, FW5, RV 5, FW3, RV3
    links.sort()
    final_links = [
        links[2],  # Fw
        links[5],  # Rv
        links[1],  # Fw 5
        links[4],  # Rv 5
        links[0],  # Fw 3
        links[3],  # Rv 3
    ]

    return final_links


def generate_bigwigs(chrom_sizes_filename, threads, hosting_location):
    # First get the normalization factors
    normalization_factors = {}

    with open(MAIN_DIRECTORY + "normalization_table.tsv") as file:
        for i, line in enumerate(file):
            if i != 0:
                sample_name, *_, normalization_factor, _ = line.split()
                normalization_factors[sample_name] = normalization_factor

    bed_directory = MAIN_DIRECTORY + "bed"
    bigwig_directory = MAIN_DIRECTORY + "bigwig"

    args = []
    for filename in os.listdir(bed_directory):
        sample_name = filename.split("/")[-1].replace("-dedup.bed", "")

        args.append(
            [
                bed_directory + "/" + filename,
                normalization_factors[sample_name],
                chrom_sizes_filename,
                hosting_location
            ]
        )

    with Pool(threads) as pool:
        all_links = pool.starmap(_generate_bigwig, args)

    # Write the links to a file
    with open(bigwig_directory + "/links.txt", 'w') as file:
        for dataset_links in all_links:
            for link in dataset_links:
                file.write(link + "\n")


def clean_directories():
    # Remove the trimmed and sam directories
    rmtree(MAIN_DIRECTORY + "trimmed")
    rmtree(MAIN_DIRECTORY + "sam")


def main(args):
    args = parse_args(args)
    url = args.url
    folder_name = args.folder_name
    threads = args.threads
    chrom_sizes_file = args.chrom_sizes_file

    server_name = args.server_name
    spike_in_genome = args.spike_in_genome
    genomes = args.genomes + [spike_in_genome]
    sample_name_conversion_filename = args.sample_name_conversion_filename

    align_params = {
        'trim5': str(args.umi_length),
        'trim3': str(args.umi_length),
        'minins': str(args.min_insert),
        'maxins': str(args.max_insert),
        'index': args.index
    }

    hosting_location = server_name + folder_name if server_name else ''

    global MAIN_DIRECTORY
    MAIN_DIRECTORY = "/home/geoff/.trial/" + folder_name + "/"

    # If given a conversion file, then only download the ones in the conversion file
    sample_name_conversion, prefixes = read_sample_name_conversion_table(sample_name_conversion_filename)

    build_directory_structure()
    download_data(url, prefixes)
    trim_adapters(threads)
    align(align_params, sample_name_conversion, threads)
    calculate_normalization_table(sample_name_conversion, spike_in_genome, genomes)
    generate_bigwigs(chrom_sizes_file, threads, hosting_location)
    clean_directories()


if __name__ == '__main__':
    main(sys.argv[1:])
