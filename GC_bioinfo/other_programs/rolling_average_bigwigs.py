import os
import sys

from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def output_line(curr_counts, bed_location, file_handles, prev_positions, chrom_sizes_dict):
    chromosome, curr_location_left, curr_location_right = bed_location
    a_file, t_file, g_file, c_file = file_handles

    for nt in ["A", "T", "G", "C"]:
        if nt == "A":
            file = a_file
        elif nt == "T":
            file = t_file
        elif nt == "G":
            file = g_file
        else:
            file = c_file

        old_bed_location = prev_positions[nt][0]
        old_chrom, old_left, old_right = old_bed_location

        # If the current position is larger than the chromosome size, just go until you hit the next chromosome
        if int(curr_location_right) > chrom_sizes_dict[chromosome]:
            continue

        # If the previous chromosome is the old location, just update to the new one and continue
        if old_chrom != chromosome:
            new_bed_location = (chromosome, curr_location_left, curr_location_right)
            prev_positions[nt] = (new_bed_location, curr_counts[nt])

        # If the previous position counts are 0 and the curr counts are not, we update prev_position and don't output to line
        if prev_positions[nt][1] == 0 and curr_counts[nt] != 0:
            new_bed_location = (chromosome, curr_location_left, curr_location_right)
            prev_positions[nt] = (new_bed_location, curr_counts[nt])
            continue

        # If they are both zero don't do anything
        if prev_positions[nt][1] == 0 and curr_counts[nt] == 0:
            continue

        # If the prev_counts are the same as the current counts, don't print the line and continue
        if prev_positions[nt][1] == curr_counts[nt]:
            continue

        # If the prev_counts are x and the current counts are not x, print the line and update prev position
        if prev_positions[nt][1] != curr_counts[nt]:
            output_counts = prev_positions[nt][1]

            file.write("\t".join(
                [old_chrom, str(old_left), str(curr_location_left), str(output_counts)]) + "\n")

            new_bed_location = (chromosome, curr_location_left, curr_location_right)
            prev_positions[nt] = (new_bed_location, curr_counts[nt])

    return prev_positions


def make_bedgraphs(fasta_file, window_size, chrom_size_dict):
    # We will do the sliding window technique with window size
    a_bedgraph = generate_random_filename(extension=".bedgraph")
    t_bedgraph = generate_random_filename(extension=".bedgraph")
    g_bedgraph = generate_random_filename(extension=".bedgraph")
    c_bedgraph = generate_random_filename(extension=".bedgraph")

    a_file = open(a_bedgraph, 'w')
    t_file = open(t_bedgraph, 'w')
    g_file = open(g_bedgraph, 'w')
    c_file = open(c_bedgraph, 'w')

    file_handles = [a_file, t_file, g_file, c_file]

    with open(fasta_file) as file:
        num_of_lines = sum(1 for _ in file)
        file.seek(0)

        curr_bases = []

        for i, line in enumerate(file):
            if i % 100_000 == 0 or i == num_of_lines - 1:
                printProgressBar(i + 1, num_of_lines, prefix='Progress:', suffix='Complete', length=50)

            if ">" in line:
                chromosome = line[1:].rstrip()
                initialized = False
                continue

            sequence = line.rstrip().upper()

            if len(curr_bases) + len(sequence) <= window_size:
                for nt in sequence:
                    curr_bases.append(nt)
                continue

            if not initialized:
                # First need to add the remaining nucleotides
                curr_position_in_sequence = 0
                while len(curr_bases) != window_size:
                    curr_bases.append(sequence[curr_position_in_sequence])
                    curr_position_in_sequence += 1

                initialized = True

                # First initialize the list to the first window_size bases
                curr_bases = list(sequence[:window_size])

                # Set the current location to the center of the window
                curr_location = int((window_size - 1) / 2)

                bed_location = [chromosome, curr_location, curr_location + 1]

                curr_counts = {
                    "A": curr_bases.count("A"),
                    "T": curr_bases.count("T"),
                    "G": curr_bases.count("G"),
                    "C": curr_bases.count("C"),
                }

                # Initialize prev_positions
                prev_positions = {
                    "A": (bed_location, curr_counts["A"]),
                    "T": (bed_location, curr_counts["T"]),
                    "G": (bed_location, curr_counts["G"]),
                    "C": (bed_location, curr_counts["C"]),
                }

            # Now we will do the sliding. This involves removing the first element then adding the next base
            for i in range(curr_position_in_sequence, len(sequence)):
                curr_location += 1
                removed_base = curr_bases.pop(0)
                curr_bases.append(sequence[i])

                # Update the curr_counts
                if removed_base in curr_counts:
                    curr_counts[removed_base] -= 1

                if sequence[i] in curr_counts:
                    curr_counts[sequence[i]] += 1

                bed_location = [chromosome, curr_location, curr_location + 1]

                prev_positions = output_line(curr_counts, bed_location, file_handles, prev_positions, chrom_size_dict)

            curr_position_in_sequence = 0

    a_file.close()
    t_file.close()
    g_file.close()
    c_file.close()

    bedgraphs = (
        a_bedgraph,
        t_bedgraph,
        g_bedgraph,
        c_bedgraph
    )

    return bedgraphs


def print_usage():
    sys.stderr.write("python3 rolling_average_bigwigs.py <window size> <chrom sizes file> <genome fasta file>\n")


def get_args(args):
    if len(args) != 3:
        print_usage()
        sys.exit(1)

    window_size, chrom_sizes_file, fasta_file = args

    try:
        window_size = int(window_size)
        if window_size < 1:
            sys.stderr.write("Window size must be positive. Exiting ...")
            sys.exit(1)
    except:
        sys.stderr.write("Window size must be an integer. Exiting ...")
        sys.exit(1)

    # Make sure the fasta file exists
    if not os.path.isfile(chrom_sizes_file):
        sys.stderr.write("Chromosome sizes file does not exist. Exiting ...")
        sys.exit(1)

    if not os.path.isfile(fasta_file):
        sys.stderr.write("Genome fasta file does not exist. Exiting ...")
        sys.exit(1)

    return window_size, chrom_sizes_file, fasta_file


def build_chrom_sizes_dict(chrom_sizes_file):
    chrom_sizes_dict = {}

    with open(chrom_sizes_file) as file:
        for line in file:
            chrom_name, size = line.split()

            chrom_sizes_dict[chrom_name] = int(size)

    return chrom_sizes_dict


def main(args):
    window_size, chrom_sizes_file, fasta_file = get_args(args)

    chrom_sizes_dict = build_chrom_sizes_dict(chrom_sizes_file)

    bedgraphs = make_bedgraphs(fasta_file, window_size, chrom_sizes_dict)

    a_bedgraph, t_bedgraph, g_bedgraph, c_bedgraph = bedgraphs

    # Run bedsort on all the files
    for bedgraph in bedgraphs:
        os.system("sort -k1,1 -k2,2n " + bedgraph + " > " + bedgraph + ".tmp")

    # Convert them to bigwig files
    os.system("bedGraphToBigWig " + " ".join([ a_bedgraph + ".tmp", chrom_sizes_file, "a_content_" + str(window_size) + "bp.bw"]))
    os.system("bedGraphToBigWig " + " ".join([ t_bedgraph + ".tmp", chrom_sizes_file, "t_content_" + str(window_size) + "bp.bw"]))
    os.system("bedGraphToBigWig " + " ".join([ g_bedgraph + ".tmp", chrom_sizes_file, "g_content_" + str(window_size) + "bp.bw"]))
    os.system("bedGraphToBigWig " + " ".join([ c_bedgraph + ".tmp", chrom_sizes_file, "c_content_" + str(window_size) + "bp.bw"]))

    remove_files(a_bedgraph, t_bedgraph, g_bedgraph, c_bedgraph)
    remove_files(a_bedgraph + ".tmp", t_bedgraph + ".tmp", g_bedgraph + ".tmp", c_bedgraph + ".tmp")


if __name__ == '__main__':
    main(sys.argv[1:])