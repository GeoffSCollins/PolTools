"""
If there is a poly A tail at the end of the sequence (at least 10 A's), convert those A's to N's
"""

import sys


final_reads = []


def print_usage():
    print("polyAToPolyN <Fastq File> <Minimum Number of A's>")


def parse_arguments():
    args = sys.argv[1:]

    if len(args) < 2:
        # Print an error
        sys.stderr.write("Not enough arguemnts")
        print_usage()
        sys.exit(1)

    fastq_file, min_as = args

    try:
        min_as = int(min_as)
    except:
        # Print an error
        sys.stderr.write("Minimum Number of A's must be an integer")
        sys.exit(1)

    return fastq_file, min_as


def convert_to_n(read_data, min_as):
    sequence = read_data[1].rstrip()
    if sequence[(min_as * -1):] == "A"*min_as:
        # The sequence ends with a poly(A) tail, so walk back until there is not an A
        curr_position = -10
        while (-1 * curr_position) - len(sequence) != 0 and sequence[curr_position] == 'A':
            curr_position -= 1

        # Now we can replace the read data's sequence with N's
        read_data[1] = read_data[1][:curr_position] + "N"*(-1*curr_position-1)
        return read_data

    return ""


def run_polyA_to_polyN(fastq_file, min_as):
    with open(fastq_file) as file:
        read_data = []
        for i, line in enumerate(file):
            if i % 4 == 0 and i != 0:
                # Do something with the data
                converted_data = convert_to_n(read_data, min_as)

                if converted_data != "":
                    final_reads.append(converted_data)

                read_data = []

            read_data.append(line)


def main():
    fastq_file, min_as = parse_arguments()
    run_polyA_to_polyN(fastq_file, min_as)


if __name__ == '__main__':
    main()
