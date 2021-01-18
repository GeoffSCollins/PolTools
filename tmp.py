
import sys
import os

from collections import defaultdict



from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.remove_files import remove_files

output_file = "output.txt"

seq_files = sys.argv[1:]

# Key of chromosome and value of dict with key of filename and value of counts
output_dict = defaultdict(lambda : defaultdict(int))

for seq_file in seq_files:
    tmp_output_file = generate_random_filename(".tmp")

    os.system(
        "cut -f1 " + seq_file + " | sort | uniq -c > " + tmp_output_file
    )

    with open(tmp_output_file) as file:
        for line in file:
            counts, chromosome = line.split()
            output_dict[chromosome][seq_file] = counts

    remove_files(tmp_output_file)


with open(output_file, 'w') as file:
    # Print the headers
    file.write("chromosome" + "\t".join([seq_file for seq_file in seq_files]) + "\n")

    for chromosome in output_dict:
        output_line = ""

        for seq_file in seq_files:
            if seq_file in output_dict[chromosome]:
                output_line += "\t" + output_dict[chromosome][seq_file]
            else:
                output_line += "\t" + "0"

        file.write(chromosome + "\t" + output_line + "\n")
