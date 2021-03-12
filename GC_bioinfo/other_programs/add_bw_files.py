import sys
from GC_bioinfo.utils.make_random_filename import generate_random_filename

filenames = sys.argv[1:-1]
output_filename = sys.argv[-1]

if len(filenames) < 2:
    sys.stderr.write("Must provide at least two bw files")


if len(filenames) == 2:
    command = "bigwigCompare -b1 " + filenames[0] + " -b2 " + filenames[1] + " -bs 1 --operation add -o " + output_filename
    print(command)

else:
    commands = []
    intermediate_files = []
    intermediate_files.append(generate_random_filename('.bw'))
    commands.append("bigwigCompare -b1 " + filenames[0] + " -b2 " + filenames[1] + " -bs 1 --operation add -o " + intermediate_files[0])

    for filename in filenames[2:-1]:
        intermediate_files.append(generate_random_filename('.bw'))
        commands.append(
            "bigwigCompare -b1 " + filename + " -b2 " + intermediate_files[-2] + " -bs 1 --operation add -o " + intermediate_files[-1]
        )

    # Add the last file and output to the output_filename
    commands.append(
        "bigwigCompare -b1 " + filenames[-1] + " -b2 " + intermediate_files[-1] + " -bs 1 --operation add -o " + output_filename
    )

    for command in commands:
        print(command)

    # Print the remove command
    print("rm " + " ".join(intermediate_files))