import sys

track_hub_text = """#############

track composite_name
compositeTrack on
type bigWig 0 1000
priority 9.5
visibility full
shortLabel composite_name
longLabel composite_name

track fw_name
parent composite_name
type bigWig 0 1000
visibility full
bigDataUrl fw_url
longLabel fw_name
shortLabel fw_name
negateValues off
color color_value
altColor color_value
alwaysZero on
autoScale on
maxHeightPixels 128:128:16

track rv_name
parent composite_name
type bigWig 0 1000
visibility full
bigDataUrl rv_url
longLabel rv_name
shortLabel rv_name
negateValues on
color color_value
altColor color_value
alwaysZero on
autoScale on
maxHeightPixels 128:128:16

"""

def output_hub_data(output_location, sequencing_files):
    for file in sequencing_files:

        # We first check which format the file is in
        if "-dedup-FW-FJ616285.1.bw" in file:
            ending_to_replace = "-dedup-FW-FJ616285.1.bw"
        elif "-dedup-FJ616285.1-FW.bw" in file:
            ending_to_replace = "-dedup-FJ616285.1-FW.bw"
        else:
            continue


        output_data = track_hub_text.replace("composite_name", file.replace(ending_to_replace, ""))

        output_data = output_data.replace("fw_name", file.replace(ending_to_replace, "_FW"))
        output_data = output_data.replace("rv_name", file.replace(ending_to_replace, "_RV"))

        output_data = output_data.replace("fw_url", "http://pricenas.biochem.uiowa.edu/" + output_location + "/" + file)
        output_data = output_data.replace("rv_url", "http://pricenas.biochem.uiowa.edu/" + output_location + "/" + file.replace("FW", "RV"))

        if "flavo" in file.lower():
            color_value = "153,58,29"
        else:
            color_value = "113,35,124"

        output_data = output_data.replace("color_value", color_value)

        print(output_data)


def print_usage():
    print("python3 bigwigs_to_track_hubs <output location on PriceNAS> <bigwigs (can use * operator)>")


def parse_args(args):
    if len(args) < 2:
        print_usage()
        sys.exit(1)

    output_location = args[0]
    sequencing_files = args[1:]
    fw_seq_files = []
    for file in sequencing_files:
        if "FW" in file and "RV" not in file:
            fw_seq_files.append(file)

    return output_location, fw_seq_files


def main(args):
    output_location, fw_seq_files = parse_args(args)
    output_hub_data(output_location, fw_seq_files)


if __name__ == '__main__':
    main(sys.argv[1:])
