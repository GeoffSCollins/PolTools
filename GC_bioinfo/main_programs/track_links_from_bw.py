import sys


def generate_links(bigwig_files):

    link_template = "track type=bigWig visibility=full name='<name>' autoScale=on alwaysZero=on " + \
                    "windowingFunction=maximum negateValues=off color=0,0,0 altColor=0,0,0 bigDataUrl=<server>/<file>"
    links = []

    for file in sorted(bigwig_files):
        curr_link = link_template

        # If the strand is negative, negate the values
        if "-RV.bw" in file or "-RV-5.bw" in file or "-RV-3.bw" in file:
            curr_link = curr_link.replace("negateValues=off", "negateValues=on")

        file_basename = file.split("/")[-1]

        # Get the sample names
        sample_name = file_basename
        sample_name = sample_name.replace(".bw", "")

        # Replace the underscores and hyphens in the sample name with spaces
        sample_name = sample_name.replace("_", " ").replace("-", " ")

        # Now put the sample name in
        curr_link = curr_link.replace("<name>", sample_name)

        # Put in the true filename
        curr_link = curr_link.replace("<file>", file_basename)

        links.append(curr_link)

    return links


def sort_links(links):
    # We want the primary sort to be the same but then we want the following order
    # Fw, Rv, Fw 5', Rv 5', Fw 3', Rv 3'
    sorted_links_list = []

    # This means the links come in sets of 6 if ends are made or 2 if no ends are made
    use_groups_of_six = False

    # First check if there are primes by getting a group of six
    if len(links) >= 6:
        trial_group = links[:6]

        # If there exists a file with the ending of -FW-5.bw, then the groups of 6 will be used

        for file in trial_group:
            if "FW-5.bw" in file:
                use_groups_of_six = True
                break


    # Now we know how big of groups to use, so we go through each group in the links and sort those
    if use_groups_of_six:
        groups = []

        for i in range(0, len(links), 6):
            groups.append(links[i:i+6])

        # The current order is Fw 3, Fw 5, Fw, Rv 3, Rv 5, Rv
        # We swap to get the desired order
        for group in groups:
            new_order = [
                group[2], # Fw
                group[5], # Rv
                group[1], # Fw 5
                group[4], # Rv 5
                group[0], # Fw 3
                group[3], # Rv 3
            ]
            sorted_links_list.extend(new_order)
    # If they are in groups of two (only fw and rv), we can just keep the sort as it goes fw then rv
    else:
        sorted_links_list = links

    return sorted_links_list


def print_usage():
    sys.stderr.write("Usage: \n")
    sys.stderr.write("GC_bioinfo track_links_from_bw <bigwig files>\n")
    sys.stderr.write(
        "More information can be found at https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/track_links_from_bw.rst\n")


def parse_args(args):
    if len(args) == 0:
        sys.stderr.write("At least one bigwig file must be provided\n")
        print_usage()
        sys.exit(1)

    return args


def main(args):
    bigwig_files = parse_args(args)
    links = generate_links(bigwig_files)
    sorted_links = sort_links(links)

    # Output the links
    for link in sorted_links:
        print(link)


if __name__ == '__main__':
    main(sys.argv[1:])
