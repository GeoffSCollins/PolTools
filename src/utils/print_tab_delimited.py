

def print_tab_delimited(input_list):
    """
    Prints out the tab delimited contents of the inputted list

    :parameter input_list: the list to be printed out
    :type input_list: list
    """
    # First convert all of the values to a string
    cleaned_list = [str(val) for val in input_list]

    # Now print out the tab delimited values
    print(*cleaned_list, sep='\t')