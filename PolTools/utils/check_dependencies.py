import subprocess
import sys
from shutil import which

import pkg_resources


class ToolNotFoundException(Exception):
    def __init__(self, tool):
        super().__init__(tool)
        self.tool = tool

    def __str__(self):
        return self.tool + " is not installed."


def _cli_tool_exists(tool):
    """
    Checks if CLI tool exists. Throws ToolNotFoundException if tool does not exist.

    :param tool: name of the tool the program will use
    :type tool: str
    :exception ToolNotFoundException: if CLI tool does not exist
    """

    if which(tool) == None:
        raise ToolNotFoundException(tool)


def _python_module_exits(module):
    """
    Returns true if the module is installed in python3

    :param module: name of the module to check
    :type module: str
    :rtype: bool
    """
    try:
        dist = pkg_resources.get_distribution(module)
    except pkg_resources.DistributionNotFound:
        install_module_response = input(module + " is not installed. Do you want to install it? y or n\n")
        if install_module_response == "y":
            # Install it
            print("Installing " + module)
            python = sys.executable
            subprocess.check_call([python, '-m', 'pip', 'install', module], stdout=subprocess.DEVNULL)

        elif install_module_response == "n":
            # Don't install and exit
            print("The program needs " + module + " to run. Exiting ...")
            sys.exit(1)
        else:
            # Invalid response and exit
            print("Invalid response. Exiting ...")
            sys.exit(1)


def check_dependencies(*program_names):
    """
    Checks if the system has bedtools, biopython, and tsrFinder
    :return:
    """
    # Check for bedtools and tsrFinder
    for program_name in program_names:
        if program_name == "bedtools":
            _cli_tool_exists("bedtools")

        if program_name == "tsrFinderPARALLEL":
            _cli_tool_exists("tsrFinderPARALLEL")

        if program_name == "tsrFinderM1I":
            _cli_tool_exists("tsrFinderM1I")

        if program_name == "biopython":
            _python_module_exits("biopython")
