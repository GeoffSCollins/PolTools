#!/usr/bin/python3
#!/bin/bash
import os
import sys


# Note! This file should be in /usr/bin or similar location so it can direct the command to the PolTools installation location
if __name__ == '__main__':
    # Note! The location of cli.py needs to be changed if not in the Price lab!
    command = 'python3 PolTools/cli.py ' + ' '.join(sys.argv[1:])
    os.system(command)
