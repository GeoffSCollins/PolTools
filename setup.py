import os

from setuptools import setup, find_packages
from setuptools.command.install import install
from pathlib import Path

from PolTools.utils.constants import tsr_finder_location

class CustomInstall(install):
    def run(self):
        # Copy and source the tab completion
        dir_path = Path(__file__).parent.absolute()
        completion_file = os.path.join(dir_path, 'PolTools-completion.bash')

        os.system('cp ' + completion_file + ' /etc/bash_completion.d/PolTools-completion.bash')

        # Write a PolTools file in /usr/bin
        cli_location =  str(os.path.join(dir_path, 'PolTools/cli.py'))
        with open("/usr/bin/PolTools", 'w') as file:
            file.write(
                "#!/usr/bin/python3\n" +
                "import os\n" +
                "import sys\n\n\n" +


                "if __name__ == '__main__':\n" +
                "\tos.system(\n" +
                """\t\t"python3 """ + cli_location + """ " + ' '.join(sys.argv[1:])\n""" +
                "\t)\n"
            )

        os.system('chmod +x /usr/bin/PolTools')

        # Compile the tsrFinder file
        os.system('g++ ' + tsr_finder_location + '.cpp -o ' + tsr_finder_location)
        install.run(self)

setup(
    cmdclass={'install': CustomInstall},
    name='PolTools',
    author='Geoff Collins',
    version='1.0',
    packages=find_packages(),
    url='https://geoffscollins.github.io/PolTools/index.html',
    python_requires='>=3.5',
    project_urls= {
        'Documentation': 'https://geoffscollins.github.io/PolTools/index.html',
        'Source Code': 'https://github.com/GeoffSCollins/PolTools'
    }
)
