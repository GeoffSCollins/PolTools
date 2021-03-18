import os

from setuptools import setup, find_packages
from setuptools.command.install import install
from pathlib import Path

from GC_bioinfo.utils.constants import tsr_finder_location

class CustomInstall(install):
    def run(self):
        # Copy and source the tab completion
        dir_path = Path(__file__).absolute()
        completion_file = os.path.join(dir_path, 'GC_bioinfo-completion.bash')

        os.system('cp ' + completion_file + ' /etc/bash_completion.d/GC_bioinfo-completion.bash')
        os.system('source /etc/bash_completion.d/GC_bioinfo-completion.bash')

        # Move GC_bioinfo.py to /usr/bin
        bin_file = os.path.join(dir_path, 'GC_bioinfo.py')
        os.system('cp ' + bin_file + ' /usr/bin')

        # Compile the tsrFinder file
        os.system('g++ ' + tsr_finder_location + '.cpp -o ' + tsr_finder_location)
        install.run(self)

setup(
    cmdclass={'install': CustomInstall},
    name='GC_bioinfo',
    author='Geoff Collins',
    version='1.0',
    packages=find_packages(),
    url='https://geoffscollins.github.io/GC_bioinfo/index.html',
    python_requires='>=3.5',
    project_urls= {
        'Documentation': 'https://geoffscollins.github.io/GC_bioinfo/index.html',
        'Source Code': 'https://github.com/GeoffSCollins/GC_bioinfo'
    }
)
