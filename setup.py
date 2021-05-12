import os

from setuptools import setup, find_packages
from setuptools.command.install import install
from pathlib import Path

from PolTools.utils.constants import tsr_finder_location


class CustomInstall(install):
    def run(self):
        # Copy the tab completion
        dir_path = Path(__file__).parent.absolute()
        completion_file = os.path.join(dir_path, 'PolTools-completion.bash')

        os.system('cp ' + completion_file + ' /etc/bash_completion.d/PolTools-completion.bash')

        # Compile the tsrFinder file
        os.system('g++ ' + tsr_finder_location + '.cpp -o ' + tsr_finder_location)
        install.run(self)


requirements_file = os.path.join(Path(__file__).parent.absolute(), 'requirements.txt')

with open(str(requirements_file)) as file:
    reqs = file.readlines()

setup(
    cmdclass={'install': CustomInstall},
    name='PolTools',
    author='Geoff Collins',
    version='1.0',
    packages=find_packages(exclude=('tests', 'PolTools/other_programs')),
    url='https://geoffscollins.github.io/PolTools/index.html',
    python_requires='>=3.5',
    project_urls={
        'Documentation': 'https://geoffscollins.github.io/PolTools/index.html',
        'Source Code': 'https://github.com/GeoffSCollins/PolTools'
    },
    install_requires=reqs,
    entry_points={
        'console_scripts': [
            'PolTools = PolTools.cli:main'
        ]
    }
)
