import os

from setuptools import setup, find_packages
from setuptools.command.install import install
from pathlib import Path


class CustomInstall(install):
    def run(self):
        # Move and source the tab completion
        dir_path = Path(__file__).absolute()
        completion_file = os.path.join(dir_path, 'GC_bioinfo-completion.bash')

        os.system('cp ' + completion_file + ' /etc/bash_completion.d/GC_bioinfo-completion.bash')
        os.system('source /etc/bash_completion.d/GC_bioinfo-completion.bash')

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
