import sys
import os

# First get the hg38 fasta file
dir_path = os.getcwd()
static_path = os.path.join(dir_path, 'GC_bioinfo/static')

if os.geteuid() != 0:
    raise PermissionError("Install should be as a root user")

os.system("wget -P " + static_path + " http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz")
os.system("gunzip " + static_path + "hg38.fa.gz -c > " + static_path + "hg38.fa")

# Move and source the tab completion
sys.stderr.write("Adding tab completion...\n")
completion_file = dir_path + '/GC_bioinfo-completion.bash'
os.system('cp ' + completion_file + ' /etc/bash_completion.d/GC_bioinfo-completion.bash')

sys.stderr.write("Finished!\n")