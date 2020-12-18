

from GC_bioinfo.utils.constants import hg38_chrom_sizes_file
from GC_bioinfo.utils.make_random_filename import generate_random_filename
from GC_bioinfo.utils.run_bedtools_getfasta import run_getfasta
from GC_bioinfo.utils.remove_files import remove_files

from collections import defaultdict

def get_nt_content_dict(nt, truQuant_output_file, content_distance):
    # Need to make a dictionary with the keys of the gene names and values of the nt content for the first distance nts

    # Get the chrom sizes
    chrom_sizes_dict = defaultdict(int)
    with open(hg38_chrom_sizes_file) as file:
        for line in file:
            chromosome, size = line.split()
            chrom_sizes_dict[chromosome] = int(size)

    # Read the genes and make the positions
    positions = defaultdict(list)

    with open(truQuant_output_file) as file:
        for i, line in enumerate(file):
            if i != 0:
                gene_name, chromosome, pause_left, pause_right, strand, total_reads, max_tss, max_tss_five_prime_reads, avg_tss, \
                std_tss, gene_body_left, gene_body_right, *_ = line.split()

                if strand == "+":
                    if int(max_tss) + content_distance <= chrom_sizes_dict[chromosome] and int(max_tss) >= 0:
                        positions[gene_name] = [chromosome, max_tss, str(int(max_tss) + content_distance), gene_name, "0", strand]
                else:
                    # If the strand is negative, we shift the right 1
                    if int(max_tss) - content_distance >= 0 and int(max_tss) + 1 <= chrom_sizes_dict[chromosome]:
                        positions[gene_name] = [chromosome, str(int(max_tss) - content_distance + 1), str(int(max_tss) + 1), gene_name, "0", strand]

    # Now we get the sequences
    sequences = defaultdict(str)

    regions_file = generate_random_filename()

    gene_names = []

    with open(regions_file, 'w') as file:
        for gene_name in positions:
            gene_names.append(gene_name)
            file.write("\t".join(positions[gene_name]) + "\n")

    fasta_file = run_getfasta(regions_file)

    with open(fasta_file) as file:
        for i, line in enumerate(file):
            if i % 2 == 1:
                sequence = line.rstrip().upper()

                curr_gene = gene_names[int(i / 2)]

                sequences[curr_gene] = sequence.count(nt)

    remove_files(regions_file, fasta_file)

    return sequences


# Now 

