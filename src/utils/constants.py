from pathlib import Path

_file_path = str(Path(__file__).parent.parent.absolute())

annotation_file = _file_path + "/static/longest_transcript_with_downstream_start_codon.txt"
rna_blacklist_file = _file_path + "/static/hg38.GencodeV27.miRNA-rRNA-scRNA-snRNA-snoRNA-rRNA-scaRNA-tRNA.padded50bp.lsu_ssu.padded50bp.tRNAscan.padded50bp.bed"
hg38_chrom_sizes_file = _file_path + "/static/hg38.chrom.sizes.GC"