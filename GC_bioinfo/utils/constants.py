from pathlib import Path

_file_path = str(Path(__file__).parent.parent.absolute())

annotation_file = _file_path + "/static/longest_transcript_with_downstream_start_codon.txt"
rna_blacklist_file = _file_path + "/static/hg38.GencodeV27.miRNA-rRNA-scRNA-snRNA-snoRNA-rRNA-scaRNA-tRNA.padded50bp.lsu_ssu.bed"
hg38_chrom_sizes_file = _file_path + "/static/hg38.chrom.sizes.GC"
tsr_finder_location = _file_path + "/utils/tsrFinder"
generate_heatmap_location = _file_path + "/utils/generate_heatmap.R"