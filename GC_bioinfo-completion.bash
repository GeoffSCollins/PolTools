#!/usr/bin/env bash

# Note! This file should be placed in /usr/bin/ or another appropriate place (next to GC_bioinfo)

_bioinfo() {
  local cur opts
 	COMPREPLY=()
 	cur="${COMP_WORDS[COMP_CWORD]}"
 	opts="base_distribution metaplot make_regions_file_centered_on_max_tss \
 	pausing_distance_distribution_from_maxTSS read_through_transcription sequence_searches \
 	tps_distance_per_gene truQuant polyAToPolyN tsrFinder gene_body_fold_change_heatmap \
 	gene_body_heatmap quantify_gene_body_fold_change_heatmap TES_heatmap TES_fold_change_heatmap \
 	multicoverage nucleotide_heatmap track_links_from_bw sequence_from_region_around_max_tss region_heatmap \
 	region_fold_change_heatmap"

 	if [[ ${cur} == -* || ${COMP_CWORD} -eq 1 ]] ; then
 	  COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
 	  return 0
 	fi
}

complete -o default -F _bioinfo GC_bioinfo