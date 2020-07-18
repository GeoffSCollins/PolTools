#!/usr/bin/env bash

# Note! This file should be placed in /usr/bin/ or another appropriate place (next to GC_bioinfo

_bioinfo() {
 	local cur opts
 	COMPREPLY=()
 	cur="${COMP_WORDS[COMP_CWORD]}"
 	opts="base_distribution divergent_pileup_metaplot five_prime_metaplot inr_reads make_regions_file_centered_on_max_tss pausing_distance_distribution_from_maxTSS read_through_transcription sequence_searches three_prime_metaplot tps_distance_per_gene truQuant"

 	if [[ ${cur} == -* || ${COMP_CWORD} -eq 1 ]] ; then
 		COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
 		return 0
 	fi
}

complete -o default -F _bioinfo GC_bioinfo