#!/bin/bash
seq=$1
cat raw_fastq/${seq}_R1all_trim.fastq raw_fastq/${seq}_R2all_trim.fastq | seqkit grep -f <(sed -n '1~4{p;n;n;n}' MATA_bowtie_raw/${seq}.fastq | sed 's/^@//' | sed 's/\/[1|2]$//' | sort | uniq) > MATA_bowtie_raw/${seq}_search.fastq
