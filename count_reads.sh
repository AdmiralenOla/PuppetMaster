#!/bin/bash

if [ $# -gt 0 ]; then
	for file in "$@"
	do
		echo "$file"
		echo $(zcat ${file} | wc -l)/4 | bc
	done
else
	echo "Counts reads in gzipped fastq. Usage: count_reads.sh file1.fastq.gz [file2.fastq.gz ...]"
fi
