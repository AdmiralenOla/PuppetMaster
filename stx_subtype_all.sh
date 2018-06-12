#!/bin/bash

# Script for performing stx subtyping on lots of strains at once

basedir=$(pwd)

for x in $(ls *R1*.fastq.gz)
do
	R1=$x
	R2=${x/R1/"R2"}
	# Find contigs. Must be in current dir and named isolatename.fasta
	contigs=${x%%_*}.fasta
	echo "Typing ${x%%_*}"
	stx_subtyping dual_stx_subtyping $R1 $R2 $contigs stx_subtyping_results
done

# Collate results into single file

echo "Collating results"

header="Isolate\tMapping\tAssembly\n"
output=""

for y in $(ls *R1*.fastq.gz)
do
	ybase=${y%.fastq.gz}
	isolatename=${y%%_*}
	MapRes=$(cat stx_subtyping_results/${ybase}.MapSNP.txt)
	MapRes=${MapRes/$ybase/}
	MapRes="$(echo -e "${MapRes}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
	MapRes="$(echo -e "${MapRes}" | tr '[:blank:]' ',')"
	AssRes=$(cat stx_subtyping_results/${isolatename}.xml.AssBLAST.txt)
	AssRes=${AssRes/${isolatename}.xml/}
	AssRes="$(echo -e "${AssRes}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
	AssRes="$(echo -e "${AssRes}" | tr '[:blank:]' ',')"
	output="${output}${isolatename}\t${MapRes}\t${AssRes}\n"
done

output=$(echo "$output" | tr -d '[:blank:]')
echo -ne "$header$output" > stx_subtyping_results.txt
echo "Finished stx subtyping"
