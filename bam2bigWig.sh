#!/usr/bin/env bash

## convert RNAseq BAMs to bigwigs

# need to point at directory containing only BAMs of interest
# this is from QUACKERS output, so they've been sorted
# may need to be tweaked depending on your use
COUNTER=1

num=$(ls *.bam|echo `wc -l`|awk '{print $1}')
echo There are ${num} bam files to convert

for file in *.bam; do
	echo Working on file ${COUNTER}...
	echo "Indexing file..."
	samtools index ${file}
	echo "Creating bigwig..."
	bamCoverage -b ${file} --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 5 --extendReads 200 -o ${file%%.*}.bigwig
	echo Finished converting file number ${COUNTER} to bigWig
	let COUNTER=COUNTER+1
done

