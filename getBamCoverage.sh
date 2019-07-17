#! /bin/bash


mkdir -p ../bigwig

for bamFile in "$@";
do	
	echo $bamFile
	fileBase=`basename $bamFile`
	bwFile=${fileBase%sorted.bam}cov.bw
	echo $bwFile
	bamCoverage -b $bamFile -o ../bigwig/${bwFile} --binSize 10  --numberOfProcessors 4
done

