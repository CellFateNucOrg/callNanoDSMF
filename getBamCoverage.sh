#! /bin/bash

mkdir -p ../aln
mkdir -p ../bigwig

samtools merge ../aln/dSMF16np_N2gw.bam ../bamFilesdb/20190411_pass_barcode01.sorted.bam ../bamFilesdb/20190411_pass_barcode05.sorted.bam 
samtools index ../aln/dSMF16np_N2gw.bam

samtools merge ../aln/dSMF20np_N2gw.bam ../bamFilesdb/20190411_pass_barcode02.sorted.bam ../bamFilesdb/20190411_pass_barcode03.sorted.bam ../bamFilesdb/20190411_pass_barcode04.sorted.bam 
samtools index ../aln/dSMF20np_N2gw.bam


allBams=(`ls ../aln/*.bam`)

for bamFile in "${allBams[@]}";
do	
	echo $bamFile
	fileBase=`basename $bamFile`
	bwFile=${fileBase%bam}cov.bw
	echo $bwFile
	bamCoverage -b $bamFile -o ../bigwig/${bwFile} --binSize 10  --numberOfProcessors 4
done

