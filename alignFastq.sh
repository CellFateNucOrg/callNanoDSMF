#! /bin/bash
## script to arrange all fastq files into batches for analysis. Combined files will be created for
## each barcode of interest for both passed and failed reads (e.g. pass_barcode01 or fail_barcode01).  

# Get variables from command line
expName=$1 	# experiment name (date of exp usually)
bc=$2 		# barcode
genomeFile=$3 	# full path to reference genome

# need absolute paths for nanopolish index. get it from the summary file.
summaryFile=`readlink -f ../sequencing_summary.txt`
expPath=`dirname $summaryFile`


####### modules to load ##########
module load vital-it
#module load R/3.5.1
module add UHTS/Analysis/minimap2/2.12;
module add UHTS/Analysis/samtools/1.8;
module add UHTS/Nanopore/nanopolish/0.10.2;
############################################

###################
## run MinIONQC ---> moved to basecalling script
##################
#
## create qc output directory
#qcDir=../fastqQC
#mkdir -p $qcDir
#
#echo "doing MinIONQC..."
## run Minion_qc
## source: https://github.com/roblanf/minion_qc
#
#MINIONQC=/home/pmeister/software/MinIONQC.R
#
#mkdir -p ${qcDir}/MinIONQC
#
#Rscript ${MINIONQC} -i ../sequencing_summary.txt -o $qcDir
#
#mv ../*.png ${qcDir}/MinIONQC/


##########################################################
# function to run Pauvre and nanoQC on individual batches
##########################################################

## activate python environment for QC programmes (Pauvre and NanoQC)
source activate albacore_env

# function for running Pauvre and NanoQC on each combined file
do_pauvre_nanoQC() {
   # get variables from arguments
   bc=$1 #barcode
   pf=$2 #pass or fail folder
   eN=$3 #expName
   qc=$4 #qcDir
   fq=../fastqFiles/${eN}_${pf}_${bc}.fastq.gz #input file to qc
   
   # assemble names for some of the output files
   outDir=${qc}/${pf}_${bc}
   mkdir -p $outDir
   outStats=${outDir}/pauvreStats.txt
   outPlot=pauvreMarginPlot.png
   
   # run QC programmes
   pauvre stats -f $fq > ${outStats}
   pauvre marginplot -f $fq  -o $outPlot
   mv pauvreMarginPlot.png $outDir  # move the png from current dir to qc output dir
   nanoQC $fq -o $outDir
}



################################################
# Collecting reads from barcodes that were used
################################################
echo "collecting reads from folder of barcodes that were used..."

# merge all reads from particular barcode into single file (pass fail separately)
mkdir -p ../fastqFiles
mkdir -p ../fastqQC

# need absolut paths for nanopolish index. get it from the summary file.
summaryFile=`readlink -f ../sequencing_summary.txt`
expPath=`dirname $summaryFile`
#expPath=/data/projects/p025/Jenny/20171027_Minion_TMP_Meth
echo $expPath

if [ -d ../workspace/pass/${bc} ];
then
    echo "passed reads..."
    cat ../workspace/pass/${bc}/*.fastq > ../fastqFiles/${expName}_pass_${bc}.fastq
    #nanopolish index -s ${summaryFile} -d ${expPath}/fast5files/ ${expPath}/fastqFiles/${expName}_pass_${bc}.fastq
    gzip ../fastqFiles/${expName}_pass_${bc}.fastq
    nanopolish index -s ${summaryFile} -d ${expPath}/fast5files/ ${expPath}/fastqFiles/${expName}_pass_${bc}.fastq.gz
    do_pauvre_nanoQC $bc pass $expName ../fastqQC 
fi

    if [ -d ../workspace/fail/${bc} ];
then
    echo "failed reads ..."
    cat ../workspace/fail/${bc}/*.fastq > ../fastqFiles/${expName}_fail_${bc}.fastq
    #nanopolish index -s ${summaryFile} -d ${expPath}/fast5files/ ${expPath}/fastqFiles/${expName}_fail_${bc}.fastq
    gzip ../fastqFiles/${expName}_fail_${bc}.fastq
    nanopolish index -s ${summaryFile} -d ${expPath}/fast5files/ ${expPath}/fastqFiles/${expName}_fail_${bc}.fastq.gz
    do_pauvre_nanoQC $bc fail $expName ../fastqQC
fi


################################################
# aligning to genome
################################################
echo "aligning to genome..."

mkdir -p ../bamFiles

# map reads to genome with minimap2
minimap2 -ax map-ont $genomeFile ../fastqFiles/${expName}_pass_${bc}.fastq.gz | samtools sort -T tmp -o ../bamFiles/${expName}_pass_${bc}.sorted.bam 
minimap2 -ax map-ont $genomeFile ../fastqFiles/${expName}_fail_${bc}.fastq.gz | samtools sort -T tmp -o ../bamFiles/${expName}_fail_${bc}.sorted.bam

echo "index bam file ..."
samtools index ../bamFiles/${expName}_pass_${bc}.sorted.bam
samtools index ../bamFiles/${expName}_fail_${bc}.sorted.bam

################################################
# identifying CmG
################################################
echo "identify CmG ..."

mkdir -p ../methylation_calls/CpG
nanopolish call-methylation -t 4 -q cpg -r ${expPath}/fastqFiles/${expName}_pass_${bc}.fastq.gz -b ${expPath}/bamFiles/${expName}_pass_${bc}.sorted.bam -g $genomeFile -w "ttTi5605:1-2378" > ${expPath}/methylation_calls/CpG/${expName}_pass_${bc}.tsv

nanopolish call-methylation -t 4 -q cpg -r ${expPath}/fastqFiles/${expName}_fail_${bc}.fastq.gz -b ${expPath}/bamFiles/${expName}_fail_${bc}.sorted.bam -g $genomeFile -w "ttTi5605:1-2378" > ${expPath}/methylation_calls/CpG/${expName}_fail_${bc}.tsv



################################################
# identifying GCm
################################################
echo "identify GCm ..."

mkdir -p ../methylation_calls/GpC
nanopolish call-methylation -t 4 -q gpc -r ${expPath}/fastqFiles/${expName}_pass_${bc}.fastq.gz -b ${expPath}/bamFiles/${expName}_pass_${bc}.sorted.bam -g $genomeFile -w "ttTi5605:1-2378" > ${expPath}/methylation_calls/GpC/${expName}_pass_${bc}.tsv

nanopolish call-methylation -t 4 -q gpc -r ${expPath}/fastqFiles/${expName}_fail_${bc}.fastq.gz -b ${expPath}/bamFiles/${expName}_fail_${bc}.sorted.bam -g $genomeFile -w "ttTi5605:1-2378" > ${expPath}/methylation_calls/GpC/${expName}_fail_${bc}.tsv
