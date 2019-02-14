#! /bin/bash
## script to arrange all fastq files into batches for analysis. Folders will be created for
## each barcode of interest for both passed and failed reads (e.g. pass_barcode01 or fail_barcode01).  
## Folders will also be created for all misclassified (barcode not in barcodesOfInterest) and 
## unclassified reads split between pass_unclass and fail_unclass. 

# Get variables from command line
expName=$1
bcOfInterest=$2

#expName="20171027"
#barcodesOfInterest=( barcode05 barcode06 barcode07 barcode08 )
#bcOfInterest=${barcodesOfInterest[0]}

####### modules to load ##########
module add UHTS/Analysis/minimap2/2.12;
module add UHTS/Analysis/samtools/1.8;
module add UHTS/Nanopore/naligning to genome
################################################nopolish/0.10.2;

##################
# run MinIONQC
#################

# create qc output directory
qcDir=../fastqQC
mkdir -p $qcDir

echo "doing MinIONQC..."
# run Minion_qc
# source: https://github.com/roblanf/minion_qc
module load vital-it
module load R/3.5.1

MINIONQC=/home/pmeister/software/MinIONQC.R

mkdir -p ${qcDir}/MinIONQC

Rscript ${MINIONQC} -i ../sequencing_summary.txt -o $qcDir

mv ../*.png ${qcDir}/MinIONQC/


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

# scroll through all barcodes you have used and take properly classified reads and move them
# to new folders
mkdir -p ../fastqFiles
for b in ${barcodesOfInterest[@]};
do
    if [ -d ../workspace/pass/${b} ];
    then
        cat ../workspace/pass/${b}/*.fastq > ../fastqFiles/${expName}_pass_${b}.fastq
	nanopolish index -d ../fast5files/*/fast5/ ../fastqFiles/${expName}_pass_${b}.fastq
        gzip ../fastqFiles/${expName}_pass_${b}.fastq
        do_pauvre_nanoQC $b pass $expName $qcDir 
    fi

        if [ -d ../workspace/fail/${b} ];
    then
        cat ../workspace/fail/${b}/*.fastq > ../fastqFiles/${expName}_fail_${b}.fastq
	nanopolish index -d fast5files/*/fast5/ ../fastqFiles/${expName}_pass_${b}.fastq
        gzip ../fastqFiles/${expName}_fail_${b}.fastq
        do_pauvre_nanoQC $b fail $expName $qcDir
    fi
done


################################################
# aligning to genome
################################################

mkdir -p ../bamFiles

# map reads to genome with minimap2
minimap2 -a -x map-ont $genomeFile ../fastqFiles/${expName}_pass_${bcOfInterest}.fastq.gz | samtools sort -T tmp -o ../bamFiles/${expName}_pass_${bcOfInterest}.sorted.bam 

samtools index ../bamFiles/${expName}_pass_${bcOfInterest}.sorted.bam

