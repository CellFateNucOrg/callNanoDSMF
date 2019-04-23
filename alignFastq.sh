#! /bin/bash
## script to arrange all fastq files into batches for analysis. Combined files will be created for
## each barcode of interest for both passed and failed reads (e.g. pass_barcode01 or fail_barcode01).  

# Get variables from command line
expName=$1 	# experiment name (date of exp usually)
bc=$2 		# barcode
genomeFile=$3 	# full path to reference genome

# need absolute paths for nanopolish index. get it from the summary file.
#summaryFile=`readlink -f ../fastqFiles/sequencing_summary.txt`
expPath=`readlink -f ../`
#expPath=/data/projects/p025/Jenny/20190411_dSMFv021-025np_N2gw
summaryFile=${expPath}/fastqFiles/sequencing_summary.txt

####### modules to load ##########
module load vital-it
#module load R/3.5.1
module add UHTS/Analysis/minimap2/2.12;
module add UHTS/Analysis/samtools/1.8;
############################################


##########################################################
# function to run Pauvre and nanoQC on individual batches
##########################################################

## activate python environment for QC programmes (Pauvre and NanoQC)
source activate albacore_env

 #th}/meth_calls/${expName}_pass_${bc}_GpCcalls_${chr[$i]}.tsv#function for running Pauvre and NanoQC on each combined file
do_pauvre_nanoQC() {
   # get variables from arguments
   bc=$1 #barcode
   pf=$2 #pass or fail folder
   eN=$3 #expName
   qc=$4 #qcDir
   fq=${expPath}/bcFastq/${eN}_${pf}_${bc}.fastq.gz #input file to qc
   
   # assemble names for some of the output files
   outDir=${qc}/${pf}_${bc}
   mkdir -p $outDir
   outStats=${outDir}/pauvreStats.txt
   
   # run QC programmes
   pauvre stats -f $fq > ${outStats}
   nanoQC $fq -o $outDir
}


#################################################
## Collecting reads from barcodes that were used
#################################################
#echo "collecting reads from folder of barcodes that were used..."
#
## merge all reads from particular barcode into single file (pass fail separately)
#echo ${expPath}/fast5Files
#
#if [ -d ${expPath}/bcFastq/pass/${bc} ];
#then
#    echo "passed reads..."
#    cat ${expPath}/bcFastq/pass/${bc}/*.fastq > ${expPath}/bcFastq/${expName}_pass_${bc}.fastq
#    gzip ${expPath}/bcFastq/${expName}_pass_${bc}.fastq
#    ${NANOPOLISH_DIR}/nanopolish index -s ${summaryFile} -d ${expPath}/fast5Files ${expPath}/bcFastq/${expName}_pass_${bc}.fastq.gz
#    do_pauvre_nanoQC $bc pass $expName ${expPath}/qc
#fi
#
#    if [ -d ${expPath}/bcFastq/fail/${bc} ];
#then
#    echo "failed reads ..."
#    cat ${expPath}/bcFastq/fail/${bc}/*.fastq > ${expPath}/bcFastq/${expName}_fail_${bc}.fastq
#    gzip ${expPath}/bcFastq/${expName}_fail_${bc}.fastq
#    ${NANOPOLISH_DIR}/nanopolish index -s ${summaryFile} -d ${expPath}/fast5Files ${expPath}/bcFastq/${expName}_fail_${bc}.fastq.gz
#    do_pauvre_nanoQC $bc fail $expName ${expPath}/qc
#fi
#
#
#################################################
## aligning to genome
#################################################${expPath}/meth_calls/${expName}_pass_${bc}_GpCcalls_${chr[$i]}.tsv
#
#echo "aligning to genome..."
#
#mkdir -p ${expPath}/bamFiles
#
## map reads to genome with minimap2
## filter reads with flag=2308: unmapped (4) + secondary alignment (256) + supplementary alignment (2048) [the latter category is the main problem]
#minimap2 -ax map-ont $genomeFile ${expPath}/bcFastq/${expName}_pass_${bc}.fastq.gz | samtools view -F 2308 -b | samtools sort -T tmp -o ${expPath}/bamFiles/${expName}_pass_${bc}.sorted.bam 
#minimap2 -ax map-ont $genomeFile ${expPath}/bcFastq/${expName}_fail_${bc}.fastq.gz | samtools view -F 2308 -b | samtools sort -T tmp -o ${expPath}/bamFiles/${expName}_fail_${bc}.sorted.bam
#
#echo "index bam file ..."
#samtools index ${expPath}/bamFiles/${expName}_pass_${bc}.sorted.bam
#samtools index ${expPath}/bamFiles/${expName}_fail_${bc}.sorted.bam


################################################
# getting chr intervals
################################################
if [ ! -f "chrom.sizes" ]
then
        faidx $genomeFile -i chromsizes > chrom.sizes
fi

chrIntervals=( `cut -f1,2 chrom.sizes | sed $'s/\t/:1-/g'` )
chr=( `cut -f1 chrom.sizes` )


################################################
# identifying CmG
################################################
echo "identify CmG ..."

mkdir -p ${expPath}/meth_calls/

for i in "${!chr[@]}"
do
${NANOPOLISH_DIR}/nanopolish call-methylation -t 4 -q cpg -w ${chrIntervals[$i]} -r ${expPath}/bcFastq/${expName}_pass_${bc}.fastq.gz -b ${expPath}/bamFiles/${expName}_pass_${bc}.sorted.bam -g $genomeFile > ${expPath}/meth_calls/${expName}_pass_${bc}_CpGcalls_${chr[$i]}.tsv

${NANOPOLISH_DIR}/nanopolish call-methylation -t 4 -q cpg -w ${chrIntervals[$i]} -r ${expPath}/bcFastq/${expName}_fail_${bc}.fastq.gz -b ${expPath}/bamFiles/${expName}_fail_${bc}.sorted.bam -g $genomeFile > ${expPath}/meth_calls/${expName}_fail_${bc}_CpGcalls_${chr[$i]}.tsv
done


#### combine separate chromosomes into single file ####

# first write header
head -1 ${expPath}/meth_calls/${expName}_pass_${bc}_CpGcalls_${chr[0]}.tsv > ${expPath}/meth_calls/${expName}_pass_${bc}_CpGcalls.tsv
head -1 ${expPath}/meth_calls/${expName}_fail_${bc}_CpGcalls_${chr[0]}.tsv > ${expPath}/meth_calls/${expName}_fail_${bc}_CpGcalls.tsv

# then combine files
for i in "${!chr[@]}"
do
        tail -n +2 ${expPath}/meth_calls/${expName}_pass_${bc}_CpGcalls_${chr[$i]}.tsv >> ${expPath}/meth_calls/${expName}_pass_${bc}_CpGcalls.tsv
        tail -n +2 ${expPath}/meth_calls/${expName}_fail_${bc}_CpGcalls_${chr[$i]}.tsv >> ${expPath}/meth_calls/${expName}_fail_${bc}_CpGcalls.tsv
        rm ${expPath}/meth_calls/${expName}_????_${bc}_CpGcalls_${chr[$i]}.tsv
done


#### caclulating frequency ######
mkdir -p ${expPath}/meth_freq
${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${expPath}/meth_calls/${expName}_pass_${bc}_CpGcalls.tsv > ${expPath}/meth_freq/${expName}_pass_${bc}_freqCmG.tsv

${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${expPath}/meth_calls/${expName}_fail_${bc}_CpGcalls.tsv > ${expPath}/meth_freq/${expName}_fail_${bc}_freqCmG.tsv


################################################
# identifying GCm
################################################
echo "identify GCm ..."

mkdir -p ${expPath}/meth_calls/

for i in "${!chr[@]}"
do
echo ${chr[$i]}
echo ${chrIntervals[$i]}
${NANOPOLISH_DIR}/nanopolish call-methylation -t 4 -q gpc -r ${expPath}/bcFastq/${expName}_pass_${bc}.fastq.gz -b ${expPath}/bamFiles/${expName}_pass_${bc}.sorted.bam -g $genomeFile -w ${chrIntervals[$i]} > ${expPath}/meth_calls/${expName}_pass_${bc}_GpCcalls_${chr[$i]}.tsv

${NANOPOLISH_DIR}/nanopolish call-methylation -t 4 -q gpc -r ${expPath}/bcFastq/${expName}_fail_${bc}.fastq.gz -b ${expPath}/bamFiles/${expName}_fail_${bc}.sorted.bam -g $genomeFile -w ${chrIntervals[$i]} > ${expPath}/meth_calls/${expName}_fail_${bc}_GpCcalls_${chr[$i]}.tsv
done


#### combine separate chromosomes into single file ####

# first write header
head -1 ${expPath}/meth_calls/${expName}_pass_${bc}_GpCcalls_${chr[0]}.tsv > ${expPath}/meth_calls/${expName}_pass_${bc}_GpCcalls.tsv
head -1 ${expPath}/meth_calls/${expName}_fail_${bc}_GpCcalls_${chr[0]}.tsv > ${expPath}/meth_calls/${expName}_fail_${bc}_GpCcalls.tsv 

# then combine files
for i in "${!chr[@]}"
do 
	tail -n +2 ${expPath}/meth_calls/${expName}_pass_${bc}_GpCcalls_${chr[$i]}.tsv >> ${expPath}/meth_calls/${expName}_pass_${bc}_GpCcalls.tsv
	tail -n +2 ${expPath}/meth_calls/${expName}_fail_${bc}_GpCcalls_${chr[$i]}.tsv >> ${expPath}/meth_calls/${expName}_fail_${bc}_GpCcalls.tsv
	rm ${expPath}/meth_calls/${expName}_????_${bc}_GpCcalls_${chr[$i]}.tsv
done


#### caclulating frequency ######
mkdir -p ${expPath}/meth_freq

${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${expPath}/meth_calls/${expName}_pass_${bc}_GpCcalls.tsv > ${expPath}/meth_freq/${expName}_pass_${bc}_freqGCm.tsv

${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${expPath}/meth_calls/${expName}_fail_${bc}_GpCcalls.tsv > ${expPath}/meth_freq/${expName}_fail_${bc}_freqGCm.tsv


