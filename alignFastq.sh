#! /bin/bash
## script to arrange all fastq files into batches for analysis. Combined files will be created for
## each barcode of interest for both passed and failed reads (e.g. pass_barcode01 or fail_barcode01).  

source ./varSettings.sh
# Get variables from command line
expName=$1 	# experiment name (date of exp usually)
bc=$2 		# barcode
genomeFile=$3 	# full path to reference genome

# need absolute paths for nanopolish index. get it from the summary file.
summaryFile=${workDir}/fastqFiles/sequencing_summary.txt

####### modules to load ##########
module load vital-it
#module load R/3.5.1
module add UHTS/Analysis/minimap2/2.12;
module add UHTS/Analysis/samtools/1.8;
############################################

source ${HOME}/.bashrc
source ${CONDA_ACTIVATE} nanopore


#################################################
## Collecting reads from barcodes that were used
#################################################
#echo "collecting reads from folder of barcodes that were used..."
#
## merge all reads from particular barcode into single file (pass fail separately)
#echo ${dataDir}/fast5Files
#mkdir -p ${workDir}/qc/NanoStat
#
#if [ -d "${workDir}/bcFastq/pass/${bc}" ];
#then
#    echo "passed reads..."
#    cat ${workDir}/bcFastq/pass/${bc}/*.fastq > ${workDir}/bcFastq/${expName}_pass_${bc}.fastq
#    gzip ${workDir}/bcFastq/${expName}_pass_${bc}.fastq
#    rm ${workDir}/bcFastq/${expName}_pass_${bc}.fastq
#    nanopolish index -s ${summaryFile} -d ${dataDir}/fast5Files ${workDir}/bcFastq/${expName}_pass_${bc}.fastq.gz
#    mkdir -p ${workDir}/qc/pass_${bc}
#    nanoQC ${workDir}/bcFastq/${expName}_pass_${bc}.fastq.gz -o ${workDir}/qc/pass_${bc}
#    NanoStat --fastq ${workDir}/bcFastq/${expName}_pass_${bc}.fastq.gz  --outdir ${workDir}/qc/NanoStat --name NanoStat_${expName}_pass_${bc}.txt --readtype 1D 
#fi
#
#if [ -d "${workDir}/bcFastq/fail/${bc}" ];
#then
#    echo "failed reads ..."
#    cat ${workDir}/bcFastq/fail/${bc}/*.fastq > ${workDir}/bcFastq/${expName}_fail_${bc}.fastq
#    gzip ${workDir}/bcFastq/${expName}_fail_${bc}.fastq
#    rm ${workDir}/bcFastq/${expName}_fail_${bc}.fastq
#    nanopolish index -s ${summaryFile} -d ${dataDir}/fast5Files ${workDir}/bcFastq/${expName}_fail_${bc}.fastq.gz
#    mkdir -p ${workDir}/qc/fail_${bc}
#    nanoQC ${workDir}/bcFastq/${expName}_fail_${bc}.fastq.gz -o ${workDir}/qc/fail_${bc}
#    NanoStat --fastq ${workDir}/bcFastq/${expName}_fail_${bc}.fastq.gz  --outdir ${workDir}/qc/NanoStat --name NanoStat_${expName}_fail_${bc}.txt --readtype 1D
#fi
#
#
#################################################
## aligning to genome
#################################################
#
#echo "aligning to genome..."
#
#mkdir -p ${workDir}/bamFiles
#
## map reads to genome with minimap2
## filter reads with flag=2308: unmapped (4) + secondary alignment (256) + supplementary alignment (2048) [the latter category is the main problem]
#minimap2 -ax map-ont $genomeFile ${workDir}/bcFastq/${expName}_pass_${bc}.fastq.gz | samtools view -F 2308 -b | samtools sort -T pass_${bc}  -o ${workDir}/bamFiles/${expName}_pass_${bc}.sorted.bam 
#minimap2 -ax map-ont $genomeFile ${workDir}/bcFastq/${expName}_fail_${bc}.fastq.gz | samtools view -F 2308 -b | samtools sort -T fail_${bc} -o ${workDir}/bamFiles/${expName}_fail_${bc}.sorted.bam
#
#echo "index bam file ..."
#samtools index ${workDir}/bamFiles/${expName}_pass_${bc}.sorted.bam
#samtools index ${workDir}/bamFiles/${expName}_fail_${bc}.sorted.bam
#

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
doFail=FALSE

mkdir -p ${workDir}/meth_calls/

for i in "${!chr[@]}"
do
	nanopolish call-methylation -t 4 -q cpg -w ${chrIntervals[$i]} -r ${workDir}/bcFastq/${expName}_pass_${bc}.fastq.gz -b ${workDir}/bamFiles/${expName}_pass_${bc}.sorted.bam -g $genomeFile > ${workDir}/meth_calls/${expName}_pass_${bc}_CpGcalls_${chr[$i]}.tsv
  if [ "$doFail" = "TRUE" ]
  then
	echo "analysing failed reads"
	nanopolish call-methylation -t 4 -q cpg -w ${chrIntervals[$i]} -r ${workDir}/bcFastq/${expName}_fail_${bc}.fastq.gz -b ${workDir}/bamFiles/${expName}_fail_${bc}.sorted.bam -g $genomeFile > ${workDir}/meth_calls/${expName}_fail_${bc}_CpGcalls_${chr[$i]}.tsv
  fi
done


#### combine separate chromosomes into single file ####

# first write header
head -1 ${workDir}/meth_calls/${expName}_pass_${bc}_CpGcalls_${chr[0]}.tsv > ${workDir}/meth_calls/${expName}_pass_${bc}_CpGcalls.tsv
if [ "$doFail" = TRUE ]
then
  head -1 ${workDir}/meth_calls/${expName}_fail_${bc}_CpGcalls_${chr[0]}.tsv > ${workDir}/meth_calls/${expName}_fail_${bc}_CpGcalls.tsv
fi


# then combine files
for i in "${!chr[@]}"
do
        tail -n +2 ${workDir}/meth_calls/${expName}_pass_${bc}_CpGcalls_${chr[$i]}.tsv >> ${workDir}/meth_calls/${expName}_pass_${bc}_CpGcalls.tsv
  if [ "$doFail" = TRUE ]
  then
        tail -n +2 ${workDir}/meth_calls/${expName}_fail_${bc}_CpGcalls_${chr[$i]}.tsv >> ${workDir}/meth_calls/${expName}_fail_${bc}_CpGcalls.tsv
  fi
        rm ${workDir}/meth_calls/${expName}_????_${bc}_CpGcalls_${chr[$i]}.tsv
done


#### caclulating frequency ######
mkdir -p ${workDir}/meth_freq
${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_calls/${expName}_pass_${bc}_CpGcalls.tsv > ${workDir}/meth_freq/${expName}_pass_${bc}_freqCmG.tsv

if [ "$doFail" = TRUE ]
then
	${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_calls/${expName}_fail_${bc}_CpGcalls.tsv > ${workDir}/meth_freq/${expName}_fail_${bc}_freqCmG.tsv
fi

################################################
# identifying GCm
################################################
echo "identify GCm ..."

mkdir -p ${workDir}/meth_calls/

for i in "${!chr[@]}"
do
	echo ${chr[$i]}
	echo ${chrIntervals[$i]}
	nanopolish call-methylation -t 4 -q gpc -r ${workDir}/bcFastq/${expName}_pass_${bc}.fastq.gz -b ${workDir}/bamFiles/${expName}_pass_${bc}.sorted.bam -g $genomeFile -w ${chrIntervals[$i]} > ${workDir}/meth_calls/${expName}_pass_${bc}_GpCcalls_${chr[$i]}.tsv

  	if [ "${doFail}" = TRUE ]
  	then
		nanopolish call-methylation -t 4 -q gpc -r ${workDir}/bcFastq/${expName}_fail_${bc}.fastq.gz -b ${workDir}/bamFiles/${expName}_fail_${bc}.sorted.bam -g $genomeFile -w ${chrIntervals[$i]} > ${workDir}/meth_calls/${expName}_fail_${bc}_GpCcalls_${chr[$i]}.tsv
	fi
done


#### combine separate chromosomes into single file ####

# first write header
head -1 ${workDir}/meth_calls/${expName}_pass_${bc}_GpCcalls_${chr[0]}.tsv > ${workDir}/meth_calls/${expName}_pass_${bc}_GpCcalls.tsv
if [ "$doFail" = TRUE ]
then
	head -1 ${workDir}/meth_calls/${expName}_fail_${bc}_GpCcalls_${chr[0]}.tsv > ${workDir}/meth_calls/${expName}_fail_${bc}_GpCcalls.tsv 
fi

# then combine files
for i in "${!chr[@]}"
do 
	tail -n +2 ${workDir}/meth_calls/${expName}_pass_${bc}_GpCcalls_${chr[$i]}.tsv >> ${workDir}/meth_calls/${expName}_pass_${bc}_GpCcalls.tsv
  	if [ "$doFail" = TRUE ]
  	then
		tail -n +2 ${workDir}/meth_calls/${expName}_fail_${bc}_GpCcalls_${chr[$i]}.tsv >> ${workDir}/meth_calls/${expName}_fail_${bc}_GpCcalls.tsv
	fi
	rm ${workDir}/meth_calls/${expName}_????_${bc}_GpCcalls_${chr[$i]}.tsv
done


#### caclulating frequency ######
mkdir -p ${workDir}/meth_freq

${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_calls/${expName}_pass_${bc}_GpCcalls.tsv > ${workDir}/meth_freq/${expName}_pass_${bc}_freqGCm.tsv

if [ "$doFail" = TRUE ]
then
	${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_calls/${expName}_fail_${bc}_GpCcalls.tsv > ${workDir}/meth_freq/${expName}_fail_${bc}_freqGCm.tsv
fi


