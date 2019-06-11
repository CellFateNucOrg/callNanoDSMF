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


################################################
# Collecting reads from barcodes that were used
################################################
echo "collecting reads from folder of barcodes that were used..."

# merge all reads from particular barcode into single file (pass fail separately)
echo ${dataDir}/fast5Files
mkdir -p ${workDir}/qc/NanoStatdb

if [ -f ${workDir}/dbbcFastq/pass/${bc}.fastq.gz ];
then
    echo "passed reads..."
    nanopolish index -s ${summaryFile} -d ${dataDir}/fast5Files ${workDir}/dbbcFastq/pass/${bc}.fastq.gz
    mkdir -p ${workDir}/qc/db_pass_${bc}
    nanoQC ${workDir}/dbbcFastq/pass/${bc}.fastq.gz -o ${workDir}/qc/db_pass_${bc}
    NanoStat --fastq ${workDir}/dbbcFastq/pass/${bc}.fastq.gz  --outdir ${workDir}/qc/NanoStatdb --name NanoStat_${expName}_pass_${bc}.txt --readtype 1D 
fi

if [ -f ${workDir}/dbbcFastq/fail/${bc}.fastq.gz ];
then
    echo "failed reads ..."
    nanopolish index -s ${summaryFile} -d ${dataDir}/fast5Files ${workDir}/dbbcFastq/fail/${bc}.fastq.gz
    mkdir -p ${workDir}/qc/db_fail_${bc}
    nanoQC ${workDir}/dbbcFastq/fail/${bc}.fastq.gz -o ${workDir}/qc/db_fail_${bc}
    NanoStat --fastq ${workDir}/dbbcFastq/fail/${bc}.fastq.gz  --outdir ${workDir}/qc/NanoStatdb --name NanoStat_${expName}_fail_${bc}.txt --readtype 1D
fi


################################################
# aligning to genome
################################################

echo "aligning to genome..."

mkdir -p ${workDir}/bamFilesdb

# map reads to genome with minimap2
# filter reads with flag=2308: unmapped (4) + secondary alignment (256) + supplementary alignment (2048) [the latter category is the main problem]
minimap2 -ax map-ont $genomeFile ${workDir}/dbbcFastq/pass/${bc}.fastq.gz | samtools view -F 2308 -b | samtools sort -T pass_${bc} -o ${workDir}/bamFilesdb/${expName}_pass_${bc}.sorted.bam 
minimap2 -ax map-ont $genomeFile ${workDir}/dbbcFastq/fail/${bc}.fastq.gz | samtools view -F 2308 -b | samtools sort -T fail_${bc} -o ${workDir}/bamFilesdb/${expName}_fail_${bc}.sorted.bam

echo "index bam file ..."
samtools index ${workDir}/bamFilesdb/${expName}_pass_${bc}.sorted.bam
samtools index ${workDir}/bamFilesdb/${expName}_fail_${bc}.sorted.bam


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

mkdir -p ${workDir}/meth_callsdb/

for i in "${!chr[@]}"
do
	nanopolish call-methylation -t 4 -q cpg -w ${chrIntervals[$i]} -r ${workDir}/dbbcFastq/pass/${bc}.fastq.gz -b ${workDir}/bamFilesdb/${expName}_pass_${bc}.sorted.bam -g $genomeFile > ${workDir}/meth_callsdb/${expName}_pass_${bc}_CpGcalls_${chr[$i]}.tsv

	nanopolish call-methylation -t 4 -q cpg -w ${chrIntervals[$i]} -r ${workDir}/dbbcFastq/fail/${bc}.fastq.gz -b ${workDir}/bamFilesdb/${expName}_fail_${bc}.sorted.bam -g $genomeFile > ${workDir}/meth_callsdb/${expName}_fail_${bc}_CpGcalls_${chr[$i]}.tsv
done


#### combine separate chromosomes into single file ####

# first write header
head -1 ${workDir}/meth_callsdb/${expName}_pass_${bc}_CpGcalls_${chr[0]}.tsv > ${workDir}/meth_callsdb/${expName}_pass_${bc}_CpGcalls.tsv
head -1 ${workDir}/meth_callsdb/${expName}_fail_${bc}_CpGcalls_${chr[0]}.tsv > ${workDir}/meth_callsdb/${expName}_fail_${bc}_CpGcalls.tsv

# then combine files
for i in "${!chr[@]}"
do
        tail -n +2 ${workDir}/meth_callsdb/${expName}_pass_${bc}_CpGcalls_${chr[$i]}.tsv >> ${workDir}/meth_callsdb/${expName}_pass_${bc}_CpGcalls.tsv
        tail -n +2 ${workDir}/meth_callsdb/${expName}_fail_${bc}_CpGcalls_${chr[$i]}.tsv >> ${workDir}/meth_callsdb/${expName}_fail_${bc}_CpGcalls.tsv
        rm ${workDir}/meth_callsdb/${expName}_????_${bc}_CpGcalls_${chr[$i]}.tsv
done


#### caclulating frequency ######
mkdir -p ${workDir}/meth_freqdb
${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_callsdb/${expName}_pass_${bc}_CpGcalls.tsv > ${workDir}/meth_freqdb/${expName}_pass_${bc}_freqCmG.tsv

${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_callsdb/${expName}_fail_${bc}_CpGcalls.tsv > ${workDir}/meth_freqdb/${expName}_fail_${bc}_freqCmG.tsv


################################################
# identifying GCm
################################################
echo "identify GCm ..."

mkdir -p ${workDir}/meth_callsdb/

for i in "${!chr[@]}"
do
echo ${chr[$i]}
echo ${chrIntervals[$i]}
	nanopolish call-methylation -t 4 -q gpc -r ${workDir}/dbbcFastq/pass/${bc}.fastq.gz -b ${workDir}/bamFilesdb/${expName}_pass_${bc}.sorted.bam -g $genomeFile -w ${chrIntervals[$i]} > ${workDir}/meth_callsdb/${expName}_pass_${bc}_GpCcalls_${chr[$i]}.tsv

	nanopolish call-methylation -t 4 -q gpc -r ${workDir}/dbbcFastq/fail/${bc}.fastq.gz -b ${workDir}/bamFilesdb/${expName}_fail_${bc}.sorted.bam -g $genomeFile -w ${chrIntervals[$i]} > ${workDir}/meth_callsdb/${expName}_fail_${bc}_GpCcalls_${chr[$i]}.tsv
done


#### combine separate chromosomes into single file ####

# first write header
head -1 ${workDir}/meth_callsdb/${expName}_pass_${bc}_GpCcalls_${chr[0]}.tsv > ${workDir}/meth_callsdb/${expName}_pass_${bc}_GpCcalls.tsv
head -1 ${workDir}/meth_callsdb/${expName}_fail_${bc}_GpCcalls_${chr[0]}.tsv > ${workDir}/meth_callsdb/${expName}_fail_${bc}_GpCcalls.tsv 

# then combine files
for i in "${!chr[@]}"
do 
	tail -n +2 ${workDir}/meth_callsdb/${expName}_pass_${bc}_GpCcalls_${chr[$i]}.tsv >> ${workDir}/meth_callsdb/${expName}_pass_${bc}_GpCcalls.tsv
	tail -n +2 ${workDir}/meth_callsdb/${expName}_fail_${bc}_GpCcalls_${chr[$i]}.tsv >> ${workDir}/meth_callsdb/${expName}_fail_${bc}_GpCcalls.tsv
	rm ${workDir}/meth_callsdb/${expName}_????_${bc}_GpCcalls_${chr[$i]}.tsv
done


#### caclulating frequency ######
mkdir -p ${workDir}/meth_freqdb

${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_callsdb/${expName}_pass_${bc}_GpCcalls.tsv > ${workDir}/meth_freqdb/${expName}_pass_${bc}_freqGCm.tsv

${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_callsdb/${expName}_fail_${bc}_GpCcalls.tsv > ${workDir}/meth_freqdb/${expName}_fail_${bc}_freqGCm.tsv


