#! /bin/bash
## script to arrange all fastq files into batches for analysis. Combined files will be created for
## each barcode of interest for both passed and failed reads (e.g. pass_barcode01 or fail_barcode01).  

# Get variables from command line
expName=$1 	# experiment name (date of exp usually)
bc=$2 		# barcode
genomeFile=$3 	# full path to reference genome
genomeName=`basename $genomeFile`
genomeName=${genomeName%.fasta}


####### modules to load ##########
module load vital-it
#module load R/3.5.1
module add UHTS/Analysis/minimap2/2.12;
module add UHTS/Analysis/samtools/1.8;
############################################




################################################
# aligning to genome
################################################

echo "aligning to genome..."

mkdir -p ../bamFilesdb/${genomeName}

# map reads to genome with minimap2
# filter reads with flag=2308: unmapped (4) + secondary alignment (256) + supplementary alignment (2048) [the latter category is the main problem]
minimap2 -ax map-ont $genomeFile ../dbbcFastq/pass/${bc}.fastq.gz | samtools view -F 2308 -b | samtools sort -T tmp -o ../bamFilesdb/${genomeName}/${expName}_pass_${bc}.sorted.bam 
minimap2 -ax map-ont $genomeFile ../dbbcFastq/fail/${bc}.fastq.gz | samtools view -F 2308 -b | samtools sort -T tmp -o ../bamFilesdb/${genomeName}/${expName}_fail_${bc}.sorted.bam

echo "index bam file ..."
samtools index ../bamFilesdb/${genomeName}/${expName}_pass_${bc}.sorted.bam
samtools index ../bamFilesdb/${genomeName}/${expName}_fail_${bc}.sorted.bam


################################################
# getting chr intervals
################################################
if [ ! -f "${genomeName}.chrom.sizes" ]
then
        faidx $genomeFile -i chromsizes > ${genomeName}.chrom.sizes
fi

chrIntervals=( `cut -f1,2 ${genomeName}.chrom.sizes | sed $'s/\t/:1-/g'` )
chr=( `cut -f1 ${genomeName}.chrom.sizes` )



################################################
# identifying CmG
################################################
echo "identify CmG ..."

mkdir -p ../meth_callsdb/${genomeName}
for i in "${!chr[@]}"
do
${NANOPOLISH_DIR}/nanopolish call-methylation -t 4 -q cpg -w ${chrIntervals[$i]} -r ${workDir}/dbbcFastq/pass/${bc}.fastq.gz -b ${workDir}/bamFilesdb/${genomeName}/${expName}_pass_${bc}.sorted.bam -g $genomeFile > ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_CpGcalls_${chr[$i]}.tsv

${NANOPOLISH_DIR}/nanopolish call-methylation -t 4 -q cpg -w ${chrIntervals[$i]} -r ${workDir}/dbbcFastq/fail/${bc}.fastq.gz -b ${workDir}/bamFilesdb/${genomeName}/${expName}_fail_${bc}.sorted.bam -g $genomeFile > ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_CpGcalls_${chr[$i]}.tsv
done


#### combine separate chromosomes into single file ####

# first write header
head -1 ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_CpGcalls_${chr[0]}.tsv > ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_CpGcalls.tsv
head -1 ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_CpGcalls_${chr[0]}.tsv > ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_CpGcalls.tsv

# then combine files
for i in "${!chr[@]}"
do
        tail -n +2 ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_CpGcalls_${chr[$i]}.tsv >> ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_CpGcalls.tsv
        tail -n +2 ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_CpGcalls_${chr[$i]}.tsv >> ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_CpGcalls.tsv
        rm ${workDir}/meth_callsdb/${genomeName}/${expName}_????_${bc}_CpGcalls_${chr[$i]}.tsv
done


#### caclulating frequency ######
mkdir -p ../meth_freq/${genomeName}
${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_CpGcalls.tsv > ${workDir}/meth_freqdb/${genomeName}/${expName}_pass_${bc}_freqCmG.tsv

${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_CpGcalls.tsv > ${workDir}/meth_freqdb/${genomeName}/${expName}_fail_${bc}_freqCmG.tsv



################################################
# identifying GCm
################################################
echo "identify GCm ..."

mkdir -p ../meth_callsdb/${genomeName}
for i in "${!chr[@]}"
do
${NANOPOLISH_DIR}/nanopolish call-methylation -t 4 -q gpc -w ${chrIntervals[$i]} -r ${workDir}/dbbcFastq/pass/${bc}.fastq.gz -b ${workDir}/bamFilesdb/${genomeName}/${expName}_pass_${bc}.sorted.bam -g $genomeFile > ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_GpCcalls_${chr[$i]}.tsv

${NANOPOLISH_DIR}/nanopolish call-methylation -t 4 -q gpc -w ${chrIntervals[$i]} -r ${workDir}/dbbcFastq/fail/${bc}.fastq.gz -b ${workDir}/bamFilesdb/${genomeName}/${expName}_fail_${bc}.sorted.bam -g $genomeFile > ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_GpCcalls_${chr[$i]}.tsv
done


#### combine separate chromosomes into single file ####

# first write header
head -1 ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_GpCcalls_${chr[0]}.tsv > ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_GpCcalls.tsv
head -1 ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_GpCcalls_${chr[0]}.tsv > ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_GpCcalls.tsv 

# then combine files
for i in "${!chr[@]}"
do 
    tail -n +2 ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_GpCcalls_${chr[$i]}.tsv >> ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_GpCcalls.tsv
    tail -n +2 ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_GpCcalls_${chr[$i]}.tsv >> ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_GpCcalls.tsv
    rm ${workDir}/meth_callsdb/${genomeName}/${expName}_????_${bc}_GpCcalls_${chr[$i]}.tsv
done



#### caclulating frequency ######
mkdir -p ../meth_freqdb/${genomeName}
${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_callsdb/${genomeName}/${expName}_pass_${bc}_GpCcalls.tsv > ${workDir}/meth_freqdb/${genomeName}/${expName}_pass_${bc}_freqGCm.tsv

${NANOPOLISH_DIR}/scripts/calculate_methylation_frequency.py -i ${workDir}/meth_callsdb/${genomeName}/${expName}_fail_${bc}_GpCcalls.tsv > ${workDir}/meth_freqdb/${genomeName}/${expName}_fail_${bc}_freqGCm.tsv

