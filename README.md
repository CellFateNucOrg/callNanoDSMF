# nanoDSMF
Analysing nanopore sequencing of dSMF data

## Setup of environment

These scripts require:
1) albacore_env python environment, as set up for the nanoHiC analysis. This includes a local installation of MinIONQC.
2) nanopolish >= v0.11.0 must be installed locally (just follow instructions on github https://github.com/jts/nanopolish/blob/master/README.md)
3) The full path of the nanopolish directory must be exported from the .bashrc file as NANOPORE_DIR

## Setup of data

Since Minion often produced multiple runs per experiment (every time you stop and start the run), it is important to organise the data properly:

Create a directory with the name of your experiment (e.g 20171027_bTMP_meth). In that directory create a subdirectory called fast5files. Copy the directories of all the separate runs for this experiment (e.g. 20171027_1121_, 20171027_1128_) into the fast5files directory.

In the expriment directory clone the code from this repository:
```
git clone https://github.com/CellFateNucOrg/nanoDSMF.git
```
Your starting directory should now look something like this:

```
+---20171027_bTMP_meth/
|   +---fast5files/
|   |   +---20171027_1121_/
|   |   +---20171027_1128_/
|   +---nanoDSMF/
```

## Setup of variables

Go into the nanoDSMF data and copy the varSettings_example.sh to a new file without the word "example". Then edit this file to have the correct settings for your experiment.
```
cd nanoDSMF
cp varSettings_example.sh verSettings.sh
nano varSettings.sh
```

You need to specify: 
- an experiment name (expName). This will be used in naming you files. 
- the barcodes for the samples you analysed (write them with exactly the same format as in the example file).
- location of the reference genome file to which you want to align your sequences.
- location of the MinIONQC program (already set for the bioinformatics cluster)

## Basecalling and MinIONQC
The *basecallAlbacore.sh* script perfom basecalling with albacore and QC with MinIONQC. To run is use the wrapper script labelled 01:

```
sbatch 01_runBasecall.sh
```

The basecalled reads will all be under the directory (relative to the inside of nanoDSMF dir), in ../workspace, separated by barcode. The MinionQC will be in ../fastqQC/MinIONQC

## Aligning to genome and calling cytosine methylation
The *alignFastq.sh* script merges all files from a single barcode (separately for pass and fail) into a single file named according to the following template: experimentName_passORfail_barcode.fastq.gz. These files are placed in the ../fastqFiles directory.

Next the script aligns the reads to the genome and places the bam files (same naming scheme) in the ../bamFiles directory.

Finally, the CpG and GpC methlation is called with nanopolish. The resulting .tsv files will be placed in the ../methylation_calls directory. 

To run this script you must first modify the wrapper to specify the task array number - there will be one task for each barcode. So if in your varSettings.sh file you specified 4 barcodes, then you should edit the SBATCH --array line in the wrapper script *02_runAlignFastq.sh* as follows:

```
#SBATCH --array=1-4
```

Then run the wrapper script:
```
sbatch 02_runAlignFastq.sh
```





