#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --output=SomaticAmplicon-%N-%j.output
#SBATCH --error=SomaticAmplicon-%N-%j.error
#SBATCH --partition=high

# this is the latest version for WREN
# Description: Somatic Amplcon Pipeline (Illumina paired-end). Not for use with other library preps.
# Author: AWMGS
# Mode: BY_SAMPLE
# Use: sbatch within sample directory

# version=2.0.3

set -euo pipefail


version="master"

# Directory structure required for pipeline
#
# /data
# └──
#	output
#     └──
# 		results
#         └── seqId
#         	├── panel1
#         	│   ├── sample1
#         	│   ├── sample2
#         	│   └── sample3
#         	└── panel2
#             	├── sample1
#             	├── sample2
#             	└── sample3
#
# Script 1 runs in sample folder, requires fastq files split by lane

######################################################################
#                            MODULES                                 #
######################################################################
module load singularity

######################################################################
#                      FUNCTIONS/VARIABLES			     #
######################################################################

countQCFlagFails() {
    #count how many core FASTQC tests failed
    # grep -E is an extended regular expression
    # grep -v is an invert match, and selects non-matching lines
    # this is basically scanning the summary.txt file from fastqc output and pulling out specfic columns and then omitting those that have PASS or WARN and leaving only FAIL if any
    # it then counts how many lines have failed
    grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    grep -v ^WARN | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}

# Define location of Singulaity SIF files
SIF="/data/resources/envs/sifs/conda.sif" # Path to sif
SIFGATK="/data/resources/envs/sifs/gatk3_3.7-0.sif" # Path to sif
SIFPISCES="/data/resources/envs/sifs/pisces.sif" # Path to sif
SIFCOVER="/data/resources/envs/sifs/coverage.sif" # Path to sif
SIFBED="/data/resources/envs/sifs/bed2hgvs.sif" # Path to sif
SIFVHOOD="/data/resources/envs/sifs/virtualhood.sif" # Path to sif

# Define Executables
PICARD="singularity exec --bind /Output,/localscratch,/data:/data $SIF picard -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/localscratch -Xmx32g" # making code look cleaner
SINGULARITY="singularity exec --bind /Output,/localscratch,/data:/data $SIF" # Initiating the singularity exec --bind /data:/data command
GATK="singularity exec --bind /Output,/localscratch,/data:/data $SIFGATK java -Djava.io.tmpdir=/localscratch -Xmx32g -jar /usr/GenomeAnalysisTK.jar -T" # Initiating the GATK singularity container and command
PISCES="singularity exec --bind /Output,/localscratch,/data:/data $SIFPISCES dotnet /app/Pisces_5.2.9.122/Pisces.dll"
AMPLICON="singularity exec --bind /Output,/localscratch,/data:/data $SIF java -jar /opt/conda/bin/AmpliconRealigner-1.1.1.jar"
SOFTCLIP="singularity exec --bind /Output,/localscratch,/data:/data $SIF java -Xmx2g -jar /opt/conda/bin/SoftClipPCRPrimer-1.1.0.jar"
COVERAGE="singularity exec --bind /Output,/localscratch,/data:/data $SIF java -Djava.io.tmpdir=/localscratch -Xmx8g -jar /opt/conda/bin/CoverageCalculator-2.0.2.jar"
VCFPARSE="singularity exec --bind /Output,/localscratch,/data:/data $SIF python /opt/conda/bin/vcf_parse-0.1.2/vcf_parse.py"
COVERCALC="singularity exec --bind /Output,/localscratch,/data:/data $SIFCOVER python /opt/conda/bin/CoverageCalculatorPy/CoverageCalculatorPy.py"
BED="singularity exec --bind /Output,/localscratch,/data:/data $SIFBED Rscript /opt/conda/bin/bed2hgvs-v0.3.0/bed2hgvs.R"
VHOOD="singularity exec --bind /Output,/localscratch,/data:/data $SIFVHOOD python"

######################################################################
#			PIPELINE				     #
######################################################################

### load sample & pipeline variables ###
. *.variables
. /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel".variables

### Preprocessing ###


### Run level steps ###
## This block should only be carried out when all samples for the panel have been processed

# Creating a marker file to then decide whether block below should be executed or not
touch move_complete.txt

# number of samples to be processed (i.e. count variables files)/ number of samples that have completed
expected=$(for i in /data/output/results/"$seqId"/"$panel"/*/*.variables; do echo $i; done | wc -l)
complete=$(for i in /data/output/results/"$seqId"/"$panel"/*/move_complete.txt; do echo $i; done | wc -l)

if [ $complete -eq $expected ]; then

    # Merge QC files
    python /data/diagnostics/scripts/merge_qc_files.py ..

    # BRCA merge report files
    if [ $merge_reports == true ]; then

        # get report headers
#        cat $(ls /data/output/results/"$seqId"/"$panel"/*/*VariantReport.txt | head -n1) | head -n1 > /data/output/results/"$seqId"/"$panel"/"$seqId"_merged_variant_report.txt
#        echo -e "Sample\tBRCA1_500X\tBRCA2_500X\tBRCA1_100X\tBRCA2_100X" > /data/output/results/"$seqId"/"$panel"/"$seqId"_merged_coverage_report.txt

        # loop over all samples and merge reports
        for sample_path in /data/output/results/"$seqId"/"$panel"/*/; do
            sample=$(basename $sample_path)
            echo "Merging coverage and variant reports for $sample"

            # merge variant report
            cat "$sample_path"/*VariantReport.txt | tail -n+2 >> /data/output/results/"$seqId"/"$panel"/"$seqId"_merged_variant_report.txt

            # rename percentagecoverage to percebtage coverage 500x and 500x gaps file
            mv "$sample_path"/"$seqId"_"$sample"_PercentageCoverage.txt "$sample_path"/"$seqId"_"$sample"_PercentageCoverage_500x.txt
            mv "$sample_path"/"$sample"_gaps.bed "$sample_path"/"$sample"_gaps_500x.bed

            # Calculate gene (clinical) percentage coverage at 100x
            $COVERAGE \
            $sample_path/"$seqId"_"$sample"_DepthOfCoverage \
            /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_genes.txt \
            /data/resources/human/refseq/ref_GRCh37.p13_top_level.gff3 \
            -p5 \
            -d100 \
            > "$sample_path"/"$seqId"_"$sample"_PercentageCoverage_100x.txt

            # rename 100x gaps file and move into sample folder
            mv "$sample"_gaps.bed "$sample_path"/"$sample"_gaps_100x.bed

            # merge 500x and 100x coverage reports into one file
            brca1_500x=$(grep BRCA1 $sample_path/"$seqId"_"$sample"_PercentageCoverage_500x.txt | cut -f3)
            brca2_500x=$(grep BRCA2 $sample_path/"$seqId"_"$sample"_PercentageCoverage_500x.txt | cut -f3)
            brca1_100x=$(grep BRCA1 $sample_path/"$seqId"_"$sample"_PercentageCoverage_100x.txt | cut -f3)
            brca2_100x=$(grep BRCA2 $sample_path/"$seqId"_"$sample"_PercentageCoverage_100x.txt | cut -f3)
            echo -e "$sample\t$brca1_500x\t$brca2_500x\t$brca1_100x\t$brca2_100x" >> /data/output/results/"$seqId"/"$panel"/"$seqId"_merged_coverage_report.txt

            # reset variables
            unset sample brca1_500x brca2_500x brca1_100x brca2_100x
        done
    fi

    # virtual hood
    if [ $generate_worksheets == true ]; then
    
        # identify name of NTC
        ntc=$(for s in /data/output/results/"$seqId"/"$panel"/*/; do echo $(basename $s);done | grep 'NTC')

        # loop over all samples and generate a report
        for sample_path in /data/output/results/"$seqId"/"$panel"/*/; do
            
            # clear previous instance
            unset referral 
            
            # set variables
            sample=$(basename $sample_path)
            # Change this path so not hardcoded 
            . /data/output/results/"$seqId"/"$panel"/"$sample"/*.variables
            echo "Generating worksheet for $sample"

            # check that referral variable is defined, if not set as NA
            if [ -z $referral ]; then referral=NA; fi

            # do not generate report where NTC is the query sample
            if [ $sample != $ntc ]; then

                if [ $referral == 'melanoma' ] || [ $referral == 'lung' ] || [ $referral == 'colorectal' ] || [ $referral == 'glioma' ] || [ $referral == 'tumour' ] || [ $referral == 'gist' ] || [ $referral == 'thyroid' ]; then
                    $VHOOD /opt/conda/bin/VirtualHood-1.2.0/CRM_report_new_referrals.py --runid $seqId --sampleid $sample --worksheet $worklistId --referral $referral --NTC_name $ntc --path /data/output/results/"$seqId"/"$panel"/ --artefacts /data/temp/artefacts_lists/
                fi
            fi
        done

    fi

fi

#load sample & pipeline variables
. *.variables
. /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel".variables

### Clean up ###

#delete unused files
rm "$seqId"_"$sampleId"_*unaligned.bam "$seqId"_"$sampleId"_aligned.bam "$seqId"_"$sampleId"_aligned.bai "$seqId"_"$sampleId"_amplicon_realigned.bam
rm "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam.bai "$seqId"_"$sampleId"_indel_realigned.intervals
rm "$seqId"_"$sampleId"_clipped.bam "$seqId"_"$sampleId"_clipped_sorted.bam "$seqId"_"$sampleId"_clipped_sorted.bam.bai "$panel"_ROI.interval_list "$panel"_ROI_b37_thick.bed
rm "$seqId"_"$sampleId"_left_aligned.vcf "$seqId"_"$sampleId"_left_aligned.vcf.idx "$seqId"_"$sampleId".bam.bai "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam
rm "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bai "$seqId"_"$sampleId"_filtered_meta.vcf "$seqId"_"$sampleId"_filtered_meta.vcf.idx "$seqId"_"$sampleId"_filtered.vcf
rm "$seqId"_"$sampleId"_filtered.vcf.idx "$seqId"_"$sampleId"_fixed.vcf "$seqId"_"$sampleId"_fixed.vcf.idx "$seqId"_"$sampleId"_indel_realigned.bam "$seqId"_"$sampleId"_indel_realigned.bai
rm "$seqId"_"$sampleId"_*_fastqc.zip "$seqId"_"$sampleId"_lcr.vcf "$seqId"_"$sampleId"_lcr.vcf.idx "$seqId"_"$sampleId"_left_aligned_annotated.vcf "$seqId"_"$sampleId"_left_aligned_annotated.vcf.idx
rm "$seqId"_"$sampleId".vcf

# create complete marker
touch 1_SomaticAmplicon.sh.e69420
