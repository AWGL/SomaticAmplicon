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

#record FASTQC pass/fail
rawSequenceQuality=PASS

#convert FASTQ to uBAM & add RGIDs
for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    # generating the lane ID e.g. L001
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    # generating the name of the R1 file
    read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
    # generating the name of the R2 file
    read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

    $SINGULARITY cutadapt \
    -a "$read1Adapter" \
    -A "$read2Adapter" \
    -m 50 \
    -o "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -p "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    "$read1Fastq" \
    "$read2Fastq"   

    #merge overlapping reads
    $SINGULARITY pear \
    -f "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -r "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    -o "$seqId"_"$sampleId"_"$laneId"_merged.fastq \
    -j 10

    #convert fastq to ubam
    $PICARD FastqToSam \
    F1="$seqId"_"$sampleId"_"$laneId"_merged.fastq.assembled.fastq \
    O="$seqId"_"$sampleId"_"$laneId"_unaligned.bam \
    QUALITY_FORMAT=Standard \
    READ_GROUP_NAME="$seqId"_"$laneId"_"$sampleId" \
    SAMPLE_NAME="$sampleId" \
    LIBRARY_NAME="$worklistId"_"$sampleId"_"$panel" \
    PLATFORM_UNIT="$seqId"_"$laneId" \
    PLATFORM="ILLUMINA" \
    SEQUENCING_CENTER="IMG" \
    PREDICTED_INSERT_SIZE="$expectedInsertSize" \
    SORT_ORDER=queryname \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR=/localscratch 

    #fastqc
    $SINGULARITY fastqc -d /localscratch --threads 20 --extract "$seqId"_"$sampleId"_"$laneId"_R1.fastq
    $SINGULARITY fastqc -d /localscratch --threads 20 --extract "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    #check FASTQC output
    if [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R1_fastqc/summary.txt) -gt 0 ] || [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R2_fastqc/summary.txt) -gt 0 ]; then
        rawSequenceQuality=FAIL
    fi

    #clean up
    rm "$seqId"_"$sampleId"_"$laneId"_R1.fastq "$seqId"_"$sampleId"_"$laneId"_R2.fastq "$seqId"_"$sampleId"_"$laneId"_merged.fastq.*

done

#merge lane bams
$PICARD MergeSamFiles \
$(ls "$seqId"_"$sampleId"_*_unaligned.bam | sed 's/^/I=/' | tr '\n' ' ') \
SORT_ORDER=queryname \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT \
USE_THREADING=true \
MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=/localscratch \
O="$seqId"_"$sampleId"_unaligned.bam

# uBam2fq, map & MergeBamAlignment -taking the merged unaligned bam then piping that into bwa using the unaligned bam as one of the inputs
$PICARD SamToFastq \
I="$seqId"_"$sampleId"_unaligned.bam \
FASTQ=/dev/stdout \
NON_PF=true \
MAX_RECORDS_IN_RAM=2000000 \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=/localscratch | \
$SINGULARITY bwa mem \
-M \
-t 20 \
-p \
/data/resources/human/mappers/b37/bwa/human_g1k_v37.fasta \
/dev/stdin | \
$PICARD MergeBamAlignment \
ATTRIBUTES_TO_RETAIN=X0 \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM="$seqId"_"$sampleId"_unaligned.bam \
OUTPUT="$seqId"_"$sampleId"_aligned.bam \
R=/data/resources/human/mappers/b37/bwa/human_g1k_v37.fasta \
PAIRED_RUN=false \
SORT_ORDER="coordinate" \
IS_BISULFITE_SEQUENCE=false \
ALIGNED_READS_ONLY=false \
CLIP_ADAPTERS=false \
MAX_RECORDS_IN_RAM=2000000 \
MAX_INSERTIONS_OR_DELETIONS=-1 \
UNMAP_CONTAMINANT_READS=false \
CLIP_OVERLAPPING_READS=false \
ALIGNER_PROPER_PAIR_FLAGS=false \
ATTRIBUTES_TO_RETAIN=XS \
INCLUDE_SECONDARY_ALIGNMENTS=true \
CREATE_INDEX=true \
TMP_DIR=/localscratch


#Realign soft clipped bases
$AMPLICON \
-I "$seqId"_"$sampleId"_aligned.bam \
-O "$seqId"_"$sampleId"_amplicon_realigned.bam \
-R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed

#sort and index BAM
$SINGULARITY samtools sort -@5 -m16G -o "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam "$seqId"_"$sampleId"_amplicon_realigned.bam
$SINGULARITY samtools index "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam

#left align indels
$GATK LeftAlignIndels \
-R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam \
-o "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam \
-dt NONE

#Identify regions requiring realignment
$GATK RealignerTargetCreator \
-R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /data/resources/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /data/resources/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-known /data/resources/human/cosmic/b37/cosmic_78.indels.b37.vcf \
-I "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam \
-o "$seqId"_"$sampleId"_indel_realigned.intervals \
-L /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip "$padding" \
-nt 10 \
-dt NONE

#Realign around indels
$GATK IndelRealigner \
-R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /data/resources/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /data/resources/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-known /data/resources/human/cosmic/b37/cosmic_78.indels.b37.vcf \
-targetIntervals "$seqId"_"$sampleId"_indel_realigned.intervals \
--maxReadsForRealignment 500000 \
--maxConsensuses 750 \
--maxReadsForConsensuses 3000 \
--maxReadsInMemory 3750000 \
-LOD 0.4 \
-I "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam \
-o "$seqId"_"$sampleId"_indel_realigned.bam \
-dt NONE

#soft clip PCR primers
$SOFTCLIP \
-I "$seqId"_"$sampleId"_indel_realigned.bam \
-O "$seqId"_"$sampleId"_clipped.bam \
-T /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed

#sort and index BAM
$SINGULARITY samtools sort -@5 -m16G -o "$seqId"_"$sampleId"_clipped_sorted.bam "$seqId"_"$sampleId"_clipped.bam
$SINGULARITY samtools index "$seqId"_"$sampleId"_clipped_sorted.bam

#fix bam tags
$PICARD SetNmMdAndUqTags \
I="$seqId"_"$sampleId"_clipped_sorted.bam \
O="$seqId"_"$sampleId".bam \
CREATE_INDEX=true \
IS_BISULFITE_SEQUENCE=false \
R=/data/resources/human/mappers/b37/bwa/human_g1k_v37.fasta

### Variant calling ###

#make bai alias for Pisces
ln -s "$seqId"_"$sampleId".bai "$seqId"_"$sampleId".bam.bai

#extract thick regions
awk '{print $1"\t"$7"\t"$8}' /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed | \
$SINGULARITY bedtools merge > "$panel"_ROI_b37_thick.bed

#Call somatic variants
$PISCES \
--rmxnfilter 5,9,0.05 \
-b ./"$seqId"_"$sampleId".bam \
-g /data/resources/human/gatk/2.8/b37/ \
--minvf 0.01 \
--ssfilter false \
--minbq 20 \
--maxvq 100 \
-c 50 \
--sbfilter 0.5 \
--minvq 20 \
--vqfilter 30 \
--gvcf false \
--callmnvs true \
--maxmnvlength 100 \
-o .

#fix VCF name
echo "$sampleId" > name
$SINGULARITY /bcftools-1.2/bcftools reheader \
-s name \
-o "$seqId"_"$sampleId"_fixed.vcf \
"$seqId"_"$sampleId".vcf
rm name

#left align and trim variants
$GATK LeftAlignAndTrimVariants \
-R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_left_aligned.vcf \
-V "$seqId"_"$sampleId"_fixed.vcf \
-L "$panel"_ROI_b37_thick.bed \
-dt NONE

#Annotate with GATK contextual information
$GATK VariantAnnotator \
-R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId".bam \
-V "$seqId"_"$sampleId"_left_aligned.vcf \
-L "$panel"_ROI_b37_thick.bed \
-o "$seqId"_"$sampleId"_left_aligned_annotated.vcf \
-A BaseQualityRankSumTest -A ChromosomeCounts -A MappingQualityRankSumTest -A MappingQualityZero -A RMSMappingQuality \
-dt NONE

#Annotate with low complexity region length using mdust
$SINGULARITY bcftools annotate \
-a /data/resources/human/gatk/2.8/b37/human_g1k_v37.mdust.v34.lpad1.bed.gz \
-c CHROM,FROM,TO,LCRLen \
-h <(echo '##INFO=<ID=LCRLen,Number=1,Type=Integer,Description="Overlapping mdust low complexity region length (mask cutoff: 34)">') \
-o "$seqId"_"$sampleId"_lcr.vcf \
"$seqId"_"$sampleId"_left_aligned_annotated.vcf

#Filter variants
$GATK VariantFiltration \
-R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_"$sampleId"_lcr.vcf \
--filterExpression "LCRLen > 8" \
--filterName "LowComplexity" \
--filterExpression "DP < 50" \
--filterName "LowDP" \
-L "$panel"_ROI_b37_thick.bed \
-o "$seqId"_"$sampleId"_filtered.vcf \
-dt NONE

### QC ###

#Convert BED to interval_list for later
$PICARD BedToIntervalList \
I="$panel"_ROI_b37_thick.bed \
O="$panel"_ROI.interval_list \
SD=/data/resources/human/gatk/2.8/b37/human_g1k_v37.dict

#HsMetrics: capture & pooling performance
$PICARD CollectHsMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_hs_metrics.txt \
R=/data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
BAIT_INTERVALS="$panel"_ROI.interval_list \
TARGET_INTERVALS="$panel"_ROI.interval_list

#Generate per-base coverage: variant detection sensitivity
$GATK DepthOfCoverage \
-R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_DepthOfCoverage \
-I "$seqId"_"$sampleId".bam \
-L "$panel"_ROI_b37_thick.bed \
--countType COUNT_FRAGMENTS \
--minMappingQuality 20 \
--minBaseQuality 20 \
--omitIntervalStatistics \
-ct "$minimumCoverage" \
-nt 10 \
-dt NONE

# generate tabix index for depth of coverage
sed 's/:/\t/g' "$seqId"_"$sampleId"_DepthOfCoverage \
    | grep -v "^Locus" \
    | sort -k1,1 -k2,2n \
    | $SINGULARITY bgzip > "$seqId"_"$sampleId"_DepthOfCoverage.gz

$SINGULARITY tabix -b 2 -e 2 -s 1 "$seqId"_"$sampleId"_DepthOfCoverage.gz


#Calculate gene (clinical) percentage coverage
$COVERAGE \
"$seqId"_"$sampleId"_DepthOfCoverage \
/data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_genes.txt \
/data/resources/human/refseq/ref_GRCh37.p13_top_level.gff3 \
-p5 \
-d"$minimumCoverage" \
> "$seqId"_"$sampleId"_PercentageCoverage.txt

#Gather QC metrics
totalReads=$(head -n8 "$seqId"_"$sampleId"_hs_metrics.txt | tail -n1 | cut -s -f6) #The total number of reads in the SAM or BAM file examine.
pctSelectedBases=$(head -n8 "$seqId"_"$sampleId"_hs_metrics.txt | tail -n1 | cut -s -f19) #On+Near Bait Bases / PF Bases Aligned.
totalTargetedUsableBases=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f2) #total number of usable bases.
meanOnTargetCoverage=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f3) #avg usable coverage
pctTargetBasesCt=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f7) #percentage panel covered with good enough data for variant detection

#Print QC metrics
echo -e "TotalReads\tRawSequenceQuality\tTotalTargetUsableBases\tPctSelectedBases\tPctTargetBasesCt\tMeanOnTargetCoverage" > "$seqId"_"$sampleId"_QC.txt
echo -e "$totalReads\t$rawSequenceQuality\t$totalTargetedUsableBases\t$pctSelectedBases\t$pctTargetBasesCt\t$meanOnTargetCoverage" >> "$seqId"_"$sampleId"_QC.txt

#Add VCF meta data to final VCF
grep '^##' "$seqId"_"$sampleId"_filtered.vcf > "$seqId"_"$sampleId"_filtered_meta.vcf
echo \#\#SAMPLE\=\<ID\="$sampleId",Tissue\=Somatic,WorklistId\="$worklistId",SeqId\="$seqId",Assay\="$panel",PipelineName\=SomaticAmplicon,PipelineVersion\="$version",RawSequenceQuality\="$rawSequenceQuality",TotalReads\="$totalReads",PctSelectedBases\="$pctSelectedBases",MeanOnTargetCoverage\="$meanOnTargetCoverage",PctTargetBasesCt\="$pctTargetBasesCt",TotalTargetedUsableBases\="$totalTargetedUsableBases",RemoteVcfFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId"_filtered_meta.vcf),RemoteBamFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId".bam)\> >> "$seqId"_"$sampleId"_filtered_meta.vcf
grep -v '^##' "$seqId"_"$sampleId"_filtered.vcf >> "$seqId"_"$sampleId"_filtered_meta.vcf

#Variant Evaluation
$GATK VariantEval \
-R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_variant_evaluation.txt \
--eval:"$seqId"_"$sampleId" "$seqId"_"$sampleId"_filtered_meta.vcf \
--comp:omni2.5 /data/resources/human/gatk/2.8/b37/1000G_omni2.5.b37.vcf \
--comp:hapmap3.3 /data/resources/human/gatk/2.8/b37/hapmap_3.3.b37.vcf \
--comp:cosmic78 /data/resources/human/cosmic/b37/cosmic_78.b37.vcf \
-L "$panel"_ROI_b37_thick.bed \
-nt 20 \
-dt NONE

$SINGULARITY vep \
--verbose \
--no_progress \
--everything \
--fork 10 \
--species homo_sapiens \
--assembly GRCh37 \
--input_file "$seqId"_"$sampleId"_filtered_meta.vcf \
--format vcf \
--output_file "$seqId"_"$sampleId"_filtered_meta_annotated.vcf \
--force_overwrite \
--no_stats \
--cache \
--dir /data/resources/human/vep-cache/refseq37_v97 \
--fasta /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
--no_intergenic \
--offline \
--cache_version 97 \
--allele_number \
--no_escape \
--shift_hgvs 1 \
--vcf \
--refseq

#check VEP has produced annotated VCF
if [ ! -e "$seqId"_"$sampleId"_filtered_meta_annotated.vcf ]; then
    cp "$seqId"_"$sampleId"_filtered_meta.vcf "$seqId"_"$sampleId"_filtered_meta_annotated.vcf
fi

#index & validate final VCF
$GATK ValidateVariants \
-R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_"$sampleId"_filtered_meta_annotated.vcf \
-dt NONE

$VCFPARSE \
--transcripts /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_PreferredTranscripts.txt \
--transcript_strictness low \
--known_variants /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_KnownVariants.vcf \
--config /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ReportConfig.txt \
--filter_non_pass \
"$seqId"_"$sampleId"_filtered_meta_annotated.vcf

mv "$sampleId"_VariantReport.txt "$seqId"_"$sampleId"_VariantReport.txt

### Optional steps ###

# custom coverage reporting
if [ $custom_coverage == true ]; then
    ################# Path change to hscoverage_outdir #################
    hscoverage_outdir=hotspot_coverage
    mkdir $hscoverage_outdir

    for bedFile in $(ls /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_coverage/*.bed); do

        #extract target name
        target=$(basename "$bedFile" | sed 's/\.bed//g')

        echo $target

        # calculate coverage
        $COVERCALC \
            -B $bedFile \
            -D "$seqId"_"$sampleId"_DepthOfCoverage.gz \
            -d $minimumCoverage \
            -p 0 \
            -g /data/diagnostics/pipelines/$pipelineName/$pipelineName-$pipelineVersion/$panel/hotspot_coverage/"$target".groups \
            -o "$seqId"_"$sampleId"_"$target" \
            -O "$hscoverage_outdir"/

        # remove header from gaps file
        if [[ $(wc -l < $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".gaps) -eq 1 ]]; then
            
            # no gaps
            touch $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".nohead.gaps
        else
            # gaps
            grep -v '^#' $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".gaps > $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".nohead.gaps
        fi

        
        rm $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".gaps
    done

    for gapsFile in $hscoverage_outdir/*nohead.gaps; do

        name=$(echo $(basename $gapsFile) | cut -d"." -f1)
        echo $name

        $BED \
        --bedfile $gapsFile \
        --outname "$name".gaps \
        --outdir "$hscoverage_outdir"/ \
        --preferred_tx /data/diagnostics/pipelines/$pipelineName/"$pipelineName"-"$pipelineVersion"/"$panel"/"$panel"_PreferredTranscripts.txt

        rm $hscoverage_outdir/"$name".nohead.gaps


    done

    # combine all total coverage files
    if [ -f $hscoverage_outdir/"$seqId"_"$sampleId"_coverage.txt ]; then rm $hscoverage_outdir/"$seqId"_"$sampleId"_coverage.txt; fi
    cat $hscoverage_outdir/*.totalCoverage | grep "FEATURE" | head -n 1 >> $hscoverage_outdir/"$seqId"_"$sampleId"_coverage.txt
    cat $hscoverage_outdir/*.totalCoverage | grep -v "FEATURE" | grep -vP "combined_\\S+_GENE" >> $hscoverage_outdir/"$seqId"_"$sampleId"_coverage.txt
    rm $hscoverage_outdir/*.totalCoverage
fi

# custom variant reporting
if [ $custom_variants == true ]; then
    mkdir hotspot_variants



    for bedFile in $(ls /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_variants/*.bed); do

        # extract target name
        target=$(basename "$bedFile" | sed 's/\.bed//g')

        # select variants
        $GATK VariantFiltration \
          -R /data/resources/human/gatk/2.8/b37/human_g1k_v37.fasta \
          -V "$seqId"_"$sampleId"_filtered_meta_annotated.vcf \
          -L "$bedFile" \
          -o hotspot_variants/"$seqId"_"$sampleId"_"$target"_filtered_meta_annotated.vcf \
          -dt NONE
    
        # write targeted dataset to table using vcf_parse python utility
        
        $VCFPARSE \
          --transcripts /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_PreferredTranscripts.txt \
          --transcript_strictness low \
          --known_variants /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_KnownVariants.vcf \
          --config /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ReportConfig.txt \
          --filter_non_pass \
          hotspot_variants/"$seqId"_"$sampleId"_"$target"_filtered_meta_annotated.vcf
        

        # move to hotspot_variants
        mv "$sampleId"_VariantReport.txt hotspot_variants/"$seqId"_"$sampleId"_"$target"_VariantReport.txt

    done



fi

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
        cat $(ls /data/output/results/"$seqId"/"$panel"/*/*VariantReport.txt | head -n1) | head -n1 > /data/output/results/"$seqId"/"$panel"/"$seqId"_merged_variant_report.txt
        echo -e "Sample\tBRCA1_500X\tBRCA2_500X\tBRCA1_100X\tBRCA2_100X" > /data/output/results/"$seqId"/"$panel"/"$seqId"_merged_coverage_report.txt

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
