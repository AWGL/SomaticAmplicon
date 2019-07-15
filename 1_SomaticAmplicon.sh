#!/bin/bash

#PBS -l walltime=20:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Somatic Amplicon Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="1.7.7"

# Directory structure required for pipeline
#
# /data
# └── results
#     └── seqId
#         ├── panel1
#         │   ├── sample1
#         │   ├── sample2
#         │   └── sample3
#         └── panel2
#             ├── sample1
#             ├── sample2
#             └── sample3
#
# Script 1 runs in sample folder, requires fastq files split by lane

countQCFlagFails() {
    #count how many core FASTQC tests failed
    grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    grep -v ^WARN | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}

#load sample & pipeline variables
. *.variables
. /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel".variables

### Preprocessing ###

#record FASTQC pass/fail
rawSequenceQuality=PASS

#convert FASTQ to uBAM & add RGIDs
for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

    #trim adapters
    /share/apps/cutadapt-distros/cutadapt-1.9.1/bin/cutadapt \
    -a "$read1Adapter" \
    -A "$read2Adapter" \
    -m 50 \
    -o "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -p "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    "$read1Fastq" \
    "$read2Fastq"

    #merge overlapping reads
    /share/apps/pear-distros/pear-0.9.10-bin-64/pear-0.9.10-bin-64 \
    -f "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -r "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    -o "$seqId"_"$sampleId"_"$laneId"_merged.fastq \
    -j 12

    #convert fastq to ubam
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.2/picard.jar FastqToSam \
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
    TMP_DIR=/state/partition1/tmpdir

    #fastqc
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R1.fastq
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    #check FASTQC output
    if [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R1_fastqc/summary.txt) -gt 0 ] || [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R2_fastqc/summary.txt) -gt 0 ]; then
        rawSequenceQuality=FAIL
    fi

    #clean up
    rm "$seqId"_"$sampleId"_"$laneId"_R1.fastq "$seqId"_"$sampleId"_"$laneId"_R2.fastq "$seqId"_"$sampleId"_"$laneId"_merged.fastq.*

done

#merge lane bams
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.2/picard.jar MergeSamFiles \
$(ls "$seqId"_"$sampleId"_*_unaligned.bam | sed 's/^/I=/' | tr '\n' ' ') \
SORT_ORDER=queryname \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT \
USE_THREADING=true \
MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=/state/partition1/tmpdir \
O="$seqId"_"$sampleId"_unaligned.bam

#uBam2fq, map & MergeBamAlignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.2/picard.jar SamToFastq \
I="$seqId"_"$sampleId"_unaligned.bam \
FASTQ=/dev/stdout \
NON_PF=true \
MAX_RECORDS_IN_RAM=2000000 \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=/state/partition1/tmpdir | \
/share/apps/bwa-distros/bwa-0.7.15/bwa mem \
-M \
-t 12 \
-p \
/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
/dev/stdin | \
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.2/picard.jar MergeBamAlignment \
ATTRIBUTES_TO_RETAIN=X0 \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM="$seqId"_"$sampleId"_unaligned.bam \
OUTPUT="$seqId"_"$sampleId"_aligned.bam \
R=/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
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
TMP_DIR=/state/partition1/tmpdir

#Realign soft clipped bases
/share/apps/jre-distros/jre1.8.0_101/bin/java -Xmx2g -jar /data/diagnostics/apps/AmpliconRealigner/AmpliconRealigner-1.1.1.jar \
-I "$seqId"_"$sampleId"_aligned.bam \
-O "$seqId"_"$sampleId"_amplicon_realigned.bam \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed

#sort and index BAM
/share/apps/samtools-distros/samtools-1.3.1/samtools sort -@8 -m8G -o "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam "$seqId"_"$sampleId"_amplicon_realigned.bam
/share/apps/samtools-distros/samtools-1.3.1/samtools index "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam

#left align indels
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T LeftAlignIndels \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam \
-o "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam \
-dt NONE

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx24g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-known /state/partition1/db/human/cosmic/b37/cosmic_78.indels.b37.vcf \
-I "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam \
-o "$seqId"_"$sampleId"_indel_realigned.intervals \
-L /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip "$padding" \
-nt 12 \
-dt NONE

#Realign around indels
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx24g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-known /state/partition1/db/human/cosmic/b37/cosmic_78.indels.b37.vcf \
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
/share/apps/jre-distros/jre1.8.0_101/bin/java -Xmx2g -jar /data/diagnostics/apps/SoftClipPCRPrimer/SoftClipPCRPrimer-1.1.0.jar \
-I "$seqId"_"$sampleId"_indel_realigned.bam \
-O "$seqId"_"$sampleId"_clipped.bam \
-T /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed

#sort and index BAM
/share/apps/samtools-distros/samtools-1.3.1/samtools sort -@8 -m8G -o "$seqId"_"$sampleId"_clipped_sorted.bam "$seqId"_"$sampleId"_clipped.bam
/share/apps/samtools-distros/samtools-1.3.1/samtools index "$seqId"_"$sampleId"_clipped_sorted.bam

#fix bam tags
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.2/picard.jar SetNmMdAndUqTags \
I="$seqId"_"$sampleId"_clipped_sorted.bam \
O="$seqId"_"$sampleId".bam \
CREATE_INDEX=true \
IS_BISULFITE_SEQUENCE=false \
R=/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta

### Variant calling ###

#make bai alias for Pisces
ln -s "$seqId"_"$sampleId".bai "$seqId"_"$sampleId".bam.bai

#extract thick regions
awk '{print $1"\t"$7"\t"$8}' /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge > "$panel"_ROI_b37_thick.bed

#load mono
. /opt/mono/env.sh

#Call somatic variants
mono /share/apps/MiSeqReporter-distros/MiSeqReporter-2.6.3/CallSomaticVariants.exe \
-B ./"$seqId"_"$sampleId".bam \
-g /data/db/human/gatk/2.8/b37 \
-f 0.01 \
-fo False \
-b 20 \
-q 100 \
-c 50 \
-s 0.5 \
-a 20 \
-F 30 \
-gVCF False \
-i false \
-PhaseSNPs true \
-MaxPhaseSNPLength 100 \
-r .

#fix VCF name
echo "$sampleId" > name
/share/apps/bcftools-distros/bcftools-1.2/bcftools reheader \
-s name \
-o "$seqId"_"$sampleId"_fixed.vcf \
$(echo "$seqId"_"$sampleId" | sed 's/_/-/g')_S999.vcf
rm name

#left align and trim variants
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T LeftAlignAndTrimVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_left_aligned.vcf \
-V "$seqId"_"$sampleId"_fixed.vcf \
-L "$panel"_ROI_b37_thick.bed \
-dt NONE

#Annotate with GATK contextual information
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId".bam \
-V "$seqId"_"$sampleId"_left_aligned.vcf \
-L "$panel"_ROI_b37_thick.bed \
-o "$seqId"_"$sampleId"_left_aligned_annotated.vcf \
-A BaseQualityRankSumTest -A ChromosomeCounts -A MappingQualityRankSumTest -A MappingQualityZero -A RMSMappingQuality \
-dt NONE

#Annotate with low complexity region length using mdust
/share/apps/bcftools-distros/bcftools-1.3.1/bcftools annotate \
-a /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.mdust.v34.lpad1.bed.gz \
-c CHROM,FROM,TO,LCRLen \
-h <(echo '##INFO=<ID=LCRLen,Number=1,Type=Integer,Description="Overlapping mdust low complexity region length (mask cutoff: 34)">') \
-o "$seqId"_"$sampleId"_lcr.vcf \
"$seqId"_"$sampleId"_left_aligned_annotated.vcf

#Filter variants
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
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
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.2/picard.jar BedToIntervalList \
I="$panel"_ROI_b37_thick.bed \
O="$panel"_ROI.interval_list \
SD=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.dict

#HsMetrics: capture & pooling performance
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.2/picard.jar CollectHsMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_hs_metrics.txt \
R=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
BAIT_INTERVALS="$panel"_ROI.interval_list \
TARGET_INTERVALS="$panel"_ROI.interval_list

#Generate per-base coverage: variant detection sensitivity
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx12g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_DepthOfCoverage \
-I "$seqId"_"$sampleId".bam \
-L "$panel"_ROI_b37_thick.bed \
--countType COUNT_FRAGMENTS \
--minMappingQuality 20 \
--minBaseQuality 20 \
--omitIntervalStatistics \
-ct "$minimumCoverage" \
-nt 12 \
-dt NONE


# generate tabix index for depth of coverage
sed 's/:/\t/g' /data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId"_DepthOfCoverage \
    | grep -v "^Locus" \
    | sort -k1,1 -k2,2n \
    | /share/apps/htslib-distros/htslib-1.4.1/bgzip > /data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId"_DepthOfCoverage.gz

/share/apps/htslib-distros/htslib-1.4.1/tabix -b 2 -e 2 -s 1 /data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId"_DepthOfCoverage.gz


#Calculate gene (clinical) percentage coverage
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /data/diagnostics/apps/CoverageCalculator-2.0.2/CoverageCalculator-2.0.2.jar \
"$seqId"_"$sampleId"_DepthOfCoverage \
/data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_genes.txt \
/state/partition1/db/human/refseq/ref_GRCh37.p13_top_level.gff3 \
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
echo -e "TotalReads\tRawSequenceQuality\tTotalTargetUsableBases\tPctSelectedBases\tPctTargetBasesCt\tMeanOnTargetCoverage" > "$seqId"_"$sampleId"_qc.txt
echo -e "$totalReads\t$rawSequenceQuality\t$totalTargetedUsableBases\t$pctSelectedBases\t$pctTargetBasesCt\t$meanOnTargetCoverage" >> "$seqId"_"$sampleId"_qc.txt

#Add VCF meta data to final VCF
grep '^##' "$seqId"_"$sampleId"_filtered.vcf > "$seqId"_"$sampleId"_filtered_meta.vcf
echo \#\#SAMPLE\=\<ID\="$sampleId",Tissue\=Somatic,WorklistId\="$worklistId",SeqId\="$seqId",Assay\="$panel",PipelineName\=SomaticAmplicon,PipelineVersion\="$version",RawSequenceQuality\="$rawSequenceQuality",TotalReads\="$totalReads",PctSelectedBases\="$pctSelectedBases",MeanOnTargetCoverage\="$meanOnTargetCoverage",PctTargetBasesCt\="$pctTargetBasesCt",TotalTargetedUsableBases\="$totalTargetedUsableBases",RemoteVcfFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId"_filtered_meta.vcf),RemoteBamFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId".bam)\> >> "$seqId"_"$sampleId"_filtered_meta.vcf
grep -v '^##' "$seqId"_"$sampleId"_filtered.vcf >> "$seqId"_"$sampleId"_filtered_meta.vcf

#Variant Evaluation
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantEval \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_variant_evaluation.txt \
--eval:"$seqId"_"$sampleId" "$seqId"_"$sampleId"_filtered_meta.vcf \
--comp:omni2.5 /state/partition1/db/human/gatk/2.8/b37/1000G_omni2.5.b37.vcf \
--comp:hapmap3.3 /state/partition1/db/human/gatk/2.8/b37/hapmap_3.3.b37.vcf \
--comp:cosmic78 /state/partition1/db/human/cosmic/b37/cosmic_78.b37.vcf \
-L "$panel"_ROI_b37_thick.bed \
-nt 12 \
-dt NONE

### Reporting ###

#annotate VCF with VEP
perl /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl \
--verbose \
--no_progress \
--everything \
--fork 12 \
--species homo_sapiens \
--assembly GRCh37 \
--input_file "$seqId"_"$sampleId"_filtered_meta.vcf \
--format vcf \
--output_file "$seqId"_"$sampleId"_filtered_meta_annotated.vcf \
--force_overwrite \
--no_stats \
--cache \
--dir /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
--fasta /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
--no_intergenic \
--offline \
--cache_version 86 \
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
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T ValidateVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_"$sampleId"_filtered_meta_annotated.vcf \
-dt NONE


#custom coverage reporting
if [ -d /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_coverage ]; then
    
    hscoverage_outdir=/data/results/$seqId/$panel/$sampleId/hotspot_coverage/
    mkdir $hscoverage_outdir

    for bedFile in $(ls /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_coverage/*.bed); do

        #extract target name
        target=$(basename "$bedFile" | sed 's/\.bed//g')

        echo $target

        # calcualte coverage
        source /home/transfer/miniconda3/bin/activate CoverageCalculatorPy

        python /home/transfer/pipelines/CoverageCalculatorPy/CoverageCalculatorPy.py \
            -B $bedFile \
            -D /data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId"_DepthOfCoverage.gz \
            --depth $minimumCoverage \
            --padding 0 \
            --groupfile /data/diagnostics/pipelines/$pipelineName/$pipelineName-$pipelineVersion/$panel/hotspot_coverage/"$target".groups \
            --outname "$seqId"_"$sampleId"_"$target" \
            --outdir $hscoverage_outdir

        # remove header from gaps file
        if [[ $(wc -l < $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".gaps) -eq 1 ]]; then
            
            # no gaps
            touch $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".nohead.gaps
        else
            # gaps
            grep -v '^#' $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".gaps > $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".nohead.gaps
        fi

        rm $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".gaps

        source /home/transfer/miniconda3/bin/deactivate

        #annotate the gaps with HGVS nomenclature using bed2hgvs.py
        source /home/transfer/miniconda3/bin/activate bed2hgvs

        python /data/diagnostics/apps/bed2hgvs/bed2hgvs-0.1.1/bed2hgvs.py --config /data/diagnostics/apps/bed2hgvs/bed2hgvs-0.1.1/configs/cluster.yaml \
            --input $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".nohead.gaps \
            --output $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".gaps \
            --transcript_map /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_PreferredTranscripts.txt

        rm $hscoverage_outdir/"$seqId"_"$sampleId"_"$target".nohead.gaps

        source /home/transfer/miniconda3/bin/deactivate

    done
fi

# combine all total coverage files
if [  -d /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_coverage ]; then 
    rm $hscoverage_outdir/"$seqId"_"$sampleId"_coverage.txt
    cat $hscoverage_outdir/*.totalCoverage | grep "FEATURE" | head -n 1 >> $hscoverage_outdir/"$sampleId"_coverage.txt
    cat $hscoverage_outdir/*.totalCoverage | grep -v "FEATURE" | grep -vP "combined_\\S+_GENE" >> $hscoverage_outdir/"$sampleId"_coverage.txt
    rm $hscoverage_outdir/*.totalCoverage
fi

#custom variant reporting
if [ -d /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_variants ]; then
    mkdir hotspot_variants

    for bedFile in $(ls /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_variants/*.bed); do

        #extract target name
        target=$(basename "$bedFile" | sed 's/\.bed//g')

        #select variants
        /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
        -V "$seqId"_"$sampleId"_filtered_meta_annotated.vcf \
        -L "$bedFile" \
        -o hotspot_variants/"$seqId"_"$sampleId"_"$target"_filtered_meta_annotated.vcf \
        -dt NONE
	
	#write targeted dataset to table using vcf_parse python utility

	source /home/transfer/miniconda3/bin/activate vcf_parse
	
	python /data/diagnostics/apps/vcf_parse/vcf_parse-0.1.2/vcf_parse.py \
	--transcripts /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_PreferredTranscripts.txt \
	--transcript_strictness low \
	--known_variants /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_KnownVariants.vcf \
	--config /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ReportConfig.txt \
        --filter_non_pass \
	hotspot_variants/"$seqId"_"$sampleId"_"$target"_filtered_meta_annotated.vcf
	
	source /home/transfer/miniconda3/bin/deactivate

        #move to hotspot_variants
        mv "$sampleId"_VariantReport.txt hotspot_variants/"$seqId"_"$sampleId"_"$target"_VariantReport.txt

    done
fi

#write big vcf dataset to table using vcf_parse python utility

source /home/transfer/miniconda3/bin/activate vcf_parse

python /data/diagnostics/apps/vcf_parse/vcf_parse-0.1.2/vcf_parse.py \
--transcripts /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_PreferredTranscripts.txt \
--transcript_strictness low \
--known_variants /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_KnownVariants.vcf \
--config /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ReportConfig.txt \
--filter_non_pass \
"$seqId"_"$sampleId"_filtered_meta_annotated.vcf

mv "$sampleId"_VariantReport.txt "$seqId"_"$sampleId"_VariantReport.txt

source /home/transfer/miniconda3/bin/deactivate

# Merge QC files
python /data/diagnostics/scripts/merge_qc_files.py ..


# Generate Analysis Workseets
# This block should only be carried out when all samples for the
# panel have been processed

# number of samples to be processed (i.e. count variables files)
expected=$(for i in /data/results/$seqId/$panel/*/*.variables; do echo $i; done | wc -l)

# number of samples that have completed
complete=$(for i in /data/results/$seqId/$panel/*/*VariantReport.txt; do echo $i; done | wc -l)

if [ $complete -eq $expected ]; then
   
   source ~/miniconda3/bin/activate VirtualHood

   # identify name of NTC
   ntc=$(for s in /data/results/$seqId/$panel/*/; do echo $(basename $s);done | grep 'NTC')

   # loop over all samples and generate a report
   for s in /data/results/$seqId/$panel/*/; do
       
       # clear previous instance
       unset referral 
       
       # set variables
       sample=$(basename $s)
       . /data/results/$seqId/$panel/$sample/*.variables

       # check that referral vraible is defined, if not set as NA
       if [ -z $referral ];then
           referral=NA
       fi

       # do not generate report where NTC is the query sample
       if [ $sample != $ntc ]; then
           
           if [ $referral == 'FOCUS4' ] || [ $referral == 'GIST' ] || [ $referral == 'iNATT' ];then
               python /data/diagnostics/apps/VirtualHood/CRM_report.py $seqId $sample $worklistId $referral $ntc
           elif [ $referral == 'Melanoma' ] || [ $referral == 'Lung' ] || [ $referral == 'Colorectal' ] || [ $referral == 'Glioma' ] || [ $referral == 'Tumour' ];then
               python /data/diagnostics/apps/VirtualHood/CRM_report_new_referrals.py $seqId $sample $worklistId $referral $ntc
           fi

       fi
   done

   source ~/miniconda3/bin/deactivate
fi



### Clean up ###

#delete unused files
rm "$seqId"_"$sampleId"_*unaligned.bam "$seqId"_"$sampleId"_aligned.bam "$seqId"_"$sampleId"_aligned.bai "$seqId"_"$sampleId"_amplicon_realigned.bam
rm "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam.bai "$seqId"_"$sampleId"_indel_realigned.intervals
rm "$seqId"_"$sampleId"_clipped.bam "$seqId"_"$sampleId"_clipped_sorted.bam "$seqId"_"$sampleId"_clipped_sorted.bam.bai "$panel"_ROI.interval_list "$panel"_ROI_b37_thick.bed
rm "$seqId"_"$sampleId"_left_aligned.vcf "$seqId"_"$sampleId"_left_aligned.vcf.idx "$seqId"_"$sampleId".bam.bai "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam
rm "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bai "$seqId"_"$sampleId"_filtered_meta.vcf "$seqId"_"$sampleId"_filtered_meta.vcf.idx "$seqId"_"$sampleId"_filtered.vcf
rm "$seqId"_"$sampleId"_filtered.vcf.idx "$seqId"_"$sampleId"_fixed.vcf "$seqId"_"$sampleId"_fixed.vcf.idx "$seqId"_"$sampleId"_indel_realigned.bam "$seqId"_"$sampleId"_indel_realigned.bai
rm "$seqId"_"$sampleId"_*_fastqc.zip "$seqId"_"$sampleId"_lcr.vcf "$seqId"_"$sampleId"_lcr.vcf.idx "$seqId"_"$sampleId"_left_aligned_annotated.vcf "$seqId"_"$sampleId"_left_aligned_annotated.vcf.idx
rm $(echo "$seqId"_"$sampleId" | sed 's/_/-/g')_S999.vcf
rm -r VariantCallingLogs
