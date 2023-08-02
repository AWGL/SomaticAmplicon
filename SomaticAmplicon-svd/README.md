# SomaticAmplicon

NGS pipeline for somatic variant calling from amplicon datasets

Supports Illumina paired-end only and must have resonable overlap (>10bp)

Uses BWA for Alignment, Pisces for Variant Calling and VEP for annotation.

Executables are stored as singularity containers (see envs/)


## Run

Directory structure required for pipeline

```
/data
└── results
     └── seqId
         ├── panel1
         │   ├── sample1
         │   ├── sample2
         │   └── sample3
         └── panel2
             ├── sample1
             ├── sample2
             └── sample3
```
Script 1 runs in sample folder, requires fastq files split by lane

```
sbatch 1_SomaticAmplicon.sh

```

## Assays

Validated for TruSightMyeloid, NGHS-101X (CRM), NGHS-102X (BRCA)

## Configuration

The variables files within each panel folder e.g. NGHS-101X should contain variables stating whether or not to execute optional section of the pipeline. For example:

```
custom_coverage=false
custom_variants=true
generate_worksheets=false
merge_reports=false

```


## Resources

Reference Genome: human_g1k_v37.fasta

SNPs/Indels for Realignment: 1000G_phase1.indels.b37.vcf, Mills_and_1000G_gold_standard.indels.b37.vcf and cosmic_78.indels.b37.vcf

SNPS/Indels for Evaluation: 1000G_omni2.5.b37.vcf, hapmap_3.3.b37.vcf cosmic_78.b37.vcf

VEP Cache: refseq37_v97

## References

[1] https://github.com/Illumina/Pisces





