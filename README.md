# SomaticAmplicon

NGS pipeline for somatic variant calling from amplicon datasets

Supports Illumina paired-end only and must have resonable overlap (>10bp)

# Run

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




