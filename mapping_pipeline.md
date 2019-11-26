

# Mapping with HiSat2

## Create HiSat2 reference index

```
/data/kelley/projects/programs/hisat2-2.1.0/hisat2-build -p 10 \
--ss /data/kelley/projects/bear_RNA/toBrownBear/ref-index/splice_sites_genomic.txt \
--exon /data/kelley/projects/bear_RNA/toBrownBear/ref-index/exons_genomic.txt \
-f /data/kelley/projects/bear_RNA/Uarctoshorribilis/GCF_003584765.1_ASM358476v1_genomic.fna /data/kelley/projects/bear_RNA/Uarctoshorribilis/Uarctoshorribilis-genomic
```

## Mapping
```
hisat2 --phred33 -k 10 --met-file stats/${file_name}.stats --rg-id ${file_name} --rg SM:${file_name} --rg PL:illumina \
-p 6 --rna-strandness RF --fr --dta --un-conc-gz unmapped/${file_name}.unmapped \
-x /data/kelley/projects/bear_RNA/Uarctoshorribilis/Uarctoshorribilis-updated -1 ${read1} -2 ${read2} -S /scratch/joanna.l.kelley_54456/test/${file_name}.sam
```

# Fix mate information 

```
java -Xmx6G -jar /data/cornejo/projects/programs/picard-tools-2.2.1/picard.jar FixMateInformation INPUT=${sam_in} \
OUTPUT=/scratch/joanna.l.kelley_54456/${file_name}.fix.sam VALIDATION_STRINGENCY=SILENT

samtools view -@ 10 -b -h /scratch/joanna.l.kelley_54456/${file_name}.fix.sam > /scratch/joanna.l.kelley_54456/${file_name}.bam

# writes output to my folder whereas the rest is written to scratch
samtools sort -m 4G -o results/${file_name}_sorted_fixed.bam -O bam -T ${file_name} -@ 10 /scratch/joanna.l.kelley_54456/${file_name}.bam
```

# Collect mapping stats

```
java -Xmx6G -jar /data/cornejo/projects/programs/picard-tools-2.2.1/picard.jar CollectRnaSeqMetrics \
REF_FLAT=~/refFlat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
I=results/${file_name}_sorted_fixed.bam \
O=rnaseqmetrics/${file_name}_rnaalignmentsummarymetrics

java -Xmx6G -jar /data/cornejo/projects/programs/picard-tools-2.2.1/picard.jar CollectAlignmentSummaryMetrics \
R=/data/kelley/projects/bear_RNA/Uarctoshorribilis/GCF_003584765.1_ASM358476v1_genomic.fna \
I=results/${file_name}_sorted_fixed.bam \
O=stats/${file_name}_alignmentsummarymetrics
```

# Stringtie version 1.3.5

```
stringtie ${sam1} -o output-forCounts/${file_name}.gtf -m 50 --rf -e -B -c 10 -p 2 \
-G /data/kelley/projects/bear_RNA/Uarctoshorribilis/GCF_003584765.1_ASM358476v1_modified_mito.gff
```

# Run prepde.py 
(distributed as part of the stringtie package) to generate gene and transcript count matrices. 
