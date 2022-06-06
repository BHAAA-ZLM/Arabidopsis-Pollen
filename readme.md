# Distinct chromatin signatures in Arabidopsis male gametophyte

## Introduction
Unlike animals, the germline in plants is established de novo during post embryonic growth and may subject to somatic modifications.
The extent to which a plant would reprogram its chromatin modifications in germline remains elusive.
Here, in our recently research "Distinct chromatin signatures in Arabidopsis male gametophyte", by adapting ultrasensitive sequencing techniques, we performed ATAC-seq, RNA-seq and histone modification ChIP-seq in the sperm and the companion vegetative nuclei of the Arabidopsis pollen, and obtained some chromatin signatures in these two kinds of cells.
Here, we provide the key analytical steps and resources, including software packages, test data, and so on, which help users reproduce our results.


## ChIP-seq pipeline (using Spm_H3K4me3 as an example)
### 1.Data quality control (fastp v0.20.0)
```bash
fastp -i Spm_H3K4me3_R1.fastq.gz -I Spm_H3K4me3_R2.fastq.gz -o Spm_H3K4me3_R1.fastp.fastq.gz -O Spm_H3K4me3_R2.fastp.fastq.gz -w 8 -l 25 --detect_adapter_for_pe -j Spm_H3K4me3.fastp.json -h Spm_H3K4me3.fastp.html
```
### 2.Genome mapping (hisat2 v2.1.0, samtools v1.9)
```bash
hisat2 -p 8 -x tair10_index --no-temp-splicesite --no-spliced-alignment --summary-file Spm_H3K4me3.summary -1 Spm_H3K4me3_R1.fastp.fastq.gz -2 Spm_H3K4me3_R2.fastp.fastq.gz | samtools view -ShuF 4 -q 30 -f 2 -@ 2 - | samtools sort -@ 2 -o Spm_H3K4me3.sorted.bam -
```

### 3.Redundancy removal (picard v2.20.8)
```bash
java -Xmx5g -XX:ParallelGCThreads=8 -jar picard.jar MarkDuplicates I= Spm_H3K4me3.sorted.bam O= Spm_H3K4me3.sorted.picardMD.bam M= Spm_H3K4me3.sorted.picardMD.txt REMOVE_DUPLICATES=true 
```

### 4.Peak calling (MACS2 v2.1.4)
```bash
macs2 callpeak -c Spm_H3K4me3_input.sorted.picardMD.bam -t Spm_H3K4me3.sorted.picardMD.bam -f BAMPE -g 1.053e8 –keep-dup all -q 0.01 -n Spm_H3K4me3 --outdir macs2/ -B –SPMR
```

PS: for H3K27me3 and H3K9me2, we using “--broad --nolambda” to replace “-q 0.01”

### 5.Filter for validated peaks (HOMER v4.11.1: mergePeaks)
```bash
mergePeaks -d given Spm_H3K4me3_rep1.narrowPeak Spm_H3K4me3_rep2.narrowPeak Spm_H3K4me3_rep 3.narrowPeak -venn mergePeaks_narrow_H3K4me3_Spm_given.venn -prefix mergePeaks_narrow_H3K4me3_Spm_given -matrix mergePeaks_narrow_H3K4me3_Spm_given
```

PS: for H3K9me2, we first stitch peaks in each replicates within 5000bp as one peaks, which using “mergePeaks -d 5000”, and then filter the validated peaks from replicates

### 6.Annotate peaks (HOMER v4.11.1: annotatePeaks.pl)
```bash
annotatePeaks.pl Spm_H3K4me3.narrowPeak tair10 -gtf TAIR10.gtf -annStats Spm_H3K4me3_tair10.annStats > Spm_H3K4me3_tair10.annotate
```
### 7.Prepare bigwig files for track visualization
```bash
bdg2bw_RemoveOverlap.sh Spm_H3K4me3_treat_pileup.bdg tair10.chrom.sizes
```

## ATAC-seq pipeline (using Spm_ATAC as an example)
### 1.Data quality control (fastp v0.20.0)
```bash
fastp -i Spm_ATAC_R1.fastq.gz -I Spm_ATAC_R2.fastq.gz -o Spm_ATAC_R1.fastp.fastq.gz -O Spm_ATAC_R2.fastp.fastq.gz -w 8 -l 25 --detect_adapter_for_pe -j Spm_ATAC.fastp.json -h Spm_ATAC.fastp.html
```
### 2.Genome mapping (hisat2 v2.1.0, samtools v1.9)
```bash
hisat2 -p 8 -x tair10_index -X 2000 --no-temp-splicesite --no-spliced-alignment --summary-file Spm_ATAC.summary -1 Spm_ATAC_R1.fastp.fastq.gz -2 Spm_ATAC_R2.fastp.fastq.gz | samtools view -ShuF 4 -q 30 -f 2 -@ 2 - | samtools sort -@ 2 -o Spm_ATAC.sorted.bam -
```
### 3.Redundancy removal (picard v2.20.8)
```bash
java -Xmx5g -XX:ParallelGCThreads=8 -jar picard.jar MarkDuplicates I= Spm_ATAC.sorted.bam O= Spm_ATAC.sorted.picardMD.bam M= Spm_ATAC.sorted.picardMD.txt REMOVE_DUPLICATES=true 
```
### 4.Convert the BAM format to BED format (bedtools bamtobed v2.29.2)
```bash
bamToBed -i Spm_ATAC.sorted.picardMD.bam > Spm_ATAC.sorted.picardMD.bed
```
### 5.Peak calling (MACS2 v2.1.4)
```bash
macs2 callpeak -c SRR4000479_gDNA_ATAC.sorted.picardMD.bed -t Spm_ATAC.sorted.picardMD.bed -f BED -g 1.053e8 –keep-dup all -n Spm_ATAC --outdir macs2/ -B –SPMR --nomodel --shift -100 --extsize 200
```
### 6.Filter for validated peaks (HOMER v4.11.1: mergePeaks)
```bash
mergePeaks -d given Spm_ATAC_rep1.narrowPeak Spm_ATAC_rep2.narrowPeak -venn mergePeaks_narrow_ATAC_Spm_given.venn -prefix mergePeaks_narrow_ATAC_Spm_given -matrix mergePeaks_narrow_ATAC_Spm_given
```
### 7.Annotate peaks (HOMER v4.11.1: annotatePeaks.pl)
```bash
annotatePeaks.pl Spm_ATAC.narrowPeak tair10 -gtf TAIR10.gtf -annStats Spm_ATAC_tair10.annStats > Spm_ATAC_tair10. annotate
```
### 8.Motif enrichment analysis (HOMER v4.11.1: findMotifsGenome.pl)
```bash
findMotifsGenome.pl Spm_ATAC.narrowPeaks tair10 Spm_ATAC_motif/ -p 8
```
### 9.Prepare bigwig files for track visualization
```bash
bdg2bw_RemoveOverlap.sh Spm_ATAC_treat_pileup.bdg tair10.chrom.sizes
```

## RNA-seq pipeline (using Spm as an example)
### 1.Data quality control (fastp v0.20.0, cutadapt v2.10)
```bash
fastp -i Spm_R1.fastq.gz -I Spm_R2.fastq.gz -o Spm_R1.fastp.fastq.gz -O Spm_R2.fastp.fastq.gz -w 8 -q 20 –u 20 -j Spm.fastp.json -h Spm.fastp.html

cutadapt -j 24 -g ^ATTGCGCAATGNNNNNNNNGGG -o Spm_R1.fastp.cutadapt.fastq.gz -p Spm_R2.fastp.cutadapt.fastq.gz Spm_R1.fastp.fastq.gz Spm_R2.fastp.fastq.gz
```
### 2.Genome mapping (hisat2 v2.1.0, samtools v1.9)
```bash
hisat2 -p 8 -x tair10_index --summary-file Spm.summary --dta-cufflinks -1 Spm_R1.fastp.cutadapt.fastq.gz -2 Spm_R2.fastp.cutadapt.fastq.gz | samtools view -ShuF 4 -q 30 -f 2 -@ 2 - | samtools sort -@ 2 -o Spm.sorted.bam -
```
### 3.Convert the BAM format to bigwig format for track visualization
```bash
genomeCoverageBed -ibam Spm.sorted.bam -bg -split > Spm.bdg

bdg2bw_RemoveOverlap.sh Spm.bdg tair10.chrom.sizes
```

### 4.Get the gene expression matrix (cuffdiff v2.2.1)
```bash
cuffdiff -o cuffdiff_results/ -L Veg,Spm -p 16 -b tair10_genome.fa -u --library-type fr-unstranded TAIR10.gtf Veg.sorted.bam Spm.sorted.bam
```


## Region RPKM (calculating peaks RPKM of Spm_H3K4me3 in Spm as the example)
### 1.SampleReads: 
```bash
samtools view -c Spm_H3K4me3.sorted.picardMD.merged.bam
```

$$
SampleReadsDepth = \frac{SampleReads}{1,000,000}
$$


### 2.RegionLength

$$
RegionLength = PeakEnd - PeakStart
$$

### 3.ReadsWithinRegion
```bash
coverageBed -a Spm_H3K4me3_peaks.bed –b Spm_H3K4me3.sorted.picardMD.merged.bam > Spm_H3K4me3_coverageBed.txt
```

PS: -a should use BED (6 column) format file as input

$ReadsWithinRegion$ : the value of 7th column in Spm_H3K4me3_coverageBed.txt

### 4.The RPKM of peaks were calculated as:
$$
PeakRPKM= \frac{ReadsWithinRegion}{\frac{RegionLength}{1000} * SampleReadsDepth} 
$$

PS: users also can calculate RPKM by using script "peak_RPKM.pl" we provide:

```bash
peak_RPKM.pl -depth <SampleReadsDepth> --col <the column of reads counts,default:7> --input Spm_H3K4me3_coverageBed.txt --output Spm_H3K4me3_coverageBed_RPKM.txt
```

The RPKM of Spm_H3K4me3 peaks in Spm were saved in output file: Spm_H3K4me3_coverageBed_RPKM.txt 

## ATAC-seq tracks normalization (using Spm_ATAC as the example)
### 1.Calculate the RPKM of all the Spm_ATAC peaks;
### 2.Pick the top 1000th peak’s RPKM as standardized value;
### 3.The value in the bdg file divided by the normalized value;
### 4.Prepare bigwig files for track visualization by normalized bdg files;
