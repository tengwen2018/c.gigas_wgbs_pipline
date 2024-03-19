# c.gigas_wgbs_pipline

**1. Data trimming of Illumina 150bp PE reads**

```bash
fastp -i sample_1.fastq.gz -I sample_2.fastq.gz -o sample_clean_1.fq.gz -O sample_clean_2.fq.gz —adapter_sequence auto —detect_adapter_for_pe —unpaired1 output_um_1.fastq.gz —unpaired2 output_um_2.fastq.gz —failed_out output_failed.fastq.gz —cut_front —cut_front_window_size=1 —cut_front_mean_quality=20 —cut_tail —cut_tail_window_size=1 —cut_tail_mean_quality=20 —cut_right —cut_right_window_size=4 —cut_right_mean_quality=20 —length_required=36 —thread 1 --trim_front1 5 --trim_front2 5
```
**2. Mapping trimmed reads to reference genome and DNA methylation calling**

```bash
bsmap -a sample_clean_1.fq.gz -b sample_clean_2.fq.gz -d ref.fa -o bsmap.sample.bam -R -p 4 -n 1 -r 0 -v 0.1 -S 1 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA > bsmap.sample.bam.log && \
samtools index bsmap.sample.bam && \
pefilter -i bsmap.sample.bam -o bsmap.sample.filter.bam > pefilter.sample.log && \
mcall -m bsmap.sample.filter.bam -p 6 -r ref.fa --sampleName sample > mcall.sample.log
```
