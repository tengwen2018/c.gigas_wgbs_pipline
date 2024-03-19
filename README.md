**The pipline of DNA methylaiton data analysis of *Crassostrea gigas***

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

**3. Differential methylation analysis**

```bash
# Male VS. female in parent cohort
mcomp -d 5 -r p8m.G.bed,p10m.G.bed -r p8f.G.bed,p10f.G.bed -m male.G.bed -m female.G.bed -p 10 -c mcomp.pmale.vs.pfemale.txt —withVariance 0

# Male VS. female in offspring cohort
mcomp -d 5 -r o10d_1.G.bed,o10d_3.G.bed,o8d_1.G.bed,o8d_3.G.bed -r o10t_1.G.bed,o10t_2.G.bed,o10t_3.G.bed,o8t_1.G.bed,o8t_2.G.bed,o8t_3.G.bed,o8d_2.G.bed,o10d_2.G.bed -m male.G.bed -m female.G.bed -p 10 -c mcomp.omale.vs.ofemale.txt —withVariance 0
```
**4. Curve plot of DNA methylation profile around gene bodies**

```bash
computeMatrix scale-regions -R gene.bed -S sample.bw -b 5000 -a 5000 --regionBodyLength 10000 --binSize 100 --startLabel "TSS" --endLabel "TTS" --skipZeros --samplesLabel c.gigas -o matrix2_gene.gz --outFileNameMatrix matrix2_gene.tab --outFileSortedRegions matrix2_gene.bed && \
plotProfile -m matrix2_gene.gz \
      -out matrix2_cgigas.pdf \
      --startLabel "TSS" \
      --endLabel "TTS" \
      --perGroup \
      --yMin 0 \
      --yMax 1 \
      --outFileNameData matrix2_gene.txt
```

**5. The correlation between estimated DNA methylation levels and sequencing depth, as well as the relationship between the coverage of CpG sites and sequencing depth**

```bash
bseqc2 -i bsmap.sample.bam -o bseqc2.sample.txt -r ref.fa -l 140
```

**6. PCA analysis**

```R
library(plyr)

o10d_1 <- read.table('o10d_1.G.bed',sep='\t',header=F)
o10d_1 <- o10d_1[c('V1','V2','V3','V4','V5')]
o10d_1 <- o10d_1[o10d_1$V5>=5,]
o10d_1 <- mutate(o10d_1, pos=paste(V1,V2,V3,sep='_'))
o10d_1 <- o10d_1[c('pos', 'V4')]
colnames(o10d_1) <- c('pos', 'o10d_1')

o10d_2 <- read.table('o10d_2.G.bed',sep='\t',header=F)
o10d_2 <- o10d_2[c('V1','V2','V3','V4','V5')]
o10d_2 <- o10d_2[o10d_2$V5>=5,]
o10d_2 <- mutate(o10d_2, pos=paste(V1,V2,V3,sep='_'))
o10d_2 <- o10d_2[c('pos', 'V4')]
colnames(o10d_2) <- c('pos', 'o10d_2')

o10d_3 <- read.table('o10d_3.G.bed',sep='\t',header=F)
o10d_3 <- o10d_3[c('V1','V2','V3','V4','V5')]
o10d_3 <- o10d_3[o10d_3$V5>=5,]
o10d_3 <- mutate(o10d_3, pos=paste(V1,V2,V3,sep='_'))
o10d_3 <- o10d_3[c('pos', 'V4')]
colnames(o10d_3) <- c('pos', 'o10d_3')

o10t_1 <- read.table('o10t_1.G.bed',sep='\t',header=F)
o10t_1 <- o10t_1[c('V1','V2','V3','V4','V5')]
o10t_1 <- o10t_1[o10t_1$V5>=5,]
o10t_1 <- mutate(o10t_1, pos=paste(V1,V2,V3,sep='_'))
o10t_1 <- o10t_1[c('pos', 'V4')]
colnames(o10t_1) <- c('pos', 'o10t_1')

o10t_2 <- read.table('o10t_2.G.bed',sep='\t',header=F)
o10t_2 <- o10t_2[c('V1','V2','V3','V4','V5')]
o10t_2 <- o10t_2[o10t_2$V5>=5,]
o10t_2 <- mutate(o10t_2, pos=paste(V1,V2,V3,sep='_'))
o10t_2 <- o10t_2[c('pos', 'V4')]
colnames(o10t_2) <- c('pos', 'o10t_2')

o10t_3 <- read.table('o10t_3.G.bed',sep='\t',header=F)
o10t_3 <- o10t_3[c('V1','V2','V3','V4','V5')]
o10t_3 <- o10t_3[o10t_3$V5>=5,]
o10t_3 <- mutate(o10t_3, pos=paste(V1,V2,V3,sep='_'))
o10t_3 <- o10t_3[c('pos', 'V4')]
colnames(o10t_3) <- c('pos', 'o10t_3')

o8d_1 <- read.table('o8d_1.G.bed',sep='\t',header=F)
o8d_1 <- o8d_1[c('V1','V2','V3','V4','V5')]
o8d_1 <- o8d_1[o8d_1$V5>=5,]
o8d_1 <- mutate(o8d_1, pos=paste(V1,V2,V3,sep='_'))
o8d_1 <- o8d_1[c('pos', 'V4')]
colnames(o8d_1) <- c('pos', 'o8d_1')

o8d_2 <- read.table('o8d_2.G.bed',sep='\t',header=F)
o8d_2 <- o8d_2[c('V1','V2','V3','V4','V5')]
o8d_2 <- o8d_2[o8d_2$V5>=5,]
o8d_2 <- mutate(o8d_2, pos=paste(V1,V2,V3,sep='_'))
o8d_2 <- o8d_2[c('pos', 'V4')]
colnames(o8d_2) <- c('pos', 'o8d_2')

o8d_3 <- read.table('o8d_3.G.bed',sep='\t',header=F)
o8d_3 <- o8d_3[c('V1','V2','V3','V4','V5')]
o8d_3 <- o8d_3[o8d_3$V5>=5,]
o8d_3 <- mutate(o8d_3, pos=paste(V1,V2,V3,sep='_'))
o8d_3 <- o8d_3[c('pos', 'V4')]
colnames(o8d_3) <- c('pos', 'o8d_3')

o8t_1 <- read.table('o8t_1.G.bed',sep='\t',header=F)
o8t_1 <- o8t_1[c('V1','V2','V3','V4','V5')]
o8t_1 <- o8t_1[o8t_1$V5>=5,]
o8t_1 <- mutate(o8t_1, pos=paste(V1,V2,V3,sep='_'))
o8t_1 <- o8t_1[c('pos', 'V4')]
colnames(o8t_1) <- c('pos', 'o8t_1')

o8t_2 <- read.table('o8t_2.G.bed',sep='\t',header=F)
o8t_2 <- o8t_2[c('V1','V2','V3','V4','V5')]
o8t_2 <- o8t_2[o8t_2$V5>=5,]
o8t_2 <- mutate(o8t_2, pos=paste(V1,V2,V3,sep='_'))
o8t_2 <- o8t_2[c('pos', 'V4')]
colnames(o8t_2) <- c('pos', 'o8t_2')

o8t_3 <- read.table('o8t_3.G.bed',sep='\t',header=F)
o8t_3 <- o8t_3[c('V1','V2','V3','V4','V5')]
o8t_3 <- o8t_3[o8t_3$V5>=5,]
o8t_3 <- mutate(o8t_3, pos=paste(V1,V2,V3,sep='_'))
o8t_3 <- o8t_3[c('pos', 'V4')]
colnames(o8t_3) <- c('pos', 'o8t_3')

res <- o10d_1
for (sample in list(o10d_2,o10d_3,o10t_1,o10t_2,o10t_3,o8d_1,o8d_2,o8d_3,o8t_1,o8t_2,o8t_3)){
res <- merge(res, sample, by='pos', all=F)
}
write.table(res, 'offspring.methy.5x.txt', quote=F, sep='\t', col.names=T, row.names=F)

library("factoextra")
library("FactoMineR")

res.pca <- PCA(t(res[,2:13]), graph=FALSE)
print(res.pca)
pdf("offspring.methy.5x.eigenvalue.pdf")
eig.val <- get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
dev.off()
ind <- get_pca_ind(res.pca)
ind_coord <- ind$coord
ind_coord <- data.frame(ind_coord)
ind_coord$group <- c("o10d","o10d","o10d","o10t","o10t","o10t","o8d","o8d","o8d","o8t","o8t","o8t")
require("ggrepel") #add text
set.seed(42)
pdf("offspring.methy.5x.pca.pdf",width=4, height=4)
ggplot(ind_coord, aes(x=Dim.1, y=Dim.2, label=rownames(ind_coord), shape=group, color=group)) +
    geom_point() +
    geom_text_repel(aes(label =rownames(ind_coord)),size = 3.5) +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_hline(yintercept = 0, linetype="dashed") +
    xlab("PCA1") + ylab("PCA2") +
    scale_color_brewer(palette="Dark2") + theme_minimal() +
    theme(legend.position="none")

dev.off()

#add parents
p10m <- read.table('p10m.G.bed',sep='\t',header=F)
p10m <- p10m[c('V1','V2','V3','V4','V5')]
p10m <- p10m[p10m$V5>=5,]
p10m <- mutate(p10m, pos=paste(V1,V2,V3,sep='_'))
p10m <- p10m[c('pos', 'V4')]
colnames(p10m) <- c('pos', 'p10m')

p10f <- read.table('p10f.G.bed',sep='\t',header=F)
p10f <- p10f[c('V1','V2','V3','V4','V5')]
p10f <- p10f[p10f$V5>=5,]
p10f <- mutate(p10f, pos=paste(V1,V2,V3,sep='_'))
p10f <- p10f[c('pos', 'V4')]
colnames(p10f) <- c('pos', 'p10f')

p8m <- read.table('p8m.G.bed',sep='\t',header=F)
p8m <- p8m[c('V1','V2','V3','V4','V5')]
p8m <- p8m[p8m$V5>=5,]
p8m <- mutate(p8m, pos=paste(V1,V2,V3,sep='_'))
p8m <- p8m[c('pos', 'V4')]
colnames(p8m) <- c('pos', 'p8m')

p8f <- read.table('p8f.G.bed',sep='\t',header=F)
p8f <- p8f[c('V1','V2','V3','V4','V5')]
p8f <- p8f[p8f$V5>=5,]
p8f <- mutate(p8f, pos=paste(V1,V2,V3,sep='_'))
p8f <- p8f[c('pos', 'V4')]
p8f <- mutate(p8f, pos=paste(V1,V2,V3,sep='_'))
p8f <- p8f[c('pos', 'V4')]
colnames(p8f) <- c('pos', 'p8f')

res <- res
for (sample in list(p10m,p10f,p8m,p8f)){
res <- merge(res, sample, by='pos', all=F)
}
write.table(res, 'offspringandparent.methy.5x.txt', quote=F, sep='\t', col.names=T, row.names=F)

res.pca <- PCA(t(res[,2:17]), graph=FALSE)
ind <- get_pca_ind(res.pca)
ind_coord <- ind$coord
ind_coord <- data.frame(ind_coord)
ind_coord$group <- c("o10d","o10d","o10d","o10t","o10t","o10t","o8d","o8d","o8d","o8t","o8t","o8t","pm","pf","pm","pf")
require("ggrepel") #add text
set.seed(42)
pdf("offspringandparent.methy.5x.pca.pdf",width=4, height=4)
ggplot(ind_coord, aes(x=Dim.1, y=Dim.2, label=rownames(ind_coord), shape=group, color=group)) +
    geom_point() +
    geom_text_repel(aes(label =rownames(ind_coord)),size = 3.5) +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_hline(yintercept = 0, linetype="dashed") +
    xlab("PCA1") + ylab("PCA2") +
    scale_color_brewer(palette="Dark2") + theme_minimal() +
    theme(legend.position="none")

dev.off()
```




