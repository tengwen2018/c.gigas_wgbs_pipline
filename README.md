**The pipline of DNA methylaiton data analysis of *Crassostrea gigas***

**1. Data trimming of Illumina 150bp PE reads**

```bash
for i in P8m P8f O8d1 O8d2 O8d3 O8t1 O8t2 O8t3 P10m P10f O10d1 O10d2 O10d3 O10t1 O10t2 O10t3
do
mkdir cleandata/$i && cd cleandata/$i
fastp -i ${i}_1.fastq.gz -I ${i}_2.fastq.gz -o ${i}_clean_1.fq.gz -O ${i}_clean_2.fq.gz\
 —adapter_sequence auto —detect_adapter_for_pe —unpaired1 output_um_1.fastq.gz\
 —unpaired2 output_um_2.fastq.gz —failed_out output_failed.fastq.gz —cut_front\
 —cut_front_window_size=1 —cut_front_mean_quality=20 —cut_tail —cut_tail_window_size=1\
 —cut_tail_mean_quality=20 —cut_right —cut_right_window_size=4 —cut_right_mean_quality=20\
 —length_required=36 —thread 1 --trim_front1 10 --trim_front2 10
done
```
**2. Mapping trimmed reads to reference genome and DNA methylation calling**

```bash
for i in P8m P8f O8d1 O8d2 O8d3 O8t1 O8t2 O8t3 P10m P10f O10d1 O10d2 O10d3 O10t1 O10t2 O10t3
do
mkdir bsmap/$i && cd bsmap/$i
ln -s ../../cleandata/$i/${i}_clean_{1..2}.fq.gz ./
bsmap -a ${i}_clean_1.fq.gz -b ${i}_clean_2.fq.gz -d ref.fa -o bsmap.$i.bam\
 -R -p 4 -n 1 -r 0 -v 0.1 -S 1 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA && \
samtools index bsmap.$i.bam && \
pefilter -i bsmap.$i.bam -o bsmap.$i.filter.bam && \
mcall -m bsmap.$i.filter.bam -p 6 -r ref.fa --sampleName $i > mcall.log
bedGraphToBigWig $i.G.bed ref.size $i.bw
done
```

**3. Differential methylation analysis**

```bash
# Male VS. female in parent cohort
mcomp -d 5 -r p8m.G.bed,p10m.G.bed -r p8f.G.bed,p10f.G.bed\
 -m male.G.bed -m female.G.bed\
 -p 10 -c mcomp.pmale.vs.pfemale.txt —withVariance 0

# Male VS. female in offspring cohort
mcomp -d 5 -r o10d_1.G.bed,o10d_3.G.bed,o8d_1.G.bed,o8d_3.G.bed\
 -r o10t_1.G.bed,o10t_2.G.bed,o10t_3.G.bed,o8t_1.G.bed,o8t_2.G.bed,o8t_3.G.bed,o8d_2.G.bed,o10d_2.G.bed\
 -m male.G.bed -m female.G.bed -p 10\
 -c mcomp.omale.vs.ofemale.txt —withVariance 0

mcomp -d 5 -r p8m.G.bed -r p8f.G.bed\
 -p 10 -c p8m.vs.p8f.txt —withVariance 0

mcomp -d 5 -r p10m.G.bed -r p10f.G.bed\
 -p 10 -c p10m.vs.p10f.txt —withVariance 0
```
**4. Curve plot of DNA methylation profile around gene bodies**

```bash
mcomp -r P8m.G.bed P8f.G.bed O8d1.G.bed O8d2.G.bed O8d3.G.bed O8t1.G.bed O8t2.G.bed O8t3.G.bed P10m.G.bed P10f.G.bed O10d1.G.bed O10d2.G.bed O10d3.G.bed O10t1.G.bed O10t2.G.bed O10t3.G.bed -m cgig.G.bed --withVariance 1 -p 4 --reference ref.fa
bedGraphToBigWig cgig.G.bed ref.size cgig.bw

computeMatrix scale-regions -R gene.bed -S cgig.bw -b 5000 -a 5000\
 --regionBodyLength 10000 --binSize 100 --startLabel "TSS" --endLabel "TTS"\
 --skipZeros --samplesLabel c.gigas -o matrix2_gene.gz\
 --outFileNameMatrix matrix2_gene.tab\
 --outFileSortedRegions matrix2_gene.bed

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
for i in P8m P8f O8d1 O8d2 O8d3 O8t1 O8t2 O8t3 P10m P10f O10d1 O10d2 O10d3 O10t1 O10t2 O10t3
do
mkdir bseqc2/$i && cd bseqc2/$i
bseqc2 -i bsmap.$i.filter.bam -o bseqc2.$i.txt -r ref.fa -l 140
done
```

**6. DNA methylation levels in variant genomic features**

```bash
for i in P8m P8f O8d1 O8d2 O8d3 O8t1 O8t2 O8t3 P10m P10f O10d1 O10d2 O10d3 O10t1 O10t2 O10t3
do
mkdir features/$i && cd features/$i
for j in `cat ../genomic_feature.list`
do
bedtools intersect -a ../../bsmap/$i/$i.G.bed -b ../$j.bed -wa | cut -f 1-5 > $j.G.bed && \
Rscript --vanilla summary.R $j.G.bed
done
done
```
```R
# summary.R
args = commandArgs(trailingOnly=TRUE)

df <- read.table(args[1],header=F)
df <- df[df$V5>=5,]
summary(df$V4)
```
```R
# Heat map visualization of DNA methylation levels in variant genomic fetures
library(ComplexHeatmap)
library(circlize)

df <- read.table("structure.G.txt",header=T,row.names=1)
df_scaled <- t(scale(t(df)))
col_fun = colorRamp2(c(-2, 0, 2), c("#1f78b4", "white", "#e31a1c"))
pdf("gene_structure_heatmap.pdf")
Heatmap(df_scaled, 
 col = col_fun,
 cluster_columns = T,
 cluster_rows = FALSE,
 row_order = c("g1stExon", "gExon", "gGene", "gIntron", "rLow_complexity", "gProm", "rLINE", "rRC_Helitron", "CGI_shore", "rLTR", "rDNA", "rSimple_repeat", "rSINE", "CGI"),
 column_order = c("p10m","o10d_1","o10d_3","o8d_1","o8d_3","p8m","p10f","p8f","o10d_2","o10t_1","o10t_2","o10t_3","o8t_1","o8t_2","o8t_3","o8d_2"),
 heatmap_legend_param=list(title = "row Z-score of mean mCG/CG"),
 width = ncol(df_scaled)*unit(5, "mm"), 
 height = nrow(df_scaled)*unit(5, "mm"),
 column_names_rot = 45
 )
dev.off()

```

**7. PCA analysis**

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

**8. The heatmap represents the Pearson correlation coefficients from pairwise comparisons of common CpG methylation levels among all samples**

```R
library(plyr)

o10d_1 <- read.table("o10d_1.G.bed", sep="\t", header=F)
o10d_1 <- o10d_1[o10d_1$V5>=5,]
o10d_1 <- o10d_1[c("V1","V2","V3","V4")]
colnames(o10d_1) <- c("chrom", "start", "end", "o10d_1")
o10d_1 <- mutate(o10d_1, pos=paste(chrom, start, end, sep="_"))
o10d_1 <- o10d_1[c("pos", "o10d_1")]

o10d_2 <- read.table("o10d_2.G.bed", sep="\t", header=F)
o10d_2 <- o10d_2[o10d_2$V5>=5,]
o10d_2 <- o10d_2[c("V1","V2","V3","V4")]
colnames(o10d_2) <- c("chrom", "start", "end", "o10d_2")
o10d_2 <- mutate(o10d_2, pos=paste(chrom, start, end, sep="_"))
o10d_2 <- o10d_2[c("pos", "o10d_2")]

o10d_3 <- read.table("o10d_3.G.bed", sep="\t", header=F)
o10d_3 <- o10d_3[o10d_3$V5>=5,]
o10d_3 <- o10d_3[c("V1","V2","V3","V4")]
colnames(o10d_3) <- c("chrom", "start", "end", "o10d_3")
o10d_3 <- mutate(o10d_3, pos=paste(chrom, start, end, sep="_"))
o10d_3 <- o10d_3[c("pos", "o10d_3")]

o10t_1 <- read.table("o10t_1.G.bed", sep="\t", header=F)
o10t_1 <- o10t_1[o10t_1$V5>=5,]
o10t_1 <- o10t_1[c("V1","V2","V3","V4")]
colnames(o10t_1) <- c("chrom", "start", "end", "o10t_1")
o10t_1 <- mutate(o10t_1, pos=paste(chrom, start, end, sep="_"))
o10t_1 <- o10t_1[c("pos", "o10t_1")]

o10t_2 <- read.table("o10t_2.G.bed", sep="\t", header=F)
o10t_2 <- o10t_2[o10t_2$V5>=5,]
o10t_2 <- o10t_2[c("V1","V2","V3","V4")]
colnames(o10t_2) <- c("chrom", "start", "end", "o10t_2")
o10t_2 <- mutate(o10t_2, pos=paste(chrom, start, end, sep="_"))
o10t_2 <- o10t_2[c("pos", "o10t_2")]

o10t_3 <- read.table("o10t_3.G.bed", sep="\t", header=F)
o10t_3 <- o10t_3[o10t_3$V5>=5,]
o10t_3 <- o10t_3[c("V1","V2","V3","V4")]
colnames(o10t_3) <- c("chrom", "start", "end", "o10t_3")
o10t_3 <- mutate(o10t_3, pos=paste(chrom, start, end, sep="_"))
o10t_3 <- o10t_3[c("pos", "o10t_3")]

o8d_1 <- read.table("o8d_1.G.bed", sep="\t", header=F)
o8d_1 <- o8d_1[o8d_1$V5>=5,]
o8d_1 <- o8d_1[c("V1","V2","V3","V4")]
colnames(o8d_1) <- c("chrom", "start", "end", "o8d_1")
o8d_1 <- mutate(o8d_1, pos=paste(chrom, start, end, sep="_"))
o8d_1 <- o8d_1[c("pos", "o8d_1")]

o8d_2 <- read.table("o8d_2.G.bed", sep="\t", header=F)
o8d_2 <- o8d_2[o8d_2$V5>=5,]
o8d_2 <- o8d_2[c("V1","V2","V3","V4")]
colnames(o8d_2) <- c("chrom", "start", "end", "o8d_2")
o8d_2 <- mutate(o8d_2, pos=paste(chrom, start, end, sep="_"))
o8d_2 <- o8d_2[c("pos", "o8d_2")]

o8d_3 <- read.table("o8d_3.G.bed", sep="\t", header=F)
o8d_3 <- o8d_3[o8d_3$V5>=5,]
o8d_3 <- o8d_3[c("V1","V2","V3","V4")]
colnames(o8d_3) <- c("chrom", "start", "end", "o8d_3")
o8d_3 <- mutate(o8d_3, pos=paste(chrom, start, end, sep="_"))
o8d_3 <- o8d_3[c("pos", "o8d_3")]

o8t_1 <- read.table("o8t_1.G.bed", sep="\t", header=F)
o8t_1 <- o8t_1[o8t_1$V5>=5,]
o8t_1 <- o8t_1[c("V1","V2","V3","V4")]
colnames(o8t_1) <- c("chrom", "start", "end", "o8t_1")
o8t_1 <- mutate(o8t_1, pos=paste(chrom, start, end, sep="_"))
o8t_1 <- o8t_1[c("pos", "o8t_1")]

o8t_2 <- read.table("o8t_2.G.bed", sep="\t", header=F)
o8t_2 <- o8t_2[o8t_2$V5>=5,]
o8t_2 <- o8t_2[c("V1","V2","V3","V4")]
colnames(o8t_2) <- c("chrom", "start", "end", "o8t_2")
o8t_2 <- mutate(o8t_2, pos=paste(chrom, start, end, sep="_"))
o8t_2 <- o8t_2[c("pos", "o8t_2")]

o8t_3 <- read.table("o8t_3.G.bed", sep="\t", header=F)
o8t_3 <- o8t_3[o8t_3$V5>=5,]
o8t_3 <- o8t_3[c("V1","V2","V3","V4")]
colnames(o8t_3) <- c("chrom", "start", "end", "o8t_3")
o8t_3 <- mutate(o8t_3, pos=paste(chrom, start, end, sep="_"))
o8t_3 <- o8t_3[c("pos", "o8t_3")]

p10m <- read.table("p10m.G.bed", sep="\t", header=F)
p10m <- p10m[p10m$V5>=5,]
p10m <- p10m[c("V1","V2","V3","V4")]
colnames(p10m) <- c("chrom", "start", "end", "p10m")
p10m <- mutate(p10m, pos=paste(chrom, start, end, sep="_"))
p10m <- p10m[c("pos", "p10m")]

p10f <- read.table("p10f.G.bed", sep="\t", header=F)
p10f <- p10f[p10f$V5>=5,]
p10f <- p10f[c("V1","V2","V3","V4")]
colnames(p10f) <- c("chrom", "start", "end", "p10f")
p10f <- mutate(p10f, pos=paste(chrom, start, end, sep="_"))
p10f <- p10f[c("pos", "p10f")]

p8m <- read.table("p8m.G.bed", sep="\t", header=F)
p8m <- p8m[p8m$V5>=5,]
p8m <- p8m[c("V1","V2","V3","V4")]
colnames(p8m) <- c("chrom", "start", "end", "p8m")
p8m <- mutate(p8m, pos=paste(chrom, start, end, sep="_"))
p8m <- p8m[c("pos", "p8m")]

p8f <- read.table("p8f.G.bed", sep="\t", header=F)
p8f <- p8f[p8f$V5>=5,]
p8f <- p8f[c("V1","V2","V3","V4")]
colnames(p8f) <- c("chrom", "start", "end", "p8f")
p8f <- mutate(p8f, pos=paste(chrom, start, end, sep="_"))
p8f <- p8f[c("pos", "p8f")]

res <- o10d_1
for (sample in list(o10d_2, o10d_3, o10t_1, o10t_2, o10t_3, o8d_1, o8d_2, o8d_3, o8t_1, o8t_2, o8t_3, p10m, p10f, p8m, p8f)){res <- merge(res, sample, by='pos', all=F)}
write.table(res, 'merged.G.txt', quote=F, sep='\t', col.names=T, row.names=F)
res <- res[c("o10d_1", "o10d_2", "o10d_3", "o10t_1", "o10t_2", "o10t_3", "o8d_1", "o8d_2", "o8d_3", "o8t_1", "o8t_2", "o8t_3", "p10m", "p10f", "p8m", "p8f")]
res <- as.matrix(res)
res <- cor(res)

library(ComplexHeatmap)
library(circlize)

pdf("cor_methy.pdf")
col_fun = colorRamp2(c(0.8, 0.9, 1), c("#1f78b4", "white", "#e31a1c"))
ht = Heatmap(res, name="cor",
col = col_fun,
cluster_rows = TRUE,
cluster_columns = TRUE,
row_names_gp = gpar(fontsize = 15),
column_names_gp = gpar(fontsize = 15),
cell_fun = function(j, i, x, y, width, height, fill) {
grid.text(sprintf("%.2f", res[i, j]), x, y, gp = gpar(fontsize = 10))
})
draw(ht)
dev.off()
```

**9. Upset plots**

```bash
ln -s ../../p8vsp10/dmc_M2_p8.G.bed_vs_p10.G.bed.txt ./
ln -s ../../o8vso10/dmc_M2_o8.G.bed_vs_o10.G.bed.txt  ./

grep -E "^#|strongHyper" dmc_M2_p8.G.bed_vs_p10.G.bed.txt | sed 's/^#//' > dmc_M2_p8.G.bed_vs_p10.G.bed.stronghyper
grep -E "^#|strongHypo" dmc_M2_p8.G.bed_vs_p10.G.bed.txt | sed 's/^#//' > dmc_M2_p8.G.bed_vs_p10.G.bed.stronghypo
grep -E "^#|strongHyper" dmc_M2_o8.G.bed_vs_o10.G.bed.txt | sed 's/^#//' > dmc_M2_o8.G.bed_vs_o10.G.bed.stronghyper
grep -E "^#|strongHypo" dmc_M2_o8.G.bed_vs_o10.G.bed.txt | sed 's/^#//' > dmc_M2_o8.G.bed_vs_o10.G.bed.stronghypo

Rscript --vanilla upset.R dmc_M2_p8.G.bed_vs_p10.G.bed.stronghypo dmc_M2_o8.G.bed_vs_o10.G.bed.stronghypo p8_vs_p10 o8_vs_o10 upset_hypo
Rscript --vanilla upset.R dmc_M2_p8.G.bed_vs_p10.G.bed.stronghyper dmc_M2_o8.G.bed_vs_o10.G.bed.stronghyper p8_vs_p10 o8_vs_o10 upset_hyper
```
```R
# upset.R

library(plyr)
library(UpSetR)

args = commandArgs(trailingOnly=TRUE)
assign(args[3], read.table(args[1], sep="\t", header=T))
assign(args[3], mutate(get(args[3]), pos=paste(chrom, start, end, sep="_")))
assign(args[3], data.frame(pos=get(args[3])$pos, name=rep(1,nrow(get(args[3])))))
x = get(args[3])
names(x)[2] <- args[3]

assign(args[4], read.table(args[2], sep="\t", header=T))
assign(args[4], mutate(get(args[4]), pos=paste(chrom, start, end, sep="_")))
assign(args[4], data.frame(pos=get(args[4])$pos, name=rep(1,nrow(get(args[4])))))
y = get(args[4])
names(y)[2] <- args[4]

res <- merge(x, y, by="pos", all=TRUE)
res[is.na(res)] = 0

write.table(res, paste(args[5], "2.txt", sep="."), quote=F, sep="\t", col.names=T, row.names=F)

pdf(paste(args[5], "2.pdf", sep="."),width=8,height=4)
upset(res, nsets=40, nintersects = 40, point.size = 2.5, line.size = 1, show.number="yes", mainbar.y.label = "intersection size", sets.x.label = "set size", sets.bar.color = "steelblue1",
 text.scale = c(2, 2, 2,1, 2, 2),
 queries = list(
  list(
    query = intersects,
    params = list(args[3], args[4]),
    color = "#Df5286",
    active = T
  )
 )
)
dev.off()
```

**10. Offspring inheriting DNA methylation states from both maternal and paternal sources**

**Family 8**
```bash
# Hyper methylated in P8f
awk -F'\t' 'BEGIN{OFS="\t"}{if($NF=="strongHyper")print $1,$2,$3}' dmr_M2_male.G.bed_vs_female.G.bed.txt > dmrhyper.bed
bedtools intersect -v -a dmrhyper.bed -b malevsfemale/dmrhyper.bed > dmrhyper2.bed # To exclude the impact of sex
bedtools intersect -a dmrhyper2.bed -b p8m.G.bed -wa -wb > p8m_dmrhyper2.txt
bedtools intersect -a dmrhyper2.bed -b p8f.G.bed -wa -wb > p8f_dmrhyper2.txt
for i in `seq 1 3`
do
bedtools intersect -a dmrhyper2.bed -b o8d$i.G.bed -wa -wb > o8d${i}_dmrhyper2.txt
bedtools intersect -a dmrhyper2.bed -b o8t$i.G.bed -wa -wb > o8t${i}_dmrhyper2.txt
done

# Hypo methylated in P8f
awk -F'\t' 'BEGIN{OFS="\t"}{if($NF=="strongHypo")print $1,$2,$3}' dmr_M2_male.G.bed_vs_female.G.bed.txt > dmrhypo.bed
bedtools intersect -v -a dmrhypo.bed -b malevsfemale/dmrhypo.bed > dmrhypo2.bed # To exclude the impact of sex
bedtools intersect -a dmrhypo2.bed -b p8m.G.bed -wa -wb > p8m_dmrhypo2.txt
bedtools intersect -a dmrhypo2.bed -b p8f.G.bed -wa -wb > p8f_dmrhypo2.txt
for i in `seq 1 3`
do
bedtools intersect -a dmrhypo2.bed -b o8d$i.G.bed -wa -wb > o8d${i}_dmrhypo2.txt
bedtools intersect -a dmrhypo2.bed -b o8t$i.G.bed -wa -wb > o8t${i}_dmrhypo2.txt
done

# No significant changes between P8f and P8m
awk -F'\t' 'BEGIN{OFS="\t"}{if($NF~"strongH")print $1,$2,$3}' dmr_M2_p8m.G.bed_vs_p8f.G.bed.txt > dmr.bed
bedtools intersect -v -a f8.G.txt -b dmr.bed -wa > f8.outof.dmr.txt
```
```R
library(ggplot2)

# Hyper methylated in P8f
p8m <- read.table("p8m_dmrhyper2.txt", header=F)
p8f <- read.table("p8f_dmrhyper2.txt", header=F)
o8d1 <- read.table("o8d1_dmrhyper2.txt", header=F)
o8d2 <- read.table("o8d2_dmrhyper2.txt", header=F)
o8d3 <- read.table("o8d3_dmrhyper2.txt", header=F)
o8t1 <- read.table("o8t1_dmrhyper2.txt", header=F)
o8t2 <- read.table("o8t2_dmrhyper2.txt", header=F)
o8t3 <- read.table("o8t3_dmrhyper2.txt", header=F)

df <- data.frame(
 methy=c(
  mean(p8m$V7), mean(p8f$V7), mean(o8d1$V7), mean(o8d2$V7), mean(o8d3$V7), mean(o8t1$V7), mean(o8t2$V7), mean(o8t3$V7)
 ), 
 group=c("p8m", "p8f", "o8d1", "o8d2", "o8d3", "o8t1", "o8t2", "o8t3")
)

bp <- ggplot(df, aes(x=factor(group, level=c("p8f", "o8d1", "o8d2", "o8d3","o8t1", "o8t2", "o8t3", "p8m")), y=methy, fill=group)) + 
 geom_bar(stat="identity")

pdf("f8.dmrhyper2.pdf", width=3,height=3)

bp + 
 scale_fill_manual(values=c("#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#D45400", "#84A3B3")) + 
 labs(title="Hyper-DMR",x="", y = "Methylation level") + 
 theme_minimal() + 
 ylim(0, 1) +
 theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

# Hypo methylated in P8f
p8m <- read.table("p8m_dmrhypo2.txt", header=F)
p8f <- read.table("p8f_dmrhypo2.txt", header=F)
o8d1 <- read.table("o8d1_dmrhypo2.txt", header=F)
o8d2 <- read.table("o8d2_dmrhypo2.txt", header=F)
o8d3 <- read.table("o8d3_dmrhypo2.txt", header=F)
o8t1 <- read.table("o8t1_dmrhypo2.txt", header=F)
o8t2 <- read.table("o8t2_dmrhypo2.txt", header=F)
o8t3 <- read.table("o8t3_dmrhypo2.txt", header=F)

df <- data.frame(
 methy=c(
  mean(p8m$V7), mean(p8f$V7), mean(o8d1$V7), mean(o8d2$V7), mean(o8d3$V7), mean(o8t1$V7), mean(o8t2$V7), mean(o8t3$V7)
 ), 
 group=c("p8m", "p8f", "o8d1", "o8d2", "o8d3", "o8t1", "o8t2", "o8t3")
)

bp <- ggplot(df, aes(x=factor(group, level=c("p8f", "o8d1", "o8d2", "o8d3","o8t1", "o8t2", "o8t3", "p8m")), y=methy, fill=group)) + 
 geom_bar(stat="identity")

pdf("dmrhypo2.pdf", width=3,height=3)

bp + 
 scale_fill_manual(values=c("#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#D45400", "#84A3B3")) + 
 labs(title="hypo2-DMR",x="", y = "Methylation level") + 
 theme_minimal() + 
 theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

# No significant changes between P8f and P8m
df <- read.table("f8.outof.dmr.txt", header=F)
colnames(df) <- c("chr", "start", "end", "o8d1", "o8d2", "o8d3", "o8t1", "o8t2", "o8t3", "p8m", "p8f")

res <- data.frame(
 methy=c(
  mean(df$o8d1), mean(df$o8d2), mean(df$o8d3), mean(df$o8t1), mean(df$o8t2), mean(df$o8t3), mean(df$p8m), mean(df$p8f)
 ), 
 group=c("o8d1", "o8d2", "o8d3", "o8t1", "o8t2", "o8t3", "p8m", "p8f")
)

bp <- ggplot(res, aes(x=factor(group, level=c("p8f", "o8d1", "o8d2", "o8d3","o8t1", "o8t2", "o8t3", "p8m")), y=methy, fill=group)) + 
 geom_bar(stat="identity")

pdf("f8.dmr.pdf", width=3,height=3)

bp + 
 scale_fill_manual(values=c("#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#D45400", "#84A3B3")) + 
 labs(title="Hyper-DMR",x="", y = "Methylation level") + 
 theme_minimal() + 
 ylim(0,1) +
 theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
```
**Family 10**
```bash
#Hyper methylated in P10f
bedtools intersect -v -a dmrhyper.bed -b malevsfemale/dmrhyper.bed > dmrhyper2.bed # To exclude the impact of sex
bedtools intersect -a dmrhyper2.bed -b p10m.G.bed -wa -wb > p10m_dmrhyper2.txt 
bedtools intersect -a dmrhyper2.bed -b p10f.G.bed -wa -wb > p10f_dmrhyper2.txt 
for i in `seq 1 3`
do
bedtools intersect -a dmrhyper2.bed -b o10d$i.G.bed -wa -wb > o10d${i}_dmrhyper2.txt
bedtools intersect -a dmrhyper2.bed -b o10t$i.G.bed -wa -wb > o10t${i}_dmrhyper2.txt
done

#Hypo methylated in P10f
bedtools intersect -v -a dmrhypo.bed -b malevsfemale/dmrhypo.bed > dmrhypo2.bed # To exclude the impact of sex
bedtools intersect -a dmrhypo2.bed -b p10m.G.bed -wa -wb > p10m_dmrhypo2.txt 
bedtools intersect -a dmrhypo2.bed -b p10f.G.bed -wa -wb > p10f_dmrhypo2.txt 
for i in `seq 1 3`
do
bedtools intersect -a dmrhypo2.bed -b o10d$i.G.bed -wa -wb > o10d${i}_dmrhypo2.txt
bedtools intersect -a dmrhypo2.bed -b o10t$i.G.bed -wa -wb > o10t${i}_dmrhypo2.txt
done

# No significant changes between P10f and P10m
awk -F'\t' 'BEGIN{OFS="\t"}{if($NF~"strongH")print $1,$2,$3}' dmr_M2_p10m.G.bed_vs_p10f.G.bed.txt > dmr.bed
bedtools intersect -v -a f10.G.txt -b dmr.bed -wa > f10.outof.dmr.txt
```
```R
library(ggplot2)
#Hyper methylated in P10f
p10m <- read.table("p10m_dmrhyper2.txt", header=F)
p10f <- read.table("p10f_dmrhyper2.txt", header=F)
o10d1 <- read.table("o10d1_dmrhyper2.txt", header=F)
o10d2 <- read.table("o10d2_dmrhyper2.txt", header=F)
o10d3 <- read.table("o10d3_dmrhyper2.txt", header=F)
o10t1 <- read.table("o10t1_dmrhyper2.txt", header=F)
o10t2 <- read.table("o10t2_dmrhyper2.txt", header=F)
o10t3 <- read.table("o10t3_dmrhyper2.txt", header=F)

df <- data.frame(
 methy=c(
  mean(p10m$V7), mean(p10f$V7), mean(o10d1$V7), mean(o10d2$V7), mean(o10d3$V7), mean(o10t1$V7), mean(o10t2$V7), mean(o10t3$V7)
 ), 
 group=c("p10m", "p10f", "o10d1", "o10d2", "o10d3", "o10t1", "o10t2", "o10t3")
)

bp <- ggplot(df, aes(x=factor(group, level=c("p10f", "o10d1", "o10d2", "o10d3", "o10t1", "o10t2", "o10t3", "p10m")), y=methy, fill=group)) + 
 geom_bar(stat="identity")

pdf("f10.dmrhyper2.pdf", width=3,height=3)

bp + 
 scale_fill_manual(values=c("#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#D45400", "#84A3B3")) + 
 labs(title="Hyper-DMR",x="", y = "Methylation level") + 
 theme_minimal() + 
 ylim(0, 1) +
 theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

#Hypo methylated in P10f
p10m <- read.table("p10m_dmrhypo2.txt", header=F)
p10f <- read.table("p10f_dmrhypo2.txt", header=F)
o10d1 <- read.table("o10d1_dmrhypo2.txt", header=F)
o10d2 <- read.table("o10d2_dmrhypo2.txt", header=F)
o10d3 <- read.table("o10d3_dmrhypo2.txt", header=F)
o10t1 <- read.table("o10t1_dmrhypo2.txt", header=F)
o10t2 <- read.table("o10t2_dmrhypo2.txt", header=F)
o10t3 <- read.table("o10t3_dmrhypo2.txt", header=F)

df <- data.frame(
 methy=c(
  mean(p10m$V7), mean(p10f$V7), mean(o10d1$V7), mean(o10d2$V7), mean(o10d3$V7), mean(o10t1$V7), mean(o10t2$V7), mean(o10t3$V7)
 ), 
 group=c("p10m", "p10f", "o10d1", "o10d2", "o10d3", "o10t1", "o10t2", "o10t3")
)

bp <- ggplot(df, aes(x=factor(group, level=c("p10f", "o10d1", "o10d2", "o10d3","o10t1", "o10t2", "o10t3", "p10m")), y=methy, fill=group)) + 
 geom_bar(stat="identity")

pdf("dmrhypo2.pdf", width=3,height=3)

bp + 
 scale_fill_manual(values=c("#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#D45400", "#84A3B3")) + 
 labs(title="hypo2-DMR",x="", y = "Methylation level") + 
 theme_minimal() + 
 theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

# No significant changes between P10f and P10m
df <- read.table("f10.outof.dmr.txt", header=F)
colnames(df) <- c("chr", "start", "end", "o10d1", "o10d2", "o10d3", "o10t1", "o10t2", "o10t3", "p10m", "p10f")

res <- data.frame(
 methy=c(
  mean(df$o10d1), mean(df$o10d2), mean(df$o10d3), mean(df$o10t1), mean(df$o10t2), mean(df$o10t3), mean(df$p10m), mean(df$p10f)
 ), 
 group=c("o10d1", "o10d2", "o10d3", "o10t1", "o10t2", "o10t3", "p10m", "p10f")
)

bp <- ggplot(res, aes(x=factor(group, level=c("p10f", "o10d1", "o10d2", "o10d3","o10t1", "o10t2", "o10t3", "p10m")), y=methy, fill=group)) + 
 geom_bar(stat="identity")

pdf("f10.dmr.pdf", width=3,height=3)

bp + 
 scale_fill_manual(values=c("#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#E0AA35", "#D45400", "#84A3B3")) + 
 labs(title="Hyper-DMR",x="", y = "Methylation level") + 
 theme_minimal() + 
 ylim(0,1) +
 theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
```


