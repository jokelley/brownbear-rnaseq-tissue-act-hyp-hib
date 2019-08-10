setwd("C:/Users/ShawnTrojahn/Desktop/Bear_rnaseq_project/Final Products/Brown Bear")

library(edgeR)
library(DESeq2)
library(WGCNA)
library(reshape2)

library(RColorBrewer)

options(stringsAsFactors = FALSE);

gene.id = read.csv("IDmap_2019_02_20.csv", header = T)

gene_all <- read.csv("gene_count_readssubsampled.csv",header=T, row.names=1)

# Load dead bear file generated from reference and merge with other gene count matrix

#deadbears <- read.csv("gene_count_matrix_deadbears.csv", header = T, row.names = 1)
#deadbears=deadbears[,which(colnames(deadbears) %in% c("B01_S2_L003","B11_S19_L003","B17_S1_L003","D01_S3_L003","D10_S30_L003","D16_S18_L003"))]
#merge the tables by gene name
#gene_all=data.frame(merge(gene_all,deadbears, by="row.names", all.x = T),row.names = 1)



bearData.all = read.table("sampleData.txt", header = T, row.names = 1)

# Count number of transcripts and number of individuals
dim(gene_all)



# Subset data with CFA removed 
gene_all = gene_all[,-which(colnames(gene_all) %in% c("CFA"))]

# Remove all annotated rRNA genes
gene_all = gene_all[-which(rownames(gene_all) %in% c("gene10974","gene17084","gene26915","gene26916","gene27720","gene29258","gene29259","gene29260","gene29262","gene29263","gene29856","gene6633","gene946","gene947","rna54509","rna54511")),]

dim(gene_all)

# Trim out any transcripts with no counts
total <- gene_all[rowSums(gene_all) > 0, ]

#sorting
gene_all =total[, order(names(gene_all))]
total = total[, order(names(total))]
# Count number of genes and individuals 
dim(gene_all)


# Remove genes that are expressed in few individuals
keep <- rowSums(cpm(gene_all)>0.5) >= 3
length(which(keep == TRUE))

# Keep data that meet the filtering criteria
gene.counts <- gene_all[keep,]

# Count number of transcripts left after trimming
dim(gene.counts)

#MDS for entire dataset
group = factor(bearData.all$shortnames)
y <- DGEList(counts=gene.counts,group=group)
y <- calcNormFactors(y)
y$samples
cols = c("blue","red","black","blue","black","red")
group.year = factor(c(1,1,1,0,0,0,0,0,0,0,0,1,1,1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
pch = c(15,16,17,0,2,1)
for.plots = factor(c("z1","z2","z3","F", "F", "L", "L", "L", "M", "M", "M", "z1","z2","z3","F", "F", "F", "L", "L", "M", "M", "M", "F", "F", "F", "L", "L", "L", "M", "M", "M", "F", "F", "F", "L", "L", "L", "M", "M", "M", "F", "F", "F", "L", "L", "L", "M", "M", "M", "F", "F", "F", "L", "L", "L", "M", "M", "M"))

pdf(file = "FigureS1-All_MDS.pdf")
plotMDS(y, top = 10000, col=cols[for.plots], pch=pch[for.plots], gene.selection = "common", xlab="MDS axis 1", ylab="MDS axis 2")
legend("bottom", as.character(c("Fat","Liver","Muscle")), col = c("blue","red","black"), pch = c(15,16,17))
# draw legend with lines and point but without labels and box. x.intersp controls horizontal distance between lines
L = legend("bottom", legend = rep(NA,6), col=c("blue","red","black"), ncol=2, bty='n', x.intersp=0.25, pch=c(15,16,17,0,1,2))

# use position data of previous legend to draw legend with invisble lines and points but with labels and box. x.intersp controls distance between lines and labels
legend(x = L$rect$left, y = L$rect$top, legend = c('Fat', 'Liver','Muscle'), col=rep(NA,3), ncol=1, x.intersp = 2.5, bg = NA)


dev.off()

#Construct MDS on top 10,000 expressed transcripts 
mds <- plotMDS(y, top=10000, gene.selection = "common", cex=0.5)



# ROAN VS PACINO ALL PLOTS



# Load gene counts file generated from reference
#setwd("C:/Users/ShawnTrojahn/Desktop/Bear_rnaseq_project/Final Products/Brown Bear")
#gene_all<-read.csv("gene_count_matrix-2019-02-20.csv",header=T, row.names = 1)
#deadbears <- read.csv("gene_count_matrix_deadbears.csv", header = T, row.names = 1)
#deadbears=deadbears[,which(colnames(deadbears) %in% c("B01_S2_L003","B11_S19_L003","B17_S1_L003","D01_S3_L003","D10_S30_L003","D16_S18_L003"))]
#merge the tables by gene name
#gene_all=data.frame(merge(gene_all,deadbears, by="row.names", all.x = T),row.names = 1)


################################################################################################
#
#           1. Remove outliers and lowly expressed genes 
#
#
################################################################################################
setwd("C:/Users/ShawnTrojahn/Desktop/Bear_rnaseq_project/Final Products/Brown Bear")
gene_all <- read.csv("gene_count_readssubsampled.csv",header=T,row.names = 1)

# Subset data with CFA removed 
gene_all = gene_all[,-which(colnames(gene_all) %in% c("CFA"))]

# Trim out any transcripts with no counts
total <- gene_all[rowSums(gene_all) > 0, ]

#sorting
total =total[, order(names(gene_all))]

# Count number of genes and individuals 
dim(total)

#write.csv(gene_all, "gene.all.csv")


################################################################################################
#
#           2. WGCNA setup
#
#
################################################################################################


# Generate information for each of the possible variables
# 1 = Fat, 2 = Liver, 3 = Muscle
tissue.group = c(1,3,2,rep(1,2), rep(2,3), rep(3,3),1,3,2,rep(1,3), rep(2,2), rep(3,3),rep(1,3), rep(2,3), rep(3,3),rep(1,3), rep(2,3), rep(3,3),rep(1,3), rep(2,3), rep(3,3),rep(1,3), rep(2,3), rep(3,3))
# 0 = Active, 2 = Hibernation, 1 = Hyperphagia
season.group = c(rep(2,3),rep(c(2,1),1), rep(c(0,2,1),2),2,2,2,0,2,1,c(2,1),rep(c(0,2,1),13))



# Create dataframe of variables 
colData.all = data.frame(cbind(tissue.group, season.group))

# Add column names
colnames(colData.all) = c("tissue","season")

# Add individual ids to each row
rownames(colData.all) = colnames(total)


#########################################################################################################################
#
#         3. WGCNA: Subset data FAT
#
#
#########################################################################################################################

# Count matrix = total 
fat.counts = total[,tissue.group == 1]
fat.counts=head(fat.counts, -8)
fat.set <- rowSums(cpm(fat.counts)>0.5) >=3
fat.common <- fat.counts[fat.set,]
dim(fat.common)

colData.fat = colData.all[colData.all$tissue == 1,]

liver.counts = total[,tissue.group == 2]
liver.counts=head(liver.counts, -8)
liver.set <- rowSums(cpm(liver.counts)>0.5) >=3
liver.common <- liver.counts[liver.set,]
dim(liver.common)

colData.liver = colData.all[colData.all$tissue == 2,]

muscle.counts = total[,tissue.group == 3]
muscle.counts=head(muscle.counts, -8)
muscle.set <- rowSums(cpm(muscle.counts)>0.5) >=3
muscle.common <- muscle.counts[muscle.set,]
dim(muscle.common)

colData.muscle = colData.all[colData.all$tissue == 3,]
# Convert season column to factors because this is required for DESeq dataset design variable for transformation
colData.fat[,2] = factor(colData.fat[,2])

colData.liver[,2] = factor(colData.liver[,2])

colData.muscle[,2] = factor(colData.muscle[,2])

# Set up DEseq data from count matrix, use drainage as design 
ddConsensus.fat = DESeqDataSetFromMatrix(countData = fat.common, colData = colData.fat, design = ~season)

ddConsensus.liver = DESeqDataSetFromMatrix(countData = liver.common, colData = colData.liver, design = ~season)

ddConsensus.muscle = DESeqDataSetFromMatrix(countData = muscle.common, colData = colData.muscle, design = ~season)


# Run transformation, blind=T to blind the transformation to the experimental design 
vsdConsensus.fat <- varianceStabilizingTransformation(ddConsensus.fat, blind=T)

vsdConsensus.liver <- varianceStabilizingTransformation(ddConsensus.liver, blind=T)

vsdConsensus.muscle <- varianceStabilizingTransformation(ddConsensus.muscle, blind=T)

# Pull out transformed data as matrix. rows = genes, columns = individuals
transformed.consensus.fat= assay(vsdConsensus.fat)

transformed.consensus.liver= assay(vsdConsensus.liver)

transformed.consensus.muscle= assay(vsdConsensus.muscle)

pdf("fatdeseqtransformedtest.pdf")
opar <- par(no.readonly = T)
par(mfrow=c(3,2))

plot(transformed.consensus.fat[,1]~transformed.consensus.fat[,15], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,1]~transformed.consensus.fat[,15])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Adipose Bear 1", col.main="Black",
      xlab="Year 1 (Normalized Count)", ylab="Year 2 (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.fat[,4]~transformed.consensus.fat[,18], ann =F, xaxt ='n', yaxt='n')
fit <- lm(transformed.consensus.fat[,4]~transformed.consensus.fat[,18])
axis(side = 1, at=c(5,10,15,20))
axis(side = 2, at=c(5,10,15,20))
title(main="Adipose Bear 2", col.main="Black",
      xlab="Year 1 (Normalized Count)", ylab="Year 2 (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,1]~transformed.consensus.liver[,15], ann =F)
fit <- lm(transformed.consensus.liver[,1]~transformed.consensus.liver[,15])
title(main="Liver Bear 1", col.main="Black",
      xlab="Year 1 (Normalized Count)", ylab="Year 2 (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,5]~transformed.consensus.liver[,18], ann =F)
fit <- lm(transformed.consensus.liver[,5]~transformed.consensus.liver[,18])
title(main="Liver Bear 2", col.main="Black",
      xlab="Year 1 (Normalized Count)", ylab="Year 2 (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,16], ann =F)
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,16])
title(main="Muscle Bear 1", col.main="Black",
      xlab="Year 1 (Normalized Count)", ylab="Year 2 (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,5]~transformed.consensus.muscle[,19], ann =F)
fit <- lm(transformed.consensus.muscle[,5]~transformed.consensus.muscle[,19])
title(main="Muscle Bear 2", col.main="Black",
      xlab="Year 1 (Normalized Count)", ylab="Year 2 (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


par(opar)
dev.off()

#Roan or Pacino vs new season

pdf("Pacino year 1 vs year 2.pdf")
opar <- par(no.readonly = T)
par(mfrow=c(3,3))

plot(transformed.consensus.fat[,1]~transformed.consensus.fat[,14], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,1]~transformed.consensus.fat[,14])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Act Year 2 Fat", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Act (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.fat[,1]~transformed.consensus.fat[,15], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,1]~transformed.consensus.fat[,15])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hib Year 2 Fat", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hib (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.fat[,1]~transformed.consensus.fat[,16], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,1]~transformed.consensus.fat[,16])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hyp Year 2 Fat", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hyp (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,1]~transformed.consensus.liver[,14], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,1]~transformed.consensus.liver[,14])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Act Year 2 Liver", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Act (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,1]~transformed.consensus.liver[,15], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,1]~transformed.consensus.liver[,15])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hib Year 2 Liver", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hib (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,1]~transformed.consensus.liver[,16], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,1]~transformed.consensus.liver[,16])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hyp Year 2 Liver", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hyp (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,15], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,15])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Act Year 2 Muscle", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Act (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,16], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,16])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hib Year 2 Muscle", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hib (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,17], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,17])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hyp Year 2 Muscle", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hyp (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

par(opar)
dev.off()

pdf("Roan year 1 vs year 2.pdf")
opar <- par(no.readonly = T)
par(mfrow=c(3,3))

plot(transformed.consensus.fat[,4]~transformed.consensus.fat[,17], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,4]~transformed.consensus.fat[,17])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Act Year 2 Fat", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Act (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.fat[,4]~transformed.consensus.fat[,18], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,4]~transformed.consensus.fat[,18])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hib Year 2 Fat", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hib (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.fat[,4]~transformed.consensus.fat[,19], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,4]~transformed.consensus.fat[,19])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hyp Year 2 Fat", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hyp (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,5]~transformed.consensus.liver[,17], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,5]~transformed.consensus.liver[,17])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Act Year 2 Liver", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Act (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,5]~transformed.consensus.liver[,18], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,5]~transformed.consensus.liver[,18])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hib Year 2 Liver", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hib (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,5]~transformed.consensus.liver[,19], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,5]~transformed.consensus.liver[,19])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hyp Year 2 Liver", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hyp (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,5]~transformed.consensus.muscle[,18], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,5]~transformed.consensus.muscle[,18])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Act Year 2 Muscle", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Act (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,5]~transformed.consensus.muscle[,19], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,5]~transformed.consensus.muscle[,19])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hib Year 2 Muscle", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hib (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,5]~transformed.consensus.muscle[,20], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,5]~transformed.consensus.muscle[,20])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Hib Year 1 vs Hyp Year 2 Muscle", col.main="Black",
      xlab="Year 1 Hib (Normalized Count)", ylab="Year 2 Hyp (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

par(opar)
dev.off()

#RoanvsPacino

pdf("roanvspacino.pdf")
opar <- par(no.readonly = T)
par(mfrow=c(3,2))

plot(transformed.consensus.fat[,1]~transformed.consensus.fat[,4], ann =F)
fit <- lm(transformed.consensus.fat[,1]~transformed.consensus.fat[,4])
title(main="Fat Year 2", col.main="Black",
      xlab="Roan (Normalized Count)", ylab="Pacino (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.fat[,15]~transformed.consensus.fat[,18], ann =F)
fit <- lm(transformed.consensus.fat[,15]~transformed.consensus.fat[,18])
title(main="Fat Year 1", col.main="Black",
      xlab="Roan (Normalized Count)", ylab="Pacino (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,1]~transformed.consensus.liver[,5], ann =F)
fit <- lm(transformed.consensus.liver[,1]~transformed.consensus.liver[,5])
title(main="Liver Year 2", col.main="Black",
      xlab="Roan (Normalized Count)", ylab="Pacino (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,15]~transformed.consensus.liver[,18], ann =F)
fit <- lm(transformed.consensus.liver[,15]~transformed.consensus.liver[,18])
title(main="Liver Year 1", col.main="Black",
      xlab="Roan (Normalized Count)", ylab="Pacino (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,5], ann =F)
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,5])
title(main="Muscle Year 2", col.main="Black",
      xlab="Roan (Normalized Count)", ylab="Pacino (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,16]~transformed.consensus.muscle[,19], ann =F)
fit <- lm(transformed.consensus.muscle[,16]~transformed.consensus.muscle[,19])
title(main="Muscle Year 1", col.main="Black",
      xlab="Roan (Normalized Count)", ylab="Pacino (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


par(opar)
dev.off()

#Plotting Fat Vs Other Tissue (means I had to normalize all samples together)

# Load gene counts file generated from reference
setwd("C:/Users/ShawnTrojahn/Desktop/Bear_rnaseq_project/Final Products/Brown Bear")
#gene_all<-read.csv("gene_count_matrix-2019-02-20.csv",header=T, row.names = 1)
#deadbears <- read.csv("gene_count_matrix_deadbears.csv", header = T, row.names = 1)
#deadbears=deadbears[,which(colnames(deadbears) %in% c("B01_S2_L003","B11_S19_L003","B17_S1_L003","D01_S3_L003","D10_S30_L003","D16_S18_L003"))]
#merge the tables by gene name
#gene_all=data.frame(merge(gene_all,deadbears, by="row.names", all.x = T),row.names = 1)


################################################################################################
#
#           1. Remove outliers and lowly expressed genes 
#
#
################################################################################################

gene_all <- read.csv("gene_count_matrix_downsample.csv", header = T, row.names = 1)
# Subset data with CFA removed 
gene_all = gene_all[,-which(colnames(gene_all) %in% c("CFA"))]

# Trim out any transcripts with no counts
total <- gene_all[rowSums(gene_all) > 0, ]

#sorting
total =total[, order(names(gene_all))]
total = head(total, -8)
# Count number of genes and individuals 
dim(total)

#write.csv(gene_all, "gene.all.csv")


################################################################################################
#
#           2. WGCNA setup
#
#
################################################################################################


# Generate information for each of the possible variables
# 1 = Fat, 2 = Liver, 3 = Muscle
tissue.group = c(1,3,2,rep(1,2), rep(2,3), rep(3,3),1,3,2,rep(1,3), rep(2,2), rep(3,3),rep(1,3), rep(2,3), rep(3,3),rep(1,3), rep(2,3), rep(3,3),rep(1,3), rep(2,3), rep(3,3),rep(1,3), rep(2,3), rep(3,3))
# 0 = Active, 2 = Hibernation, 1 = Hyperphagia
season.group = c(rep(2,3),rep(c(2,1),1), rep(c(0,2,1),2),2,2,2,0,2,1,c(2,1),rep(c(0,2,1),13))



# Create dataframe of variables 
colData.all = data.frame(cbind(tissue.group, season.group))

# Add column names
colnames(colData.all) = c("tissue","season")

# Add individual ids to each row
rownames(colData.all) = colnames(total)

##################################################################################
#
#
#                 Down Sampled Read Plots
#
#
##################################################################################

pdf("Pacino year 1 vs year 2 downsampled reads.pdf")
opar <- par(no.readonly = T)
par(mfrow=c(3,3), oma=c(0,3,3,1.25))

plot(transformed.consensus.fat[,1]~transformed.consensus.fat[,15], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,1]~transformed.consensus.fat[,15])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Adipose", side=2, line=5, cex=1.5, font=2)
mtext("Bear 1", side=3, outer=T, cex=2, font=2)
mtext("Hibernation vs \nHibernation",side=3,cex=1,font=2,line=1)
title(xlab="Year 2 Hibernation", ylab="Year 1 Hibernation",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.fat[,1]~transformed.consensus.fat[,14], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,1]~transformed.consensus.fat[,14])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Hibernation vs \nActive",side=3,cex=1,font=2, line=1)
title(xlab="Year 2 Hibernation", ylab="Year 1 Active",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.fat[,1]~transformed.consensus.fat[,16], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,1]~transformed.consensus.fat[,16])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Hibernation vs \nHyperphagia",side=3,cex=1,font=2, line=1)
title(xlab="Year 2 Hibernation", ylab="Year 1 Hyperphagia",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.liver[,1]~transformed.consensus.liver[,15], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,1]~transformed.consensus.liver[,15])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Liver", side=2, line=5, cex=1.5, font=2)
title(xlab="Year 2 Hibernation", ylab="Year 1 Hibernation",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,1]~transformed.consensus.liver[,14], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,1]~transformed.consensus.liver[,14])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 2 Hibernation", ylab="Year 1 Active",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,1]~transformed.consensus.liver[,16], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,1]~transformed.consensus.liver[,16])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 2 Hibernation", ylab="Year 1 Hyperphagia",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,16], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,16])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Muscle", side=2, line=5, cex=1.5, font=2)
title(xlab="Year 2 Hibernation", ylab="Year 2 Hibernation",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,15], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,15])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 2 Hibernation", ylab="Year 1 Active",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,17], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,17])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 2 Hibernation", ylab="Year 1 Hyperphagia",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

par(opar)
dev.off()

pdf("Roan year 1 vs year 2 downsampled reads.pdf")
opar <- par(no.readonly = T)
par(mfrow=c(3,3), oma=c(0,3,3,1.25))

plot(transformed.consensus.fat[,4]~transformed.consensus.fat[,18], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,4]~transformed.consensus.fat[,18])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Adipose", side=2, line=5, cex=1.5, font=2)
mtext("Bear 2", side=3, outer=T, cex=2, font=2)
mtext("Hibernation vs \nHibernation",side=3,cex=1,font=2,line=1)
title(xlab="Year 2 Hibernation", ylab="Year 1 Hibernation",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.fat[,4]~transformed.consensus.fat[,17], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,4]~transformed.consensus.fat[,17])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Hibernation vs \nActive",side=3,cex=1,font=2, line=1)
title(xlab="Year 2 Hibernation", ylab="Year 1 Active",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.fat[,4]~transformed.consensus.fat[,19], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,4]~transformed.consensus.fat[,19])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Hibernation vs \nHyperphagia",side=3,cex=1,font=2, line=1)
title(xlab="Year 2 Hibernation", ylab="Year 1 Hyperphagia",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.liver[,5]~transformed.consensus.liver[,18], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,5]~transformed.consensus.liver[,18])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Liver", side=2, line=5, cex=1.5, font=2)
title(xlab="Year 2 Hibernation", ylab="Year 1 Hibernation",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,5]~transformed.consensus.liver[,17], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,5]~transformed.consensus.liver[,17])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 2 Hibernation", ylab="Year 1 Active",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,5]~transformed.consensus.liver[,19], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,5]~transformed.consensus.liver[,19])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 2 Hibernation", ylab="Year 1 Hyperphagia",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,5]~transformed.consensus.muscle[,19], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,19])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Muscle", side=2, line=5, cex=1.5, font=2)
title(xlab="Year 2 Hibernation", ylab="Year 2 Hibernation",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,18], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,18])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 2 Hibernation", ylab="Year 1 Active",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,20], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,1]~transformed.consensus.muscle[,20])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 2 Hibernation", ylab="Year 1 Hyperphagia",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

par(opar)
dev.off()


pdf("Roan year 1 vs year 1 downsampled reads.pdf")
opar <- par(no.readonly = T)
par(mfrow=c(3,3), oma=c(0,3,3,1.25))

plot(transformed.consensus.fat[,18]~transformed.consensus.fat[,18], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,18]~transformed.consensus.fat[,18])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Adipose", side=2, line=5, cex=1.5, font=2)
mtext("Roan", side=3, outer=T, cex=2, font=2)
mtext("Hibernation vs \nHibernation",side=3,cex=1,font=2,line=1)
title(xlab="Year 1 Hibernation", ylab="Year 1 Hibernation",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.fat[,18]~transformed.consensus.fat[,17], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,18]~transformed.consensus.fat[,17])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Hibernation vs \nActive",side=3,cex=1,font=2, line=1)
title(xlab="Year 1 Hibernation", ylab="Year 1 Active",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.fat[,18]~transformed.consensus.fat[,19], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.fat[,18]~transformed.consensus.fat[,19])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Hibernation vs \nHyperphagia",side=3,cex=1,font=2, line=1)
title(xlab="Year 1 Hibernation", ylab="Year 1 Hyperphagia",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.liver[,18]~transformed.consensus.liver[,18], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,18]~transformed.consensus.liver[,18])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Liver", side=2, line=5, cex=1.5, font=2)
title(xlab="Year 1 Hibernation", ylab="Year 1 Hibernation",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,18]~transformed.consensus.liver[,17], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,18]~transformed.consensus.liver[,17])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 1 Hibernation", ylab="Year 1 Active",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.liver[,18]~transformed.consensus.liver[,19], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.liver[,18]~transformed.consensus.liver[,19])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 1 Hibernation", ylab="Year 1 Hyperphagia",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.muscle[,19]~transformed.consensus.muscle[,19], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,19]~transformed.consensus.muscle[,19])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
mtext("Muscle", side=2, line=5, cex=1.5, font=2)
title(xlab="Year 1 Hibernation", ylab="Year 1 Hibernation",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.muscle[,19]~transformed.consensus.muscle[,18], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,19]~transformed.consensus.muscle[,18])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 1 Hibernation", ylab="Year 1 Active",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))


plot(transformed.consensus.muscle[,19]~transformed.consensus.muscle[,20], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.muscle[,19]~transformed.consensus.muscle[,20])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(xlab="Year 1 Hibernation", ylab="Year 1 Hyperphagia",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

par(opar)
dev.off()


#########################################################################################################################
#
#         3. WGCNA: Subset data FAT
#
#
#########################################################################################################################

# Count matrix = total 
total.counts = total
total.set <- rowSums(cpm(total.counts)>0.5)>=3
total.common <- total.counts[total.set,]
dim(total.common)

colData.total = colData.all



# Convert season column to factors because this is required for DESeq dataset design variable for transformation
colData.total[,2] = factor(colData.total[,2])


# Set up DEseq data from count matrix, use drainage as design 
ddConsensus.total = DESeqDataSetFromMatrix(countData = total.common, colData = colData.total, design = ~season)


# Run transformation, blind=T to blind the transformation to the experimental design 
vsdConsensus.total <- varianceStabilizingTransformation(ddConsensus.total, blind=T)


# Pull out transformed data as matrix. rows = genes, columns = individuals
transformed.consensus.total= assay(vsdConsensus.total)



#plotting tissues against each other in one season

pdf("Fat vs other tissues.pdf")
opar <- par(no.readonly = T)
par(mfrow=c(2,2))

plot(transformed.consensus.total[,1]~transformed.consensus.total[,3], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.total[,1]~transformed.consensus.total[,3])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Pacino Hib Year 1 Fat vs Hib Year 1 Liver", col.main="Black",
      xlab="Year 1 Fat (Normalized Count)", ylab="Year 1 Liver (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.total[,1]~transformed.consensus.total[,2], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.total[,1]~transformed.consensus.total[,2])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Pacino Hib Year 1 Fat vs Hib Year 1 Muscle", col.main="Black",
      xlab="Year 1 Fat (Normalized Count)", ylab="Year 1 Muscle (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.total[,12]~transformed.consensus.total[,14], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.total[,12]~transformed.consensus.total[,14])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Roan Hib Year 1 Fat vs Hib Year 1 Liver", col.main="Black",
      xlab="Year 1 Fat (Normalized Count)", ylab="Year 1 Liver (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

plot(transformed.consensus.total[,12]~transformed.consensus.total[,13], ann =F, xaxt ='n', yaxt ='n')
fit <- lm(transformed.consensus.total[,12]~transformed.consensus.total[,13])
axis(side =1, at=c(5,10,15,20))
axis(side=2, at=c(5,10,15,20))
title(main="Roan Hib Year 1 Fat vs Hib Year 1 Muscle", col.main="Black",
      xlab="Year 1 Fat (Normalized Count)", ylab="Year 1 Muscle (Normalized Count)",
      col.lab="black",cex.lab=1.25)
legend("topleft", bty="n", legend = paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))

par(opar)
dev.off()

####LIVER GENE PLOTS

setwd("C:/Users/ShawnTrojahn/Desktop/Bear_rnaseq_project/Final Products/Brown Bear")


library(edgeR)
library(reshape2)

idmap <- read.csv("IDmap_2019_02_20.csv", header = T)
gene_all<- read.csv("gene_count_matrix-2019-02-20.csv", header =T, row.names = 1)
gene_all = gene_all[,-which(colnames(gene_all) %in% c("CFA"))]
gene_all <- gene_all[which(colnames(gene_all) %in% c("CLA","CLHi","CLHy","FLHi", "FLHy","JLA",  "JLHi", "JLHy", "OLA" , "OLHi" ,"OLHy", "PLA",  "PLHi", "PLHy", "RLA",  "RLHi", "RLHy"))]

geneofinterest <- "gene11320"

gene_boxplot_liver("gene3191",gene_all,idmap)
gene_boxplot_liver("gene17136",gene_all,idmap)
gene_boxplot_liver("gene17747",gene_all,idmap)
gene_boxplot_liver("gene18770",gene_all,idmap)
gene_boxplot_liver("gene11320",gene_all,idmap)


# n is count matrix
# d is gene ID
# i is ID matrix to map gene iD to gene name
gene_boxplot_liver = function(d, n, i){
  cpm.genes = as.data.frame(cpm(n))
  
  group = factor(c("Liver Act", "Liver Hib", "Liver Hyp", "Liver Hib", "Liver Hyp", "Liver Act", "Liver Hib", "Liver Hyp", "Liver Act", "Liver Hib", "Liver Hyp", "Liver Act", "Liver Hib", "Liver Hyp", "Liver Act", "Liver Hib", "Liver Hyp"))
  #colgroup = factor(c("Fat", "Fat", "Liver", "Liver", "Liver", "Muscle", "Muscle", "Muscle", "Fat", "Fat", "Fat", "Liver", "Liver", "Muscle", "Muscle", "Muscle", "Fat", "Fat", "Fat", "Liver", "Liver", "Liver", "Muscle", "Muscle", "Muscle", "Fat", "Fat", "Fat", "Liver", "Liver", "Liver", "Muscle", "Muscle", "Muscle", "Fat", "Fat", "Fat", "Liver", "Liver", "Liver", "Muscle", "Muscle", "Muscle", "Fat", "Fat", "Fat", "Liver", "Liver", "Liver", "Muscle", "Muscle", "Muscle"))
  #num_colors <- nlevels(colgroup)
  #group_colors <- color_pallete_function(num_colors)
  
  top.gene = cpm.genes[d,]
  top.title = i[i$Gene_name == d,]$Annotation
  
  colnames(top.gene) = group
  ttgene = t(top.gene)
  tmgene = melt(as.matrix(ttgene))
  # the next line is new to reorder the plot
  tmgene$Var1 = factor(tmgene$Var1 , levels=levels(tmgene$Var1)[c(1,3,2)])
  
  a = boxplot(tmgene$value ~ tmgene$Var1,las=1, names.arg="", main = top.title, ylab = "Normalized Counts (cpm)", cex.axis = 0.5)
  stripchart(tmgene$value ~ tmgene$Var1, vertical=T,add=T,pch=20,cex=.75)
  pdf(file = paste(as.character(top.title),".pdf",sep=""))
  a = boxplot(tmgene$value ~ tmgene$Var1,las=1, names.arg="", main = top.title, ylab = "Normalized Counts (cpm)", cex.axis = 0.5)
  stripchart(tmgene$value ~ tmgene$Var1, vertical=T,add=T,pch=20,cex=.75)
  dev.off()
}




























