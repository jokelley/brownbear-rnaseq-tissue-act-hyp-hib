---
title: "RNAseq bears hibernation vs active"
output:
  html_notebook:
    fig_caption: yes
    toc: yes
    toc_depth: 3
    toc_float: yes
  pdf_document:
    toc: yes
    toc_depth: '3'
  html_document:
    keep_md: TRUE
---


```{r, eval=FALSE,results='hide'}
install.packages("BiocManager")
BiocManager::install(c("edgeR","DESeq2","AnnotationDbi", "impute", "GO.db", "preprocessCore"))
BiocManager::install("EDAseq")
install.packages("WGCNA")
install.packages("RColorBrewer")
install.packages("gplots")
```

# Initate the project

### Load libraries

```{r, results='hide'}
library(edgeR)
library(DESeq2)
library(WGCNA)
library(reshape2)
library(EDASeq)
library(RColorBrewer)
library(gplots)
options(stringsAsFactors = FALSE);
```

Source plot functions and ID map
```{r}
source("bear_plotting_functions.R")
gene.id = read.csv("IDmap_2019_02_20.csv", header = T)
bearData.all = read.table("sampleData.txt", header = T, row.names = 1)
```

Load counts file, remove outliers and lowly expressed genes 
```{r}

gene_all <- read.csv(gzfile("gene_count_matrix-2019-02-20.csv.gz"),header=T, row.names=1)
 
# Count number of transcripts and number of individuals
dim(gene_all)

# Subset data with CFA removed 
gene_all = gene_all[,-which(colnames(gene_all) %in% c("CFA"))]

# Remove all annotated rRNA genes
gene_all = gene_all[-which(rownames(gene_all) %in% c("gene10974","gene17084","gene26915","gene26916","gene27720","gene29258","gene29259","gene29260","gene29262","gene29263","gene29856","gene6633","gene946","gene947","rna54509","rna54511")),]
                           
dim(gene_all)

# Trim out any transcripts with no counts
total <- gene_all[rowSums(gene_all) > 0, ]

# Count number of genes and individuals
dim(total)

# Remove genes that are expressed in few individuals
keep <- rowSums(cpm(gene_all)>0.5) >= 3
length(which(keep == TRUE))

# Keep data that meet the filtering criteria
gene.counts <- gene_all[keep,]

# Count number of transcripts left after trimming
dim(gene.counts)
```


## MDS for entire dataset
```{r}
group = factor(bearData.all$shortnames)
y <- DGEList(counts=gene.counts,group=group)
y <- calcNormFactors(y)
y$samples
cols = c("blue","red","black")
pch = c(15,16,17)
for.plots = factor(c("F", "F", "L", "L", "L", "M", "M", "M", "F", "F", "F", "L", "L", "M", "M", "M", "F", "F", "F", "L", "L", "L", "M", "M", "M", "F", "F", "F", "L", "L", "L", "M", "M", "M", "F", "F", "F", "L", "L", "L", "M", "M", "M", "F", "F", "F", "L", "L", "L", "M", "M", "M"))

pdf(file = "FigureS1-All_MDS.pdf")
plotMDS(y, top = 10000, col=cols[for.plots], pch=pch[for.plots], gene.selection = "common", xlab="MDS axis 1", ylab="MDS axis 2")
legend("bottomright", as.character(c("Fat","Liver","Muscle")), col = cols[unique(for.plots)], pch = pch[unique(for.plots)])
dev.off()

#Construct MDS on top 10,000 expressed transcripts 
mds <- plotMDS(y, top=10000, gene.selection = "common", cex=0.5)
```


# Adipose

```{r}
fat.counts = total[,bearData.all$tissue == 1]
fat.set <- rowSums(cpm(fat.counts)>0.5) >=3
fat.common <- fat.counts[fat.set,]
dim(fat.common)

bearData.fat = bearData.all[bearData.all$tissue == 1,c(2,4,6)]

# Convert season column to factors because this is required for DESeq dataset design variable for transformation
bearData.fat$season = factor(bearData.fat$season)

# Set up DEseq data from count matrix, use drainage as design 
ddConsensus = DESeqDataSetFromMatrix(countData = fat.common, colData = bearData.fat, design = ~season)

# Check DEseq matrix
ddConsensus

# Run transformation, blind=T to blind the transformation to the experimental design 
vsdConsensus <- varianceStabilizingTransformation(ddConsensus, blind=T)

# Pull out transformed data as matrix. rows = genes, columns = individuals
transformed.consensus= assay(vsdConsensus)

```

## Fat MDS
```{r}
group.fat = factor(c("FHi", "FHy", "FA", "FHi", "FHy", "FA", "FHi", "FHy", "FA", "FHi", "FHy", "FA", "FHi", "FHy", "FA", "FHi", "FHy"))
group.sex = factor(c(0,0,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1))

cols = c("red","blue","yellow")
#pch = c(21,22,24)
pch = c(21,22)

pdf(file = "FigureS2a-Fat_MDS.pdf")
plotMDS(transformed.consensus, col="black", pch=pch[group.sex], bg = cols[unique(group.fat)], cex = 2, top = 10000, gene.selection = "common", main="Adipose", xlab = "MDS axis 1", ylab = "MDS axis 2")
dev.off()

pdf(file = "FigureS2a-Fat_MDS-names.pdf")
plotMDS(transformed.consensus, col="black", cex = 1, top = 10000, gene.selection = "common", main="Adipose", xlab = "MDS axis 1", ylab = "MDS axis 2")
dev.off()
```

## Fat WGCNA
```{r}
# Transpose and assign to new variable for transformed expression data
datExpr = t(assay(vsdConsensus))

# Check for genes and samples with too many missing values, if allOK returns TRUE, everything has passed 
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

# Plot relationship among samples 
sampleTree = hclust(dist(datExpr), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
```

```{r}
# Set up Color matrix
datTraits = data.frame(bearData.all[,10:13])
datTraits.fat = datTraits[datTraits$colortissue == 1,]
rownames(datTraits.fat) = colnames(fat.common)
datTraits.fat = datTraits.fat[,c(2,4,3)]
```

```{r}
pdf("SampleDendogram-Fat.pdf")
# Add trait heatmap to sample tree plot 
plotDendroAndColors(sampleTree, datTraits.fat,
                    groupLabels = names(datTraits.fat), 
                    main = "Adipose - sample dendrogram and trait map")
dev.off()
```

### Choose a set of soft-thresholding powers
```{r}
powers = c(1:15)

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
pdf(file = "NetworkTopologyAnalysis-Fat.pdf", width = 12, height = 7)
par(mfrow = c(1,2));
cex1 = 0.75;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",pch=16,
     main = paste("Scale independence"));

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
#abline(h=0.80,col="blue")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", pch=16,
     main = paste("Mean connectivity"))
dev.off()
```

```{r}
net = blockwiseModules(datExpr, power = 8, networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE, 
                       #saveTOMFileBase = "fatTOM",
                       maxBlockSize = 25000,
                       verbose = 5)

table(net$colors)
```

```{r}
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

pdf("DendroTree-Fat.pdf")
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram")
dev.off()
```

```{r}
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
```

```{r}
datTraits = data.frame(bearData.all$individual, bearData.all$tissue, bearData.all$sex, bearData.all$hibernation, bearData.all$active, bearData.all$hyperphagia)
colnames(datTraits) = c("individual","tissue","sex","hibernation","active","hyperphagia")
rownames(datTraits) = rownames(bearData.all)
datTraits = datTraits[datTraits$tissue == 1,]
datTraits = datTraits[,c(1,3,4,5,6)]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
```

```{r}
pdf("TraitCorrelationMatrix-Fat.pdf")
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Adipose module-trait relationships"))
dev.off()
```

```{r}
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$hibernation);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for hibernation",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

colnames(datExpr)[moduleColors=="turquoise"]

table(moduleColors)
```

```{r}
# Need to add code here to pull out modules, etc
annot = gene.id
dim(annot)
names(annot)
probes = colnames(datExpr)
probes2annot = match(probes, annot$Gene_name)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(geneName = probes,
                       geneSymbol = annot$Annotation[probes2annot],
                       transcriptName = annot$SeqName[probes2annot],
                       LocusLinkID = annot$Description[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "Fat-WGCNA-geneInfo.csv")

# Write one csv per module
for (modCol in unique(geneInfo0$moduleColor)) {
  cat(modCol,"\t")
  tmp = subset(geneInfo0, moduleColor == modCol) #should i use drop=T)
  tmp2 <- paste("Adipose/Fat-",modCol,".csv", sep="")
  cat(dim(tmp), "\n")
  write.table(tmp, file=tmp2, sep=",", row.names=F)
}

# Write one text file with transcripts per module
for (modCol in unique(geneInfo0$moduleColor)) {
  cat(modCol,"\t")
  tmp = subset(geneInfo0, moduleColor == modCol) #should i use drop=T)
  tmp2 <- paste("Adipose/Fat-transcripts-",modCol,".txt", sep="")
  cat(dim(tmp), "\n")
  write.table(tmp$transcriptName, file=tmp2, row.names=F, col.names=F, quote=F)
}
```

```{r}
# Set up calculation of Module Eigengenes 
PCs1A    = moduleEigengenes(datExpr,  colors=moduleColors) 
ME_1A    = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a") 
MDS_1A   = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(moduleColors))

# Plot Module Eigengene values, the below is best plotted into a pdf all together
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="",sub="")
plot(MDS_1A, col= colorsA1,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTree$order
#plotMat(scale(log(t(datExpr)[ordergenes,])) , rlabels= moduleColors[ordergenes], clabels= colnames(t(datExpr)), rcols=moduleColors[ordergenes])

pdf("AdiposeModules-Fat.pdf")
for (which.module in names(table(moduleColors))){
  ME = ME_1A[, paste("ME",which.module, sep="")] 
  barplot(ME, col=which.module, main=which.module, cex.main=2, 
          ylab="eigengene expression",xlab="array sample") 
}; 
dev.off()
```

## Visualize with VISANT - Adipose

```{r}
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 8, networkType = "signed", TOMType = "signed");

moduleColors = labels2colors(net$colors)

# Select module probes
probes = colnames(datExpr)

# Select module
module = "blue";

inModule = (moduleColors==module);
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

probeToGene = data.frame(annot$Gene_name, annot$Annotation)
pSet = ifelse(is.na(probeToGene$annot.Annotation), as.character(probeToGene$annot.Gene_name), as.character(probeToGene$annot.Annotation))
probeToGeneAgain = cbind(probeToGene,pSet)
probeToGeneAgain[1:20,] 
tail(probeToGeneAgain)
probeToGene = probeToGeneAgain[,c(1,3)]
tail(probeToGene)
#probeToGene = cbind(as.character(probeToGeneAgain$combinedAnnot.gene),as.character(probeToGeneAgain$pSet))
#probeToGene[1:20,]


# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("Adipose/fat-VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)


# Because the  module is rather large, we can restrict the genes in the output to say the 30 top hub genes in the module:
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("Adipose/fat-VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)



# Select module turquoise (highly positively correlated with hibernation)
module = "turquoise";

inModule = (moduleColors==module);
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

probeToGene = data.frame(annot$Gene_name, annot$Annotation)
pSet = ifelse(is.na(probeToGene$annot.Annotation), as.character(probeToGene$annot.Gene_name), as.character(probeToGene$annot.Annotation))
probeToGeneAgain = cbind(probeToGene,pSet)
probeToGene = probeToGeneAgain[,c(1,3)]

# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("Adipose/fat-VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)


# Because the  module is rather large, we can restrict the genes in the output to say the 30 top hub genes in the module:
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("Adipose/fat-VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)

# names of top probes
modProbes[top]

write.table(modProbes[top], file = paste("Adiposefat-", module, "-top30.txt", sep=""), quote=F)

annot = gene.id
# Loop over all modules 
for (module in names(table(moduleColors))){
  moduleColors = labels2colors(net$colors)
  probes = colnames(datExpr)

  inModule = (moduleColors==module);
  modProbes = probes[inModule];

  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)

probeToGene = data.frame(annot$Gene_name, annot$Annotation)
pSet = ifelse(is.na(probeToGene$annot.Annotation), as.character(probeToGene$annot.Gene_name), as.character(probeToGene$annot.Annotation))
probeToGeneAgain = cbind(probeToGene,pSet)
probeToGene = probeToGeneAgain[,c(1,3)]

# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("Adipose/fat-VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)

# Because the  module is rather large, we can restrict the genes in the output to say the 30 top hub genes in the module:
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("Adipose/fat-VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)
}; 

```


# Liver

```{r}
liver.counts = total[,bearData.all$tissue == 2]
liver.set <- rowSums(cpm(liver.counts)>0.5) >=3
liver.common <- liver.counts[liver.set,]
dim(liver.common)

bearData.liver = bearData.all[bearData.all$tissue == 2,c(2,4,6)]

# Convert season column to factors because this is required for DESeq dataset design variable for transformation
bearData.liver$season = factor(bearData.liver$season)

# Set up DEseq data from count matrix, use drainage as design 
ddConsensus = DESeqDataSetFromMatrix(countData = liver.common, colData = bearData.liver, design = ~season)

# Check DEseq matrix
ddConsensus

# Run transformation, blind=T to blind the transformation to the experimental design 
vsdConsensus <- varianceStabilizingTransformation(ddConsensus, blind=T)

# Pull out transformed data as matrix. rows = genes, columns = individuals
transformed.consensus= assay(vsdConsensus)
```

## MDS

```{r}
group.liver = factor(c("LA", "LHi", "LHy", "LHi", "LHy", "LA", "LHi", "LHy", "LA", "LHi", "LHy", "LA", "LHi", "LHy", "LA", "LHi", "LHy"))
group.sex = factor(c(0,0,0,1,1,1,1,1,0,0,0,1,1,1,1,1,1))
cols = c("red","blue","yellow")
#pch = c(21,22,24)
pch = c(21,22)

pdf(file = "FigureS2b-Liver_MDS.pdf")
plotMDS(transformed.consensus, col="black", bg=cols[group.liver], pch=pch[group.sex], cex = 2, top = 10000, gene.selection = "common", main="Liver", xlab = "MDS axis 1", ylab = "MDS axis 2")
dev.off()

pdf(file = "FigureS2b-Liver_MDS-names.pdf")
plotMDS(transformed.consensus, col="black", cex = 1, top = 10000, gene.selection = "common", main="Liver", xlab = "MDS axis 1", ylab = "MDS axis 2")
dev.off()
```


## WGCNA 
```{r}
# Transpose and assign to new variable for transformed expression data
datExpr = t(assay(vsdConsensus))

# Check for genes and samples with too many missing values, if allOK returns TRUE, everything has passed 
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

# Plot relationship among samples 
sampleTree = hclust(dist(datExpr), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
```

```{r}
# Set up Color matrix
datTraits = data.frame(bearData.all[,10:13])
datTraits.liver = datTraits[datTraits$colortissue == 2,]
rownames(datTraits.liver) = colnames(liver.common)
datTraits.liver = datTraits.liver[,c(2,4,3)]
```


```{r}
pdf("SampleDendogram-liver.pdf")
# Add trait heatmap to sample tree plot 
plotDendroAndColors(sampleTree, datTraits.liver,
                    groupLabels = names(datTraits.liver), 
                    main = "Liver - sample dendrogram and trait map")
dev.off()
```

```{r}
# Choose a set of soft-thresholding powers
powers = c(1:15)

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
pdf(file = "NetworkTopologyAnalysis-liver.pdf", width = 12, height = 7) 
par(mfrow = c(1,2));
cex1 = 0.75;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",pch=16,
     main = paste("Scale independence"));

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", main = paste("Mean connectivity"))
dev.off()
```



```{r}
net = blockwiseModules(datExpr, power = 5, networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, 
                       saveTOMFileBase = "Liver/liverTOM",
                       maxBlockSize = 25000,
                       verbose = 5)
table(net$colors)
```


```{r}
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

pdf("DendroTree-liver.pdf")
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram")
dev.off()
```

```{r}
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
```

```{r}
datTraits = data.frame(bearData.all$individual, bearData.all$tissue, bearData.all$sex, bearData.all$hibernation, bearData.all$active, bearData.all$hyperphagia)
colnames(datTraits) = c("individual","tissue","sex","hibernation","active","hyperphagia")
rownames(datTraits) = rownames(bearData.all)
datTraits = datTraits[datTraits$tissue == 2,]
datTraits = datTraits[,c(1,3,4,5,6)]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


pdf("TraitCorrelationMatrix-liver.pdf")
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Liver module-trait relationships"))
dev.off()
```


```{r}
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$hibernation);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for hibernation",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

colnames(datExpr)[moduleColors=="turquoise"]

table(moduleColors)
```

```{r}

probes = colnames(datExpr)
probes2annot = match(probes, annot$Gene_name)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
# Create the starting data frame
geneInfo0 = data.frame(geneName = probes,
                       geneSymbol = annot$Annotation[probes2annot],
                       transcriptName = annot$SeqName[probes2annot],
                       LocusLinkID = annot$Description[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "liver-WGCNA-geneInfo.csv")

# Write one csv per module
for (modCol in unique(geneInfo0$moduleColor)) {
  cat(modCol,"\t")
  tmp = subset(geneInfo0, moduleColor == modCol) #should i use drop=T)
  tmp2 <- paste("Liver/liver-",modCol,".csv", sep="")
  cat(dim(tmp), "\n")
  write.table(tmp, file=tmp2, sep=",", row.names=F)
}

# Write one text file with transcripts per module
for (modCol in unique(geneInfo0$moduleColor)) {
  cat(modCol,"\t")
  tmp = subset(geneInfo0, moduleColor == modCol) #should i use drop=T)
  tmp2 <- paste("Liver/Liver-transcripts-",modCol,".txt", sep="")
  cat(dim(tmp), "\n")
  write.table(tmp$transcriptName, file=tmp2, row.names=F, col.names=F, quote=F)
}


# Set up calculation of Module Eigengenes 
PCs1A    = moduleEigengenes(datExpr,  colors=moduleColors) 
ME_1A    = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a") 
MDS_1A   = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(moduleColors))

# Plot Module Eigengene values, the below is best plotted into a pdf all together
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="",sub="")
plot(MDS_1A, col= colorsA1,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTree$order
#plotMat(scale(log(t(datExpr)[ordergenes,])) , rlabels= moduleColors[ordergenes], clabels= colnames(t(datExpr)), rcols=moduleColors[ordergenes])

pdf("Liver/Modules-liver.pdf")
for (which.module in names(table(moduleColors))){
  ME = ME_1A[, paste("ME",which.module, sep="")] 
  barplot(ME, col=which.module, main="", cex.main=2, 
          ylab="eigengene expression",xlab="array sample") 
}; 
dev.off()
```

## Visualize with VISTANT - Liver

```{r}
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 5, networkType = "signed", TOMType = "signed");

# Select module
module = "blue";

moduleColors = labels2colors(net$colors)

# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

probeToGene = data.frame(annot$Gene_name, annot$Annotation)
pSet = ifelse(is.na(probeToGene$annot.Annotation), as.character(probeToGene$annot.Gene_name), as.character(probeToGene$annot.Annotation))
probeToGeneAgain = cbind(probeToGene,pSet)
probeToGeneAgain[1:20,] 
tail(probeToGeneAgain)
probeToGene = probeToGeneAgain[,c(1,3)]
tail(probeToGene)

# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("Liver/liver-VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)


# Because the  module is rather large, we can restrict the genes in the output to say the 30 top hub genes in the module:
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("Liver/liver-VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)


# names of top probes
modProbes[top]

write.table(modProbes[top], file = paste("Liver/liver-", module, "-top30.txt", sep=""), quote=F)


# Loop over all modules 
for (module in names(table(moduleColors))){
  moduleColors = labels2colors(net$colors)
  probes = colnames(datExpr)
  annot = gene.id

  inModule = (moduleColors==module);
  modProbes = probes[inModule];

  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)

probeToGene = data.frame(annot$Gene_name, annot$Annotation)
pSet = ifelse(is.na(probeToGene$annot.Annotation), as.character(probeToGene$annot.Gene_name), as.character(probeToGene$annot.Annotation))
probeToGeneAgain = cbind(probeToGene,pSet)
probeToGene = probeToGeneAgain[,c(1,3)]

# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("Liver/liver-VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)

# Because the  module is rather large, we can restrict the genes in the output to say the 30 top hub genes in the module:
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("Liver/liver-VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)
}; 

```


# Muscle

```{r}
muscle.counts = total[,bearData.all$tissue == 3]
muscle.set <- rowSums(cpm(muscle.counts)>0.5) >=3
muscle.common <- muscle.counts[muscle.set,]
dim(muscle.common)

bearData.muscle = bearData.all[bearData.all$tissue == 3,c(2,4,6)]

# Convert season column to factors because this is required for DESeq dataset design variable for transformation
bearData.muscle$season = factor(bearData.muscle$season)

# Set up DEseq data from count matrix, use drainage as design 
ddConsensus = DESeqDataSetFromMatrix(countData = muscle.common, colData = bearData.muscle, design = ~season)

# Check DEseq matrix
ddConsensus

# Run transformation, blind=T to blind the transformation to the experimental design 
vsdConsensus <- varianceStabilizingTransformation(ddConsensus, blind=T)

# Pull out transformed data as matrix. rows = genes, columns = individuals
transformed.consensus= assay(vsdConsensus)
```


## MDS

```{r}
group.muscle = factor(c("MA", "MHi", "MHy", "MA", "MHi", "MHy", "MA", "MHi", "MHy", "MA", "MHi", "MHy", "MA", "MHi", "MHy", "MA", "MHi", "MHy"))
group.sex = factor(c(0,0,0,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1))
cols = c("red","blue","yellow")
pch = c(21,22)

pdf(file = "FigureS2c-Muscle_MDS.pdf")
plotMDS(transformed.consensus, col="black", bg=cols[group.muscle], pch=pch[group.sex], cex = 2, top = 10000, gene.selection = "common", main="Muscle", xlab = "MDS axis 1", ylab = "MDS axis 2")
#legend("topleft", as.character(c("Active","Hibernation","Hyperphagia")), fill = cols[unique(group.muscle)], bty = "n")
dev.off()

pdf(file = "FigureS2c-Muscle_MDS-names.pdf")
plotMDS(transformed.consensus, col="black", cex = 1, top = 10000, gene.selection = "common", main="Muscle", xlab = "MDS axis 1", ylab = "MDS axis 2")
dev.off()
```


## WGCNA

```{r}
# Transpose and assign to new variable for transformed expression data
datExpr = t(assay(vsdConsensus))

# Check for genes and samples with too many missing values, if allOK returns TRUE, everything has passed 
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

# Plot relationship among samples 
sampleTree = hclust(dist(datExpr), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
```

```{r}
# Set up Color matrix
datTraits = data.frame(bearData.all[,10:13])
datTraits.muscle = datTraits[datTraits$colortissue == 3,]
rownames(datTraits.muscle) = colnames(muscle.common)
datTraits.muscle = datTraits.muscle[,c(2,4,3)]
```


```{r}
pdf("SampleDendogram-muscle.pdf")
# Add trait heatmap to sample tree plot 
plotDendroAndColors(sampleTree, datTraits.muscle,
                    groupLabels = names(datTraits.muscle), 
                    main = "Muscle - sample dendrogram and trait map")
dev.off()
```

```{r}
# Choose a set of soft-thresholding powers
powers = c(1:15)

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
pdf(file = "NetworkTopologyAnalysis-muscle.pdf", width = 12, height = 7) 
par(mfrow = c(1,2));
cex1 = 0.75;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",pch=16,
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
#abline(h=0.80,col="blue")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
```

```{r}
net = blockwiseModules(datExpr, power = 8, networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE, 
                       #saveTOMFileBase = "muscleTOM",
                       maxBlockSize = 25000,
                       verbose = 5)


table(net$colors)
```

```{r}
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

pdf("DendroTree-muscle.pdf")
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram")
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
```

```{r}
# Set up Trait matrix again 
datTraits = data.frame(bearData.all$individual, bearData.all$tissue, bearData.all$sex, bearData.all$hibernation, bearData.all$active, bearData.all$hyperphagia)
colnames(datTraits) = c("individual","tissue","sex","hibernation","active","hyperphagia")
rownames(datTraits) = rownames(bearData.all)
datTraits = datTraits[datTraits$tissue == 3,]
datTraits = datTraits[,c(1,3,4,5,6)]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

pdf("TraitCorrelationMatrix-muscle.pdf")
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Muscle module-trait relationships"))
dev.off()
```

```{r}
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$hibernation);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for hibernation",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

colnames(datExpr)[moduleColors=="turquoise"]

table(moduleColors)
```

```{r}
probes = colnames(datExpr)
probes2annot = match(probes, annot$Gene_name)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(geneName = probes,
                       geneSymbol = annot$Annotation[probes2annot],
                       transcriptName = annot$SeqName[probes2annot],
                       LocusLinkID = annot$Description[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "Muscle/muscle-WGCNA-geneInfo.csv")

# Write one csv per module
for (modCol in unique(geneInfo0$moduleColor)) {
  cat(modCol,"\t")
  tmp = subset(geneInfo0, moduleColor == modCol) #should i use drop=T)
  tmp2 <- paste("Muscle/muscle-",modCol,".csv", sep="")
  cat(dim(tmp), "\n")
  write.table(tmp, file=tmp2, sep=",", row.names=F)
}


# Write one text file with transcripts per module
for (modCol in unique(geneInfo0$moduleColor)) {
  cat(modCol,"\t")
  tmp = subset(geneInfo0, moduleColor == modCol) #should i use drop=T)
  tmp2 <- paste("Muscle/muscle-transcripts-",modCol,".txt", sep="")
  cat(dim(tmp), "\n")
  write.table(tmp$transcriptName, file=tmp2, row.names=F, col.names=F, quote=F)
}


# Set up calculation of Module Eigengenes 
PCs1A    = moduleEigengenes(datExpr,  colors=moduleColors) 
ME_1A    = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a") 
MDS_1A   = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(moduleColors))

# Plot Module Eigengene values, the below is best plotted into a pdf all together
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="",sub="")
plot(MDS_1A, col= colorsA1,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTree$order
#plotMat(scale(log(t(datExpr)[ordergenes,])) , rlabels= moduleColors[ordergenes], clabels= colnames(t(datExpr)), rcols=moduleColors[ordergenes])

pdf("Muscle/Modules-muscle.pdf")
for (which.module in names(table(moduleColors))){
  ME = ME_1A[, paste("ME",which.module, sep="")] 
  barplot(ME, col=which.module, main="", cex.main=2, 
          ylab="eigengene expression",xlab="array sample") 
}; 
dev.off()

```

## Visualize with VISTANT - Muscle

```{r}
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 8, networkType = "signed", TOMType = "signed");

# Select module
module = "brown";

moduleColors = labels2colors(net$colors)

# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

probeToGene = data.frame(annot$Gene_name, annot$Annotation)
pSet = ifelse(is.na(probeToGene$annot.Annotation), as.character(probeToGene$annot.Gene_name), as.character(probeToGene$annot.Annotation))
probeToGeneAgain = cbind(probeToGene,pSet)
probeToGeneAgain[1:20,] 
tail(probeToGeneAgain)
probeToGene = probeToGeneAgain[,c(1,3)]
tail(probeToGene)

# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("Muscle/muscle-VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)


# Because the  module is rather large, we can restrict the genes in the output to say the 30 top hub genes in the module:
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("Muscle/muscle-VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)

# names of top probes
modProbes[top]

write.table(modProbes[top], file = paste("Muscle/muscle-", module, "-top30.txt", sep=""), quote=F)

# Loop over all modules 
for (module in names(table(moduleColors))){
  moduleColors = labels2colors(net$colors)
  probes = colnames(datExpr)
  annot = gene.id

  inModule = (moduleColors==module);
  modProbes = probes[inModule];

  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)

probeToGene = data.frame(annot$Gene_name, annot$Annotation)
pSet = ifelse(is.na(probeToGene$annot.Annotation), as.character(probeToGene$annot.Gene_name), as.character(probeToGene$annot.Annotation))
probeToGeneAgain = cbind(probeToGene,pSet)
probeToGene = probeToGeneAgain[,c(1,3)]

# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("Muscle/muscle-VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)

# Because the  module is rather large, we can restrict the genes in the output to say the 30 top hub genes in the module:
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("Muscle/muscle-VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = probeToGene)
}; 

```



# Differential gene expression analysis

## Create DGE list for each tissue
Filter the genes with low counts out of the data and calculate normalization factors
The calcNormFactors function normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes.
```{r}
fat.counts = total[,bearData.all$tissue == 1]
group.fat = factor(c("FHi", "FHy", "FA", "FHi", "FHy", "FA", "FHi", "FHy", "FA", "FHi", "FHy", "FA", "FHi", "FHy", "FA", "FHi", "FHy"))
y.fat <- DGEList(counts=fat.counts,group=group.fat)
keep.fat <- rowSums(cpm(y.fat)>0.5) >=3
keep.fat <- y.fat[keep.fat, , keep.lib.sizes=FALSE]
keep.fat <- calcNormFactors(keep.fat)

liver.counts = total[,bearData.all$tissue == 2]
group.liver = factor(c("LA", "LHi", "LHy", "LHi", "LHy", "LA", "LHi", "LHy", "LA", "LHi", "LHy", "LA", "LHi", "LHy", "LA", "LHi", "LHy"))
y.liver <- DGEList(counts=liver.counts,group=group.liver)
keep.liver <- rowSums(cpm(y.liver)>0.5) >=3
keep.liver <- y.liver[keep.liver, , keep.lib.sizes=FALSE]
keep.liver <- calcNormFactors(keep.liver)

muscle.counts = total[,bearData.all$tissue == 3]
group.muscle = factor(c("MA", "MHi", "MHy", "MA", "MHi", "MHy", "MA", "MHi", "MHy", "MA", "MHi", "MHy", "MA", "MHi", "MHy", "MA", "MHi", "MHy"))
y.muscle <- DGEList(counts=muscle.counts,group=group.muscle)
keep.muscle <- rowSums(cpm(y.muscle)>0.5) >=3
keep.muscle <- y.muscle[keep.muscle, , keep.lib.sizes=FALSE]
keep.muscle <- calcNormFactors(keep.muscle)
```

## Relative log expression plots and PCA from DEseq tools

```{r}
filtered <- as.matrix(gene.counts)
uq <- betweenLaneNormalization(filtered, which="median")
colors <- brewer.pal(3, "Set1")
plotRLE(uq, outline=FALSE, las=3, ylab="Relative Log Expression", cex.axis=1, cex.lab=1)
plotPCA(uq, cex=1, cex.axis=1, cex.lab=1)

z <- as.factor(group.fat)
names(z) <- colnames(matrix)
colors <- brewer.pal(3, "Set1")
colLib <- colors[z]
uq.fat = betweenLaneNormalization(as.matrix(fat.common), which="full")
plotRLE(uq.fat, col=colLib, outline=FALSE, las=3, ylab="Relative Log Expression", cex.axis=1, cex.lab=1)
plotPCA(uq.fat, col=colLib, cex=1, cex.axis=1, cex.lab=1)
```

```{r}
z <- as.factor(group.liver)
names(z) <- colnames(matrix)
colors <- brewer.pal(3, "Set1")
colLib <- colors[z]
uq.liver = betweenLaneNormalization(as.matrix(liver.common), which="median")
plotRLE(uq.liver, col=colLib, outline=FALSE, las=3, ylab="Relative Log Expression", cex.axis=1, cex.lab=1)
plotPCA(uq.liver, col=colLib, cex=1, cex.axis=1, cex.lab=1)

z <- as.factor(group.muscle)
names(z) <- colnames(matrix)
colors <- brewer.pal(3, "Set1")
colLib <- colors[z]
uq.muscle = betweenLaneNormalization(as.matrix(muscle.common), which="median")
plotRLE(uq.muscle, col=colLib, outline=FALSE, las=3, ylab="Relative Log Expression", cex.axis=1, cex.lab=1)
plotPCA(uq.muscle, col=colLib, cex=1, cex.axis=1, cex.lab=1)
```

## Set up design matrices for GLM 

```{r}
#Design Matrices with bear as a blocking factor
Bears.fat <- factor(c("Cooke","Cooke","Frank","Frank","Frank","John","John","John","Oakley","Oakley","Oakley","Pacino","Pacino","Pacino","Roan","Roan","Roan"))
Season.fat <- factor(c("Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia"))
fat.design <- model.matrix(~0+Season.fat+Bears.fat, data = keep.fat$samples)
fat.design

Bears.liver <- factor(c("Cooke","Cooke","Cooke","Frank","Frank","John","John","John","Oakley","Oakley","Oakley","Pacino","Pacino","Pacino","Roan","Roan","Roan"))
Season.liver <- factor(c("Active","Hibernation","Hyperphagia","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia"))
liver.design <- model.matrix(~0+Season.liver+Bears.liver, data = keep.liver$samples)
liver.design

Bears.muscle <- factor(c("Cooke","Cooke","Cooke","Frank","Frank","Frank","John","John","John","Oakley","Oakley","Oakley","Pacino","Pacino","Pacino","Roan","Roan","Roan"))
Season.muscle <- factor(c("Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia","Active","Hibernation","Hyperphagia"))
muscle.design <- model.matrix(~0+Season.muscle+Bears.muscle, data = keep.muscle$samples)
muscle.design
```

```{r}
#Estimate Dispersions
keep.fat <- estimateDisp(keep.fat, fat.design, robust=TRUE)
keep.fat$common.dispersion

keep.liver <- estimateDisp(keep.liver, liver.design, robust=TRUE)
keep.liver$common.dispersion

keep.muscle <- estimateDisp(keep.muscle, muscle.design, robust=TRUE)
keep.muscle$common.dispersion

#Plot the BCV for each tissue
plotBCV(keep.fat)
plotBCV(keep.liver)
plotBCV(keep.muscle)
```

#Fit the data
```{r}
fat.fit <- glmQLFit(keep.fat, fat.design, robust=TRUE)
head(fat.fit$coefficients)
plotQLDisp(fat.fit)

liver.fit <- glmQLFit(keep.liver, liver.design, robust=TRUE)
head(liver.fit$coefficients)
plotQLDisp(liver.fit)

muscle.fit <- glmQLFit(keep.muscle, muscle.design, robust=TRUE)
head(muscle.fit$coefficients)
plotQLDisp(muscle.fit)
```

## Hibernation vs non-hibernation comparison 

```{r}
#F-test
#Find numer of genes significantly up or down regulated in each tissue at a given P-value (actually FDR)
qlf.fat <- glmQLFTest(fat.fit, contrast=c(-0.5,1,-0.5,0,0,0,0,0))
topTags(qlf.fat)
summary(dt_fat <- decideTestsDGE(qlf.fat, p.value = 0.01))
summary(dt_fat <- decideTestsDGE(qlf.fat, p.value = 0.05))

qlf.liver <- glmQLFTest(liver.fit, contrast=c(-0.5,1,-0.5,0,0,0,0,0))
topTags(qlf.liver)
summary(dt_liver <- decideTestsDGE(qlf.liver, p.value = 0.01))
summary(dt_liver <- decideTestsDGE(qlf.liver, p.value = 0.05))

qlf.muscle <- glmQLFTest(muscle.fit, contrast=c(-0.5,1,-0.5,0,0,0,0,0))
topTags(qlf.muscle)
summary(dt_muscle <- decideTestsDGE(qlf.muscle, p.value = 0.01))
summary(dt_muscle <- decideTestsDGE(qlf.muscle, p.value = 0.05))
```

### DE genes

```{r}
isDE_liver <- as.logical(dt_liver)
DE_liver <- rownames(keep.liver)[isDE_liver]
isDE_muscle <- as.logical(dt_muscle)
DE_muscle <- rownames(keep.muscle)[isDE_muscle]
isDE_fat <- as.logical(dt_fat)
DE_fat <- rownames(keep.fat)[isDE_fat]
```

### Plotsmear for each tissue

```{r}
pdf(file = "Liver_plotsmear_FedvsFast.pdf")
plotSmear(qlf.liver, de.tags = DE_liver, main = "Liver", ylim = c(-8.1,8.1))
abline(h=c(-1,1), col="blue")
dev.off()
pdf(file = "Muscle_plotsmear_FedvsFast.pdf")
plotSmear(qlf.muscle, de.tags = DE_muscle, main = "Muscle", ylim = c(-5.1,5.1))
abline(h=c(-1,1), col="blue")
dev.off()

min(qlf.fat$table$logFC)
max(qlf.fat$table$logFC)

pdf(file = "Fat_plotsmear_FedvsFast.pdf")
plotSmear(qlf.fat, de.tags = DE_fat, main = "Fat", ylim = c(-11,11))
abline(h=c(-1,1), col="blue")
dev.off()
```

### Create a list of all genes with their associated FDR

```{r}
tt_liver <- topTags(qlf.liver, n = dim(qlf.liver)[[1]])$table
tt_muscle <- topTags(qlf.muscle, n = dim(qlf.muscle)[[1]])$table
tt_fat <- topTags(qlf.fat, n = dim(qlf.fat)[[1]])$table
```
### Table of all genes and significance

```{r}
merge <- merge(tt_fat, tt_liver, all = TRUE, by.x = 0, by.y = 0)
merge2 <- merge(merge, tt_muscle, all = TRUE, by.x = 1, by.y = 0)
merge3 <- merge(merge2, annot, by.x = 1, by.y = 1, all = TRUE)
write.table(merge3[,c(1:18,29,30)], file = "gene_DE_master_list-02-20-2019.csv", quote = FALSE, sep = "\t", row.names = F)
```

### Add median counts per tissue to the final spreadsheet
```{r}
f1 = cpm(fat.counts)
l1 = cpm(liver.counts)
m1 = cpm(muscle.counts)

f1.median = as.data.frame(rowMedians(f1))
rownames(f1.median) = rownames(f1)
l1.median = as.data.frame(rowMedians(l1))
rownames(l1.median) = rownames(l1)
m1.median = as.data.frame(rowMedians(m1))
rownames(m1.median) = rownames(m1)

mg1 = merge(f1.median,l1.median, by.x = 0, by.y = 0, all = T) 
mg2 = merge(mg1,m1.median, by.x = 1, by.y = 0, all = T)
colnames(mg2)  = c("gene", "medianCPM_adipose", "medianCPM_liver", "medianCPM_muscle")
alldata = merge(merge3, mg2, by.x = 1, by.y = 1, all = T)
write.table(alldata[,c(1:18,29,30,31:33)], file = "gene_DE_master_list-02-21-2019.csv", quote = FALSE, sep = "\t", row.names = F)
```


```{r}
tt_liver.DE <- topTags(qlf.liver, n = summary(dt_liver <- decideTestsDGE(qlf.liver, p.value = 0.01))[1]+summary(dt_liver <- decideTestsDGE(qlf.liver, p.value = 0.01))[3], p.value = 0.01)
liver_UP    = rownames(tt_liver.DE[tt_liver.DE$table$logFC > 0,])
liver_DOWN  = rownames(tt_liver.DE[tt_liver.DE$table$logFC < 0,])

write.table(liver_UP, file = "Liver-UP-FastedvsFed.txt", quote=F, row.names = F, col.names = F)
write.table(liver_DOWN, file = "Liver-DOWN-FastedvsFed.txt", quote=F, row.names = F, col.names = F)

tt_fat.DE <- topTags(qlf.fat, n = summary(dt_fat <- decideTestsDGE(qlf.fat, p.value = 0.01))[1]+summary(dt_fat <- decideTestsDGE(qlf.fat, p.value = 0.01))[3], p.value = 0.01)
fat_UP    = rownames(tt_fat.DE[tt_fat.DE$table$logFC > 0,])
fat_DOWN  = rownames(tt_fat.DE[tt_fat.DE$table$logFC < 0,])

write.table(fat_UP, file = "fat-UP-FastedvsFed.txt", quote=F, row.names = F, col.names = F)
write.table(fat_DOWN, file = "fat-DOWN-FastedvsFed.txt", quote=F, row.names = F, col.names = F)

tt_muscle.DE <- topTags(qlf.muscle, n = summary(dt_muscle <- decideTestsDGE(qlf.muscle, p.value = 0.01))[1]+summary(dt_muscle <- decideTestsDGE(qlf.muscle, p.value = 0.01))[3], p.value = 0.01)
muscle_UP    = rownames(tt_muscle.DE[tt_muscle.DE$table$logFC > 0,])
muscle_DOWN  = rownames(tt_muscle.DE[tt_muscle.DE$table$logFC < 0,])

write.table(muscle_UP, file = "muscle-UP-FastedvsFed.txt", quote=F, row.names = F, col.names = F)
write.table(muscle_DOWN, file = "muscle-DOWN-FastedvsFed.txt", quote=F, row.names = F, col.names = F)
```

```{r}
tt_liver.DE <- topTags(qlf.liver, n = summary(dt_liver <- decideTestsDGE(qlf.liver, p.value = 0.05))[1]+summary(dt_liver <- decideTestsDGE(qlf.liver, p.value = 0.05))[3], p.value = 0.05)
liver_UP    = rownames(tt_liver.DE[tt_liver.DE$table$logFC > 0,])
liver_DOWN  = rownames(tt_liver.DE[tt_liver.DE$table$logFC < 0,])

liver_UP.1    = rownames(tt_liver.DE[tt_liver.DE$table$logFC > 1,])
liver_DOWN.1  = rownames(tt_liver.DE[tt_liver.DE$table$logFC < -1,])
length(liver_UP.1)
length(liver_DOWN.1)

write.table(liver_UP, file = "Liver-UP-FastedvsFed-0.05.txt", quote=F, row.names = F, col.names = F)
write.table(liver_DOWN, file = "Liver-DOWN-FastedvsFed-0.05.txt", quote=F, row.names = F, col.names = F)

tt_fat.DE <- topTags(qlf.fat, n = summary(dt_fat <- decideTestsDGE(qlf.fat, p.value = 0.05))[1]+summary(dt_fat <- decideTestsDGE(qlf.fat, p.value = 0.05))[3], p.value = 0.05)
fat_UP    = rownames(tt_fat.DE[tt_fat.DE$table$logFC > 0,])
fat_DOWN  = rownames(tt_fat.DE[tt_fat.DE$table$logFC < 0,])

fat_UP.1    = rownames(tt_fat.DE[tt_fat.DE$table$logFC > 1,])
fat_DOWN.1  = rownames(tt_fat.DE[tt_fat.DE$table$logFC < -1,])
length(fat_UP.1)
length(fat_DOWN.1)

write.table(fat_UP, file = "fat-UP-FastedvsFed-0.05.txt", quote=F, row.names = F, col.names = F)
write.table(fat_DOWN, file = "fat-DOWN-FastedvsFed-0.05.txt", quote=F, row.names = F, col.names = F)

tt_muscle.DE <- topTags(qlf.muscle, n = summary(dt_muscle <- decideTestsDGE(qlf.muscle, p.value = 0.05))[1]+summary(dt_muscle <- decideTestsDGE(qlf.muscle, p.value = 0.05))[3], p.value = 0.05)
muscle_UP    = rownames(tt_muscle.DE[tt_muscle.DE$table$logFC > 0,])
muscle_DOWN  = rownames(tt_muscle.DE[tt_muscle.DE$table$logFC < 0,])

muscle_UP.1    = rownames(tt_muscle.DE[tt_muscle.DE$table$logFC > 1,])
muscle_DOWN.1  = rownames(tt_muscle.DE[tt_muscle.DE$table$logFC < -1,])
length(muscle_UP.1)
length(muscle_DOWN.1)

write.table(muscle_UP, file = "muscle-UP-FastedvsFed-0.05.txt", quote=F, row.names = F, col.names = F)
write.table(muscle_DOWN, file = "muscle-DOWN-FastedvsFed-0.05.txt", quote=F, row.names = F, col.names = F)
```

### Venn Diagrams 

```{r}

# Another way to make a Venn
universe.UP <- unique(c(fat_UP, muscle_UP, liver_UP))
GroupA <- universe.UP %in% fat_UP
GroupB <- universe.UP %in% muscle_UP
GroupC <- universe.UP %in% liver_UP
input.df <- data.frame(Fat=GroupA, Muscle=GroupB, Liver=GroupC)
colnames(input.df)=c("Fat","Muscle", "Liver")
head(input.df)
a <- vennCounts(input.df)
pdf(file = "FastvsFed_up_venn-0.05.pdf")
vennDiagram(a)
dev.off()
shared.up <- universe.UP[which(input.df["Fat"] == T & input.df["Muscle"] == T & input.df["Liver"] == T)]
write.table(shared.up, file = "Shared-UpRegulated-FastedvsFed-0.05.txt", quote=F, row.names = F, col.names = F)

universe.DOWN <- unique(c(fat_DOWN, muscle_DOWN, liver_DOWN))
GroupA <- universe.DOWN %in% fat_DOWN
GroupB <- universe.DOWN %in% muscle_DOWN
GroupC <- universe.DOWN %in% liver_DOWN
input.df <- data.frame(Fat=GroupA, Muscle=GroupB, Liver=GroupC)
colnames(input.df)=c("Fat","Muscle", "Liver")
head(input.df)
a <- vennCounts(input.df)
pdf(file = "FastvsFed_down_venn-0.05.pdf")
vennDiagram(a)
dev.off()
shared.DOWN <- universe.DOWN[which(input.df["Fat"] == T & input.df["Muscle"] == T & input.df["Liver"] == T)]
write.table(shared.DOWN, file = "Shared-DownRegulated-FastedvsFed-0.05.txt", quote=F, row.names = F, col.names = F)
```

# Active vs Hyperphagia comparison

```{r}
qlf.acvshy.liver <- glmQLFTest(liver.fit, contrast=c(1,0,-1,0,0,0,0,0))
topTags(qlf.acvshy.liver)
dt_acvshy.liver <- decideTestsDGE(qlf.acvshy.liver, p.value = 0.01)
dt_acvshy.liver <- decideTestsDGE(qlf.acvshy.liver, p.value = 0.05)
summary(dt_acvshy.liver)
topTags(qlf.acvshy.liver)

qlf.acvshy.fat <- glmQLFTest(fat.fit, contrast=c(1,0,-1,0,0,0,0,0))
topTags(qlf.acvshy.fat)
dt_acvshy.fat <- decideTestsDGE(qlf.acvshy.fat, p.value = 0.01)
dt_acvshy.fat <- decideTestsDGE(qlf.acvshy.fat, p.value = 0.05)
summary(dt_acvshy.fat)

tt_acvshy.fat.DE <- topTags(qlf.acvshy.fat, n = summary(dt_qlf.acvshy.fat <- decideTestsDGE(qlf.acvshy.fat, p.value = 0.05))[1]+summary(dt_qlf.acvshy.fat <- decideTestsDGE(qlf.acvshy.fat, p.value = 0.05))[3], p.value = 0.05)

mergeInfo <- merge(tt_acvshy.fat.DE, annot, by.x = 0, by.y = 1)
write.table(mergeInfo[,c(1:8,19,20)], file = "gene_list-acvshy-03-06-2019.csv", quote = FALSE, sep = "\t", row.names = F)


qlf.acvshy.muscle <- glmQLFTest(muscle.fit, contrast=c(1,0,-1,0,0,0,0,0))
topTags(qlf.acvshy.muscle)
dt_acvshy.muscle <- decideTestsDGE(qlf.acvshy.muscle, p.value = 0.01)
dt_acvshy.muscle <- decideTestsDGE(qlf.acvshy.muscle, p.value = 0.05)
summary(dt_acvshy.muscle)
```


```{r}
tt_fat.acvshy.DE <- topTags(qlf.acvshy.fat, n = summary(dt_acvshy.fat)[1]+summary(dt_acvshy.fat)[3], p.value = 0.05)
tt_fat.acvshy.names = rownames(tt_fat.acvshy.DE)

acvshy_UP    = rownames(tt_fat.acvshy.DE[tt_fat.acvshy.DE$table$logFC > 0,])
acvshy_DOWN  = rownames(tt_fat.acvshy.DE[tt_fat.acvshy.DE$table$logFC < 0,])

write.table(acvshy_UP, file = "fat-UP-ActivevsHyper-0.05.txt", quote=F, row.names = F, col.names = F)
write.table(acvshy_UP, file = "fat-DOWN-ActivevsHyper-0.05.txt", quote=F, row.names = F, col.names = F)

```
