#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# Continuous human uterine NK cell differentiation in response to endometrial 
# regeneration and pregnancy
#
# Link to publication:
# Strunz, B., Bister, J., Hamilton, R.S., Jönsson, H., Crona-Guterstam, Y., 
# Kvedaraite, E.,  Filipovic, I., Sleiers, N., Dumitrescu, B., Friberg, D., 
# Brännström, M., Lentini, A., Reinius, B., Cornillet, M., Willinger, T., 
# Gidlöf, S., Ivarsson, M.A., & Björkström, N.K. 
# [[<s>JOURNAL</s>]](https://) [[<s>DOI</s>]](https://doi.org/)
#
# Script available from:
# https://github.com/darogan/KI_Strunz_Bjorkstrom
#
#
# Analysis Performed by Russell S. Hamilton
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

library("tidyr")
library("dplyr")
library("methods")
library("utils")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("Matrix")
library("matrixStats")
library("useful")
library("reshape")
library("reshape2")
library("DESeq2")
library("biomaRt")
library("ggforce")
library("Cairo")
library("pheatmap")
library('RColorBrewer')
library("eulerr")





NUMCORES        <- 2
Project         <- "RSH_KI_0001"
elementTextSize <- 12

l2fc_cut  <- 1.5 #3

#lfcSE_cut <- 2
padj_cut  <- 0.05
CV2_cut   <- 1.95**2
reads_cut <- 10

topN      <- 25

baseDir <- "/Users/rhamilto/Documents/Repositories/KI_Strunz_Bjorkstrom/"
setwd(baseDir)
print(baseDir)

# Heatmap colours
# blue (low expression) over white (median) to red (high) with the following:
#  "dodgerblue4","white","firebrick"


message("+-------------------------------------------------------------------------------")
message("+ Set Up the Sample Table")
message("+-------------------------------------------------------------------------------")

sampleFiles <- list.files(baseDir, pattern='*featureCounts_counts.txt', recursive = TRUE)
sampleNames <- gsub("_R1_001_trimmed_GRCh38.star.bam_featureCounts_counts.txt", "", sampleFiles)
sampleNames <- gsub("_[ATGC]*_L00[0-9]", "", sampleNames)
sampleNames <- gsub("FeatureCounts\\/", "", sampleNames)
sampleNames <- gsub("Sample", "s", sampleNames)
sampleNames <- gsub("_S[0-9]*", "", sampleNames)

sampleID    <- gsub("_.*", "", sampleNames)
sampleNames <- gsub("^[0-9][0-9]_", "", sampleNames)
sampleNames <- paste0(sampleNames, "_", sampleID)

sampleOrigin <- gsub("^[0-9][0-9]_", "", sampleNames)
sampleOrigin <- gsub("_P.*", "", sampleOrigin)
sampleOrigin <- gsub("[0-9]*", "", sampleOrigin)
sampleOrigin <- gsub("x_", "", sampleOrigin)
sampleOrigin <- gsub("HU", "Hu", sampleOrigin)

sampleGroup <- sampleNames
sampleGroup <- gsub("HU12_P11_01",	 "CD39+KIR+", sampleGroup)
sampleGroup <- gsub("HU12_P9_02",	   "CD39-KIR-", sampleGroup)
sampleGroup <- gsub("HU12_P10_03",   "CD39-KIR+", sampleGroup)
sampleGroup <- gsub("HU13_P11_05",   "CD39+KIR+", sampleGroup)
sampleGroup <- gsub("HU13_P9_06",    "CD39-KIR-", sampleGroup)
sampleGroup <- gsub("HU13_P10_07",   "CD39-KIR+", sampleGroup)
sampleGroup <- gsub("Hu05_P11_09",   "CD39+KIR+", sampleGroup)
sampleGroup <- gsub("Hu05_P9_10",    "CD39-KIR-", sampleGroup)
sampleGroup <- gsub("Hu05_P10_11",   "CD39-KIR+", sampleGroup)
sampleGroup <- gsub("Hu7_P12_14",    "CD39+KIR+", sampleGroup)
sampleGroup <- gsub("Hu7_P13_15",    "CD39-KIR+", sampleGroup)
sampleGroup <- gsub("Hu7_P9_16",     "CD39-KIR-", sampleGroup)
sampleGroup <- gsub("^x", "", sampleGroup)
sampleGroup <- gsub("\\-", "m", sampleGroup)
sampleGroup <- gsub("\\+", "p", sampleGroup)

sampleTable              <- data.frame(sampleNames=sampleNames, fileNameDGE=sampleFiles, Origin=sampleOrigin, Group=sampleGroup)
sampleTable$origin_group <- paste0(sampleOrigin, "-", sampleGroup)
print(sampleTable)


message("+-------------------------------------------------------------------------------")
message("+ Use ensEMBL Annotations")
message("+-------------------------------------------------------------------------------")

ensembl    <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description','chromosome_name', 'gene_biotype'), mart = ensembl)  
head(ensEMBL2id)
nrow(ensEMBL2id)


message("+-------------------------------------------------------------------------------")
message("+ Read in the featureCount gene count files per sample")
message("+-------------------------------------------------------------------------------")

DESeqDataSetFromFeatureCounts <- function (sampleTable, directory = ".", design, ignoreRank = FALSE, ...) 
{
  # From https://www.biostars.org/p/277316/
  if (missing(design)) 
    stop("design is missing")
  l <- lapply(as.character(sampleTable[, 2]), function(fn) read.table(file.path(directory, fn), skip=2))
  if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1)))) 
    stop("Gene IDs (first column) differ between files.")
  tbl <- sapply(l, function(a) a$V7)
  colnames(tbl) <- sampleTable[, 1]
  rownames(tbl) <- l[[1]]$V1
  rownames(sampleTable) <- sampleTable[, 1]
  dds <- DESeqDataSetFromMatrix(countData = tbl, colData = sampleTable[, -(1:2), drop = FALSE], design = design, ignoreRank, ...)
  return(dds)
}

dds.Samples <- DESeqDataSetFromFeatureCounts(sampleTable=sampleTable, directory=baseDir, design= ~ Group )
dds.Samples <- DESeq(dds.Samples, parallel=TRUE)


#
# We want to compare the samples to another set of data at a later stage, so sizefactors have been 
# pre-calculated together with this larger set. So here we provide the sizefactors to reproduce the 
# figures in the paper. They can also be calculated just on the samples in the paper by skipping this step
#
sf.Hu.man <- c("HU12_P11_01"=1.0996484, "HU12_P9_02"=1.0192777, "HU12_P10_03"=1.3543006,
               "HU13_P11_05"=1.2475651, "HU13_P9_06"=1.3392764, "HU13_P10_07"=1.2173126,
               "Hu05_P11_09"=1.6901837, "Hu05_P9_10"=0.5744079, "Hu05_P10_11"=1.3150025,
               "Hu7_P12_14"=2.3282260,  "Hu7_P13_15"=0.4262548, "Hu7_P9_16"=1.0750926 )
sizeFactors(dds.Samples) <- sf.Hu.man


colData(dds.Samples)


normCounts                    <- counts(dds.Samples,      normalized=TRUE)
normCounts.df                 <- as.data.frame(normCounts)
normCounts.df$ensembl_gene_id <- rownames(normCounts.df)
normCounts.df.annot           <- merge(normCounts.df,ensEMBL2id, by="ensembl_gene_id" )
outDir <- paste0(baseDir, "/DifferentialGeneExpression/")
write.csv(normCounts.df.annot, file=paste0(outDir, Project, "_DESeq2_NormalisedReadCounts.ann.csv"))


message("+-------------------------------------------------------------------------------")
message("+ Perform Sex Assignment")
message("+-------------------------------------------------------------------------------")
# Samples should all be female - but lets just double check

normCounts                    <- counts(dds.Samples,      normalized=TRUE)
normCounts.df                 <- as.data.frame(normCounts)
normCounts.df$ensembl_gene_id <- rownames(normCounts.df)
normCounts.df.annot           <- merge(normCounts.df,ensEMBL2id, by="ensembl_gene_id" )
normCounts.df.annot.sex       <- normCounts.df.annot
normCounts.df.annot.sex       <- subset(normCounts.df.annot.sex, 
                                        (external_gene_name=="XIST"  | external_gene_name=="DDX3Y" | 
                                           external_gene_name=="KDM5D" | external_gene_name=="UTY" | 
                                           external_gene_name=="ZFY1"))

rownames(normCounts.df.annot.sex) <- normCounts.df.annot.sex$external_gene_name
normCounts.df.annot.sex           <- normCounts.df.annot.sex[ , !(names(normCounts.df.annot.sex) %in% 
                                                              c("ensembl_gene_id","description","entrezgene","gene_biotype"))]
normCounts.df.annot.sex.mlt <- melt(normCounts.df.annot.sex)
head(normCounts.df.annot.sex.mlt)

plt.sex <- ggplot(normCounts.df.annot.sex.mlt, aes(x=external_gene_name, y=value, fill=chromosome_name)) + 
           geom_bar(stat="identity", position=position_dodge()) +
           facet_wrap(~ variable, ncol=4) +
           scale_fill_manual(name="Sex Linked Gene Sets", values=(c("X"="pink", "Y"="blue"))) +
           ylab("Normalised Read Counts") +
           xlab("Sex Linked Genes") +
           theme(text=element_text(size=8,  family="sans"),
                 axis.text.y = element_text(size=8),
                 axis.text.x = element_text(size=8, angle = 45, hjust = 1),
                 legend.position="none")
plt.sex

sampleTable$Sex <- "F"


message("+-------------------------------------------------------------------------------")
message("+ Perform Custom PCA")
message("+-------------------------------------------------------------------------------")

vst.Samples       <- vst(dds.Samples,     blind=F)



customPCA <- function(sampleTBL, RLD, TOPNUM, model, ensEMBL2id, expFactorXneg, expFactorXpos, expFactorYneg, expFactorYpos, scaleTicks) {
  
  elementTextSize <- 12
  rv              <- rowVars(RLD)
  select          <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca             <- prcomp(t(RLD[select, ]))
  
  pc1var          <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var          <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab          <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab          <- paste0("PC2 (",as.character(pc2var),"%)")
  
  sampleTBL$sampleName <- gsub("^x", "", sampleTBL$sampleName)
  sampleTBL$Label <- sampleTBL$Group
  sampleTBL$Label <- gsub("p", "+", sampleTBL$Label)
  sampleTBL$Label <- gsub("m", "-", sampleTBL$Label)
  
  scores    <- data.frame(sampleName=sampleTBL$sampleName, pca$x, Origin=sampleTBL$Origin, Group=sampleTBL$Group, Label=sampleTBL$Label)
  PC1min    <- min(scores$PC1); PC1max <- max(scores$PC1)
  PC2min    <- min(scores$PC2); PC2max <- max(scores$PC2)
  
  plt.pca.nl.spl <- ggplot(scores, aes(x = PC1, y = PC2, color=paste0(Group,"_",Origin), label=sampleName) ) +
                    geom_mark_ellipse(aes(fill = NULL, color=paste0(Group,"_",Origin), group=paste0(Group), 
                                          label=Label), alpha=0.1, label.fontsize = 6, label.buffer = unit(1, 'mm')) +
                    geom_point(size = 3, alpha=0.75 ) + 
                    scale_color_manual(name="Group", 
                                       values=(c("CD39pKIRp_Hu"="orange","CD39pKIRm_Hu"="grey","CD39mKIRm_Hu"="grey50", "CD39mKIRp_Hu"="darkviolet",
                                                "CD39pKIRp_Hys"="orange","CD39pKIRm_Hys"="black", "CD39mKIRm_Hys"="grey50","CD39mKIRp_Hys"="darkviolet"))) +
                    xlab(pc1lab) + ylab(pc2lab) + coord_fixed() +
                    scale_x_continuous(breaks=seq(signif(PC1min*expFactorXneg,digits=1), signif(PC1max*expFactorXpos,digits=1), scaleTicks)) +
                    scale_y_continuous(breaks=seq(signif(PC2min*expFactorYneg,digits=1), signif(PC2max*expFactorYpos,digits=1), scaleTicks)) +
                    expand_limits(x = c(signif(PC1min,digits=1)*expFactorXneg, signif(PC1max,digits=1)*expFactorXpos), 
                                  y = c(signif(PC2min,digits=1)*expFactorYneg, signif(PC2max,digits=1)*expFactorYpos)) +
                    facet_grid( ~ Origin ) +
                    ggtitle(paste0("PCA (Top ", TOPNUM, " Most Variable Genes)")) +             
                    theme(text = element_text(size=elementTextSize), legend.position="none") 
  
  loadings                 <- as.data.frame(pca$rotation)
  loadings$ensembl_gene_id <- rownames(loadings)
  loadings                 <- merge(loadings, ensEMBL2id, by="ensembl_gene_id")
  
  pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
  pca.1.25      <- pca.1[c(1:25),]
  pca.1.max     <- max(pca.1.25$PC1)
  pca.1.min     <- min(pca.1.25$PC1)
  
  pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
  pca.2.25      <- pca.2[c(1:25),]
  pca.2.max     <- max(pca.2.25$PC2)
  pca.2.min     <- min(pca.2.25$PC2)
  
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(pca.1.25$external_gene_name,levels=unique(pca.1.25$external_gene_name)), y=PC1)) + 
                   geom_point(size=1.5) + ylim(min(pca.1.min,pca.2.min),max(pca.1.max,pca.2.max)) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(pca.2.25$external_gene_name,levels=unique(pca.2.25$external_gene_name)), y=PC2)) + 
                   geom_point(size=1.5) + ylim(min(pca.1.min,pca.2.min),max(pca.1.max,pca.2.max)) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  return(list(plt.pca.nl.spl, pca.1.25.plot, pca.2.25.plot) )
}



pca.plt.500     <- customPCA(sampleTable, assay(vst.Samples), 500,"vsd.500", ensEMBL2id, 
                                expFactorXneg=1.2, expFactorXpos=1.4, expFactorYneg=1.2, expFactorYpos=1.1, scaleTicks=10)
pc_ex           <- plot_grid(pca.plt.500[[2]], pca.plt.500[[3]], nrow=2)
pdf(paste0(baseDir, "/QC/", Project, "_Fig.PCA.Final.pdf"),width=8,height=5)
par(bg=NA)
plot_grid(pca.plt.500[[1]], pc_ex, ncol=2, labels=c("A","B"), rel_widths = c(1.25,1), align = "h", axis = "bt")
dev.off()


message("+-------------------------------------------------------------------------------")
message("+ Exploring DEG analysis filtering ")
message("+-------------------------------------------------------------------------------")


resultsNames(dds.Samples)
as.data.frame(colData(dds.Samples))

plotDispEsts(dds.Samples)


res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm  <- lfcShrink(dds=dds.Samples, contrast=c("Group", "CD39pKIRp", "CD39mKIRm"),  type="normal", parallel=TRUE)
res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp  <- lfcShrink(dds=dds.Samples, contrast=c("Group", "CD39pKIRp", "CD39mKIRp"),  type="normal", parallel=TRUE)
res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm  <- lfcShrink(dds=dds.Samples, contrast=c("Group", "CD39mKIRp", "CD39mKIRm"),  type="normal", parallel=TRUE)


message("+-------------------------------------------------------------------------------")
message("+ Annotate RESULTS using ensEMBL                                                ")
message("+-------------------------------------------------------------------------------")

functionGetSigDEG <- function(result, Project, ensEMBL2id, TITLE, l2fc_cut, padj_cut, outDir) {
  # Filter by adjusted p-value
  result.df     <- as.data.frame(result)
  result.df.sub <- as.data.frame(subset(result, padj <= padj_cut  & abs(log2FoldChange) >= l2fc_cut))  
  # Print out the number of significant hits
  message(paste0("+   ", TITLE, "        padj    = ", nrow(result.df.sub)))
  # Annotate Genes
  result.df$ensembl_gene_id <- rownames(result.df)
  # Annotate from the live ensEMBL table generated above
  result.ann <- merge(result.df,ensEMBL2id, by="ensembl_gene_id")
  # Tidy up the description field
  result.ann$description            <- gsub("..Source.*", "", result.ann$description)
  # Order the results object
  result.ann <- result.ann[order(result.ann$log2FoldChange,decreasing=TRUE),]
  # Write the results to file
  write.csv(result.ann[order(abs(result.ann$log2FoldChange),decreasing=TRUE),], 
            file=paste0(outDir, Project, "_DESeq2_DEGs_padj_", padj_cut, '_l2fc_', l2fc_cut, "_", TITLE, '.ann.csv'))
  return(result.ann)
}


outDir <- paste0(baseDir, "/DifferentialGeneExpression/")

res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.ann  <- functionGetSigDEG(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm,  Project, ensEMBL2id, "Hu_CD39pKIRp_vs_Hu_CD39mKIRm",  l2fc_cut, padj_cut, outDir)
res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.ann  <- functionGetSigDEG(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp,  Project, ensEMBL2id, "Hu_CD39pKIRp_vs_Hu_CD39mKIRp",  l2fc_cut, padj_cut, outDir)
res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.ann  <- functionGetSigDEG(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm,  Project, ensEMBL2id, "Hu_CD39mKIRp_vs_Hu_CD39mKIRm",  l2fc_cut, padj_cut, outDir)




message("+-------------------------------------------------------------------------------")
message("+ Exlore Thresholds                                                             ")
message("+-------------------------------------------------------------------------------")


functionApplyCV2Filtering <- function(results.ann, dds, ORIGIN1, ORIGIN2, GROUP1, GROUP2){

t.results                         <- results.ann
t.counts                          <- as.data.frame(counts(dds,normalized=TRUE))
t.counts$ensembl_gene_id          <- rownames(t.counts)
t.table                           <- merge(t.results, t.counts, by="ensembl_gene_id")
t.table.undup                     <- t.table[!duplicated(t.table$ensembl_gene_id),]
rownames(t.table.undup)           <- t.table.undup$ensembl_gene_id

t.table.undup.filt                <- t.table.undup[, c(1,2,3,4,7,12:23)]
t.table.undup.filt.m              <- melt(t.table.undup.filt, id=c("baseMean","log2FoldChange","lfcSE","padj", "ensembl_gene_id"))
colnames(t.table.undup.filt.m)[6] <- "sampleNames"
t.table.undup.filt.m.ann          <- merge(t.table.undup.filt.m, sampleTable, by="sampleNames")

print( colnames(t.table.undup.filt.m.ann))

t.table.undup.filt.m.ann          <- t.table.undup.filt.m.ann[, c(1:7,9,10)]
t.table.undup.filt.m.ann.agg      <- do.call(data.frame, aggregate(t.table.undup.filt.m.ann$value, 
                                             by = list(t.table.undup.filt.m.ann$ensembl_gene_id, 
                                                       t.table.undup.filt.m.ann$Origin,
                                                       t.table.undup.filt.m.ann$Group),
                                             FUN=function(x) c(mean=mean(x), sd=sd(x))) )

colnames(t.table.undup.filt.m.ann.agg) <- c("ensembl_gene_id", "Origin", "Group", "Mean", "SD")
t.table.undup.filt.m.ann.agg.mrg       <- merge(t.table.undup.filt.m.ann, t.table.undup.filt.m.ann.agg, 
                                                by=c("ensembl_gene_id", "Origin", "Group"))

test.data               <- subset(t.table.undup.filt.m.ann.agg.mrg, 
                                 (Origin==ORIGIN1 | Origin==ORIGIN2) & (Group==GROUP1 | Group==GROUP2))
test.data.u             <- test.data[,c(1,2,3,6,7,8,10,11)]
test.data.u             <- test.data.u[!duplicated(test.data.u),]
test.data.u$CV2         <- test.data.u$SD/test.data.u$Mean
test.data.u             <- merge(test.data.u,ensEMBL2id, by="ensembl_gene_id")
test.data.u$description <- gsub("..Source.*", "", test.data.u$description)

#test.data.x             <- subset(t.table.undup.filt.m.ann.agg.mrg, Origin=="Hu" & (Group==GROUP1 | Group==GROUP2))
test.data.x             <- subset(t.table.undup.filt.m.ann.agg.mrg, ((Origin==ORIGIN1 & Group==GROUP1) | (Origin==ORIGIN2 &Group==GROUP2)))

test.data.ux            <- test.data.x[,c(1,2,3,6,7,8,10,11)]
test.data.ux            <- test.data.ux[!duplicated(test.data.ux),]
test.data.ux$CV2        <- test.data.ux$SD/test.data.ux$Mean

test.data.ux.lab        <- merge(test.data.ux,ensEMBL2id, by="ensembl_gene_id")
test.data.ux.lab.u      <- test.data.ux.lab[order(test.data.ux$log2FoldChange, decreasing=T),c(1,4,6,10)]
test.data.ux.lab.u      <- test.data.ux.lab.u[!duplicated(test.data.ux.lab.u),]
test.data.ux.lab.d      <- test.data.ux.lab[order(test.data.ux$log2FoldChange, decreasing=F),c(1,4,6,10)]
test.data.ux.lab.d      <- test.data.ux.lab.d[!duplicated(test.data.ux.lab.d),]

return(list(test.data.u, test.data.ux, test.data.ux.lab.u, test.data.ux.lab.d))
}

res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_filt  <- functionApplyCV2Filtering(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.ann,  dds.Samples, 
                                                                            ORIGIN1="Hu", ORIGIN2="Hu",  
                                                                            GROUP1="CD39pKIRp", GROUP2="CD39mKIRm")
res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_filt  <- functionApplyCV2Filtering(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.ann,  dds.Samples, 
                                                                            ORIGIN1="Hu", ORIGIN2="Hu",  
                                                                            GROUP1="CD39pKIRp", GROUP2="CD39mKIRp")
res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_filt  <- functionApplyCV2Filtering(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.ann,  dds.Samples, 
                                                                            ORIGIN1="Hu", ORIGIN2="Hu",  
                                                                            GROUP1="CD39mKIRp", GROUP2="CD39mKIRm")




head(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_filt[[2]])  
head(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_filt[[2]]) 
head(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_filt[[2]])  




#
# FUNCTION: Produce a final DGE table for the comparison
#
functionMakeCV2FilteredTable <- function(xx.filt, Project, ensEMBL2id, ORIGIN1, ORIGIN2, GROUP1, GROUP2, padj_cut, l2fc_cut, CV2_cut){

 # xx.filt[is.na(xx.filt)]   <- 0

  xx.filt$GROUP1_CV2  <- with(xx.filt, ifelse(Group==GROUP1, CV2, NA))
  xx.filt$GROUP2_CV2  <- with(xx.filt, ifelse(Group==GROUP2, CV2, NA))
  xx.filt$GROUP1_Mean <- with(xx.filt, ifelse(Group==GROUP1, Mean, NA))
  xx.filt$GROUP2_Mean <- with(xx.filt, ifelse(Group==GROUP2, Mean, NA))
  xx.filt             <- xx.filt[,c("ensembl_gene_id", "Origin", "Group", "log2FoldChange", "lfcSE", "padj", "GROUP1_CV2", "GROUP2_CV2", "GROUP1_Mean", "GROUP2_Mean" )]
  colnames(xx.filt)   <- c("ensembl_gene_id", "Origin", "Group", "log2FoldChange", "lfcSE", "padj", paste0(GROUP1,"_CV2"), paste0(GROUP2,"_CV2"), paste0(GROUP1,"_Mean"), paste0(GROUP2,"_Mean") )
  
  xx.1                    <- subset(xx.filt, Origin ==ORIGIN1 & Group==GROUP1)
  xx.2                    <- subset(xx.filt, Origin ==ORIGIN2 & Group==GROUP2)
  xxx                     <- merge(xx.1, xx.2, by=c("ensembl_gene_id", "log2FoldChange","lfcSE","padj"))
  xx.filt.ann             <- merge(xxx,ensEMBL2id, by="ensembl_gene_id")
  xx.filt.ann$description <- gsub("..Source.*", "", xx.filt.ann$description)
  xx.filt.ann             <- xx.filt.ann[,c(1:4,5,6,7,9,11,12,14,16,17,18,19)]
  
  write.csv(xx.filt.ann[order(abs(xx.filt.ann$log2FoldChange),decreasing=TRUE),], 
            file=paste0("DifferentialGeneExpression/", Project, "_DESeq2_DEGs_CV2Filtered_padj_", padj_cut, '_l2fc_', l2fc_cut, "_CV2_", CV2_cut, ".", ORIGIN1, "_", GROUP1, "_vs_", ORIGIN2, "_", GROUP2, '.ann.csv'))
  
  return(xx.filt.ann)
}

res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_filt.tab  <- functionMakeCV2FilteredTable(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_filt[[2]], 
                                                                                   Project, ensEMBL2id, 
                                                                                   ORIGIN1="Hu", ORIGIN2="Hu", GROUP1="CD39pKIRp", GROUP2="CD39mKIRm",
                                                                                   padj_cut, l2fc_cut, CV2_cut)
res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_filt.tab  <- functionMakeCV2FilteredTable(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_filt[[2]], 
                                                                                   Project, ensEMBL2id, 
                                                                                   ORIGIN1="Hu", ORIGIN2="Hu", GROUP1="CD39pKIRp", GROUP2="CD39mKIRp",
                                                                                   padj_cut, l2fc_cut, CV2_cut)
res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_filt.tab  <- functionMakeCV2FilteredTable(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_filt[[2]], 
                                                                                   Project, ensEMBL2id,
                                                                                   ORIGIN1="Hu", ORIGIN2="Hu", GROUP1="CD39mKIRp", GROUP2="CD39mKIRm",
                                                                                   padj_cut, l2fc_cut, CV2_cut)

head(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_filt.tab)
head(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_filt.tab)
head(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_filt.tab)  


#
# FUNCTION: Add a column for CV2 value the group with the highest mean
#
functionMakeCV2MaxMean <- function(res.table){

  res.table$MaxMeanCV2 <- ifelse(res.table[,c(8)] > res.table[,c(12)], res.table[,c(7)],
                                 ifelse(res.table[,c(8)] < res.table[,c(12)], res.table[,c(11)],
                                        ifelse(res.table[,c(8)] == res.table[,c(12)], 0, NA)))
  
return(res.table)
}

res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.tab  <- functionMakeCV2MaxMean(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_filt.tab)
res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_maxmean.tab  <- functionMakeCV2MaxMean(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_filt.tab)
res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.tab  <- functionMakeCV2MaxMean(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_filt.tab)


write.csv(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.tab[order(abs(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.tab$log2FoldChange),decreasing=TRUE),], 
          file=paste0("DifferentialGeneExpression/", Project, "_DESeq2_DEGs_MaxMeanCV2Filtered", ".", "Hu_CD39pKIRp_vs_Hu_CD39mKIRm", '.ann.csv'))
write.csv(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_maxmean.tab[order(abs(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_maxmean.tab$log2FoldChange),decreasing=TRUE),], 
          file=paste0("DifferentialGeneExpression/", Project, "_DESeq2_DEGs_MaxMeanCV2Filtered", ".", "Hu_CD39pKIRp_vs_Hu_CD39mKIRp", '.ann.csv'))
write.csv(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.tab[order(abs(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.tab$log2FoldChange),decreasing=TRUE),], 
          file=paste0("DifferentialGeneExpression/", Project, "_DESeq2_DEGs_MaxMeanCV2Filtered", ".", "Hu_CD39mKIRp_vs_Hu_CD39mKIRm", '.ann.csv'))


#
#
#
test.data.ux                 <- res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_filt[[2]]
test.data.ux.ann             <- merge(test.data.ux,ensEMBL2id, by="ensembl_gene_id")
test.data.ux.ann$description <- gsub("..Source.*", "", test.data.ux.ann$description)
head(test.data.ux.ann)

test.data.ux.ann.sig.filt <- subset(test.data.ux.ann, CV2**2 <= CV2_cut & abs(log2FoldChange) >= l2fc_cut & padj <= padj_cut)
head(test.data.ux.ann.sig.filt)

#
# Plot log2FC vs CV2
#
pdf(paste0("ExploreDEG/", Project, "_DESeq2_DEGs_l2fc_", l2fc_cut, '_CV2_cut_', CV2_cut, "_l2fc_vs_cv2.pdf"),width=5,height=5, onefile=T)
par(bg=NA)
ggplot(test.data.ux.ann, aes(x=log2FoldChange, y=CV2**2) ) +
  geom_hline(yintercept = CV2_cut,     colour="firebrick", linetype = "dashed", alpha=1) +
  geom_vline(xintercept = l2fc_cut,    colour="firebrick", linetype = "dashed", alpha=1) +
  geom_vline(xintercept = -(l2fc_cut), colour="firebrick", linetype = "dashed", alpha=1) +
  geom_point(data=subset(test.data.ux.ann, CV2**2 <= CV2_cut & abs(log2FoldChange) >= l2fc_cut & padj > padj_cut), 
             colour="grey", alpha=0.5, size=0.5) +
  geom_point(data=subset(test.data.ux.ann, CV2**2 <= CV2_cut & abs(log2FoldChange) >= l2fc_cut & padj <= padj_cut), 
             colour="firebrick", alpha=0.5, size=0.5) +
  geom_point(data=subset(test.data.ux.ann, CV2**2 <= CV2_cut & abs(log2FoldChange) < l2fc_cut), colour="grey", alpha=0.5, size=0.5) +
  geom_point(data=subset(test.data.ux.ann, CV2**2 > CV2_cut & padj >  padj_cut & abs(log2FoldChange) <  l2fc_cut), 
             colour="grey", alpha=0.5, size=0.5) +
  geom_point(data=subset(test.data.ux.ann, CV2**2 > CV2_cut & padj <= padj_cut & abs(log2FoldChange) >= l2fc_cut), 
             colour="black", alpha=0.5, size=0.5) +
  geom_label_repel(data=subset(test.data.ux.ann, ensembl_gene_id=="ENSG00000173482" | ensembl_gene_id=="ENSG00000076706" | 
                                            ensembl_gene_id=="ENSG00000166105" | ensembl_gene_id=="ENSG00000196208" | 
                                            ensembl_gene_id=="ENSG00000197702"), 
                   aes(x=(log2FoldChange), y=CV2**2, label=external_gene_name),
                   fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                   size=3,  segment.color = 'firebrick', nudge_x = 0, nudge_y=0.5) +
  ylim(0,4.25) +
  xlab(bquote("log"[2]~"(Fold Change)")) +
  ylab(bquote("CV"^2)) 
dev.off()


#
# Plot log2FC vs padj
#
pdf(paste0("ExploreDEG/", Project, "_DESeq2_DEGs_l2fc_", l2fc_cut, '_CV2_cut_', CV2_cut, "_volcano_cv2.pdf"),width=7.5,height=7.5, onefile=T)
par(bg=NA)
ggplot(test.data.ux, aes(x=log2FoldChange, y=-log(padj)) ) +
  geom_vline(xintercept = l2fc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
  geom_vline(xintercept = -(l2fc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
  geom_hline(yintercept = -log(padj_cut), colour="black", linetype = "dashed", alpha=0.5) +
  geom_point(data=subset(test.data.ux, ( abs(log2FoldChange) < l2fc_cut | padj > 0.05 )), colour="grey", alpha=0.5, size=0.5) +
  geom_point(data=subset(test.data.ux, CV2**2 <= CV2_cut & abs(log2FoldChange) >= l2fc_cut & padj <= padj_cut), 
             colour="firebrick", alpha=0.5, size=0.5) +
  geom_point(data=subset(test.data.ux, CV2**2 > CV2_cut & abs(log2FoldChange) >= l2fc_cut & padj <= padj_cut), 
             colour="black", alpha=0.5, size=0.5) +
  
  geom_label_repel(#data=test.data.ux.lab.u[c(1:topN),], 
                   data=subset(test.data.ux.ann, ensembl_gene_id=="ENSG00000160219"),
                   aes(label=external_gene_name), fill='purple', colour='white', point.padding = unit(0.01, "lines"),  
                   size=2,  segment.color = 'black', nudge_x = 0, nudge_y=50) +
  
  geom_label_repel(data=test.data.ux.ann[c(1:topN),], 
                   aes(label=external_gene_name), fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                   size=2,  segment.color = 'firebrick', nudge_x = 0, nudge_y=5) +
  geom_label_repel(data=test.data.ux.ann[c(1:topN),], 
                   aes(label=external_gene_name), fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                   size=2,  segment.color = 'firebrick', nudge_x = 0, nudge_y=5) +
  ggtitle("Hu_CD39pKIRp_vs_Hu_CD39mKIRm") +
  xlab(bquote("log"[2]~"(Fold Change)")) + 
  ylab(bquote("-log"[10]~"(adj.p.value)")) 
dev.off()  


#
# Plot Rank vs CV2
#
test.data.u.ranked      <- test.data.ux.ann[order(test.data.ux.ann$CV2, decreasing=T),]

#subset(test.data.u.ranked, ensembl_gene_id=="ENSG00000173482" | ensembl_gene_id=="ENSG00000076706" | 
#         ensembl_gene_id=="ENSG00000166105" | ensembl_gene_id=="ENSG00000196208" | 
#         ensembl_gene_id=="ENSG00000197702")

test.data.u.ranked      <- subset(test.data.u.ranked, padj <= padj_cut & abs(log2FoldChange) > 2 )
test.data.u.ranked$Rank <- 1:nrow(test.data.u.ranked)
head(test.data.u.ranked)
tail(test.data.u.ranked)

CV2_rank_cut <- head(subset(test.data.u.ranked, CV2 < 2), 1)$Rank
pdf(paste0("ExploreDEG/", Project, "_DESeq2_DEGs_l2fc_", 2, '_CV2_cut_', CV2_cut, "_Rank_vs_cv2.pdf"),width=5,height=5, onefile=T)
par(bg=NA)
ggplot(test.data.u.ranked, aes(x=Rank, y=CV2**2) ) +
  geom_vline(xintercept = CV2_rank_cut, colour="firebrick", linetype = "dashed", alpha=0.5) +
  annotate("text", x = CV2_rank_cut, y = 0.1, label = paste0("rank=",CV2_rank_cut), color="firebrick") +
  geom_line(alpha=0.5, size=0.5) +
  geom_label_repel(data=subset(test.data.u.ranked, ensembl_gene_id=="ENSG00000173482" | ensembl_gene_id=="ENSG00000076706" | 
                                 ensembl_gene_id=="ENSG00000166105" | ensembl_gene_id=="ENSG00000196208" | 
                                 ensembl_gene_id=="ENSG00000197702"), 
                   aes(x=Rank, y=CV2**2, label=paste0(Group, "_", external_gene_name)),
                   fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                   size=3,  segment.color = 'firebrick', nudge_x = 500, nudge_y=-.5) +
  ylab(bquote("CV"^2)) +
  ggtitle(paste0("l2fc > ", 2, '; padj < ', padj_cut, "; removed ", round(100*(CV2_rank_cut/nrow(test.data.u.ranked)),2), "%" ))
dev.off()

#
# Calculate the range of count values
#

head(test.data.u.ranked)
head(t.table.undup.filt.m.ann)

ttt                       <- merge(test.data.u.ranked, t.table.undup.filt.m.ann, by=c("ensembl_gene_id", "Origin", "Group"))
ttt                       <- ttt[order(ttt$Rank, decreasing=F), ]
ttt[is.nan(ttt$CV2),]$CV2 <- 0
ttt.agg                   <- do.call(data.frame, aggregate(ttt$value, 
                                     by = list(ttt$ensembl_gene_id, ttt$Origin, ttt$Group),
                                     FUN=function(x) c(min=min(x), max=max(x))) )
colnames(ttt.agg)         <- c("ensembl_gene_id", "Origin", "Group", "min", "max")
head(ttt.agg)


#
# Plot macCV2Rank vs count
#
ttt.cv.agg            <- aggregate(ttt$CV2, by = list(ttt$ensembl_gene_id), max)
colnames(ttt.cv.agg)  <- c("ensembl_gene_id", "maxCV2")
ttt.cv.agg            <- ttt.cv.agg[order(ttt.cv.agg$maxCV2, decreasing=T), ]
ttt.cv.agg$maxCV2Rank <- 1:nrow(ttt.cv.agg)

head(ttt)
head(ttt.cv.agg)

ttt.cv                <- merge(ttt, ttt.cv.agg, by=c("ensembl_gene_id"))
ttt.cv                <- ttt.cv[order(ttt.cv$maxCV2, decreasing=T), ]

head(ttt.agg)
head(ttt.cv)

ttt.agg2              <- merge(ttt.agg, ttt.cv, by=c("ensembl_gene_id", "Origin", "Group"))
ttt.agg2              <- ttt.agg2[,c(1:5,12, 24)]
ttt.agg2              <- ttt.agg2[!duplicated(ttt.agg2),]
head(ttt.agg2)

pdf(paste0("ExploreDEG/", Project, "_DESeq2_DEGs_l2fc_", l2fc_cut, '_CV2_cut_', CV2_cut, "_maxCV2Rank_vs_counts.pdf"),width=25,height=5, onefile=T)
par(bg=NA)
ggplot(data=ttt.cv,
       aes(x=maxCV2Rank, y=log2(value), color=Group) ) +
  geom_vline(xintercept = CV2_rank_cut, colour="black", linetype = "dashed", size=1, alpha=0.75 ) +
  geom_point(size=0.5) +
  geom_segment( data=ttt.agg2, 
                aes(x=maxCV2Rank, xend=maxCV2Rank, y=log2(min), yend=log2(max), color=Group ) ) +
  geom_label_repel(data=subset(ttt.agg2, ensembl_gene_id=="ENSG00000173482" | ensembl_gene_id=="ENSG00000076706" | 
                                 ensembl_gene_id=="ENSG00000166105" | ensembl_gene_id=="ENSG00000196208" | 
                                 ensembl_gene_id=="ENSG00000197702"), 
                   aes(x=maxCV2Rank, y=log2(max), label=paste0(Group, "_", ensembl_gene_id)),
                   fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                   size=3,  segment.color = 'firebrick', nudge_x = 0, nudge_y=0) +
 # ylim(0,12.5) + 
  #xlim(0,500) +
  xlab(bquote("Rank (max CV"^2*")")) +
  ylab(bquote("log"[2]~"(normalised counts)")) 
dev.off()
  



message("+-------------------------------------------------------------------------------")
message("+ Gene Count plots per individual Gene specified                                ")
message("+-------------------------------------------------------------------------------")

makeGeneCountPlot <- function(DDS, ensEMBL2id, CONDITION, TITLE, gene2plot,outdir) {
  #
  # Plot the normalised read counts for a specified gene
  #
  if(missing(outdir)){ outdir = "" }
  else( outdir <- paste(outdir, "/", sep=""))
  
  genename2plot <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == gene2plot, ]$external_gene_name
  t2            <- plotCounts(DDS, gene=gene2plot, intgroup=c(CONDITION), normalized=TRUE, returnData=TRUE)
  colnames(t2) <- c("count", "condition")
  
  plt.col <- ggplot(t2, aes(x=condition, y=count, colour=condition, fill=condition)) + 
             geom_boxplot(width = 0.1, colour='black', fill='white', outlier.alpha = 0) + 
             geom_point(position=position_jitter(w=0.1,h=0), alpha=0.5, size=1) +
             ggtitle(paste0(Project, " ::: ", gene2plot, " ::: ", genename2plot)) + 
             xlab("") + ylab("Normalised count") +
             scale_fill_manual(name="Genotype", values = c("purple","purple4", "lightgreen", "green", "pink", "red", "cyan", "blue")) +
             theme(text = element_text(size=elementTextSize), axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") 
  
  pdf(paste0(outdir, Project, "-DGE_", gene2plot, "_", TITLE, "_collated.pdf"),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ plt.col })
  dev.off()
  
  png(paste0(outdir, Project, "-DGE_", gene2plot, "_", TITLE, "_collated.png"),units="cm", width=8, height=6, res=300)
  par(bg=NA)
  print({ 
     plt.col +
      theme(text=element_text(size=4), plot.title=element_text(size=6),
            axis.text.x = element_text(size=4), axis.text.y = element_text(size=4),
            legend.key.height=unit(0.5, "lines"), legend.key.width=unit(0.2, "lines") )
       })
  dev.off()
  
  
  t2$samples   <- rownames(t2)
  colnames(t2) <- c("count", "condition", "samples")
  t2           <- t2[order(t2$condition),]
  t2$samples2 <- factor(t2$samples, as.character(t2$samples))
  
  plt.ind <- ggplot(t2, aes(x=samples2, y=count, fill=condition, group=condition)) + 
             geom_bar(stat="identity", alpha=.5) + 
             scale_fill_manual(name="Comparison", values = c("purple","purple4", "lightgreen", "green", "pink", "red", "cyan", "blue")) +
             theme(text=element_text(size=elementTextSize), axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") +
             ggtitle(paste0(Project, " ::: ", gene2plot, " ::: ", genename2plot))
  
  pdf(paste0(outdir, Project, "-DGE_", gene2plot, "_", TITLE, "_individual.pdf"),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ plt.ind })
  dev.off()
  
  png(paste0(outdir, Project, "-DGE_", gene2plot, "_", TITLE, "_individual.png"),units="cm", width=8, height=6, res=300)
  par(bg=NA)
  print({ 
      plt.ind +
      theme(text=element_text(size=4), plot.title=element_text(size=6),
            axis.text.x = element_text(size=4), axis.text.y = element_text(size=4),
            legend.key.height=unit(0.5, "lines"), legend.key.width=unit(0.2, "lines") ) })
  dev.off()
  print(paste("Created plot for", gene2plot), sep=" ")
}

makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000173482', paste0(baseDir,"/ExploreDEG")) # PTPRM
makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000076706', paste0(baseDir,"/ExploreDEG")) # AP002956.1
makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000166105', paste0(baseDir,"/ExploreDEG")) # GLB1L3
makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000196208', paste0(baseDir,"/ExploreDEG")) # GREB1
makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000197702', paste0(baseDir,"/ExploreDEG")) # PARVA

makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000115523', paste0(baseDir,"/ExploreDEG")) 

makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000000938', paste0(baseDir,"/ExploreDEG"))
makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000243772', paste0(baseDir,"/ExploreDEG"))

makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000176083', paste0(baseDir,"/ExploreDEG"))
makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000161249', paste0(baseDir,"/ExploreDEG"))


makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000176083', paste0(baseDir,"/ExploreDEG")) # ZNF683
makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000242473', paste0(baseDir,"/ExploreDEG")) # 


makeGeneCountPlot(dds.Samples, ensEMBL2id, "Group", "Group", 'ENSG00000138185', paste0(baseDir,"/ExploreDEG")) # ENTPD1


message("+-------------------------------------------------------------------------------")
message("+ Make Custom MA and Volcano Plots")
message("+-------------------------------------------------------------------------------")

#head(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_filt.tab)

make_MA_Vcn_Plots <- function(res.table.ann, Project, l2fc_cut, CV2_cut, padj_cut, topN, TITLE, outDir,xrange,yrange) {
  
   res.table.ann$GRP1_CV2 <- res.table.ann[,c(7)]
   res.table.ann$GRP2_CV2 <- res.table.ann[,c(11)]
   
   res.table.ann.ord.u    <- res.table.ann[order(res.table.ann$log2FoldChange,decreasing=T),] 
   res.table.ann.ord.u    <- subset(res.table.ann.ord.u, log2FoldChange >= l2fc_cut & GRP1_CV2 < CV2_cut &  GRP2_CV2 < CV2_cut & padj <= padj_cut)
   res.table.ann.ord.d    <- res.table.ann[order(res.table.ann$log2FoldChange,decreasing=F),] 
   res.table.ann.ord.d    <- subset(res.table.ann.ord.d, log2FoldChange <= l2fc_cut & GRP1_CV2 < CV2_cut &  GRP2_CV2 < CV2_cut & padj <= padj_cut)

   plt.vlc <- ggplot(data=res.table.ann, aes(x=log2FoldChange, y=-log(padj), label=external_gene_name)) +
                     geom_vline(xintercept = l2fc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
                     geom_vline(xintercept = -(l2fc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
                     geom_hline(yintercept = -log(padj_cut), colour="black", linetype = "dashed", alpha=0.5) +
              geom_point(data=subset(res.table.ann, abs(log2FoldChange) < l2fc_cut | padj > padj_cut), 
                         alpha=0.75, size=1.0, colour="grey85") +
              geom_point(data=subset(res.table.ann, padj<=padj_cut & abs(log2FoldChange) >= l2fc_cut & (GRP1_CV2 >= CV2_cut |  GRP2_CV2 >= CV2_cut)),      
                         alpha=0.75, size=1.0, colour="grey50") +
              geom_point(data=subset(res.table.ann, padj<=padj_cut & abs(log2FoldChange) >= l2fc_cut & GRP1_CV2 < CV2_cut &  GRP2_CV2 < CV2_cut),      
                         alpha=0.75, size=1.0, colour="firebrick") +
              geom_label_repel(data=res.table.ann.ord.u[c(1:topN),], fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                               size=3,  segment.color = 'firebrick', nudge_x = 0, nudge_y=5) +
              geom_label_repel(data=res.table.ann.ord.d[c(1:topN),], fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                               size=3,  segment.color = 'firebrick', nudge_x = 0, nudge_y=5) +
              xlab(bquote("log"[2]~"(Fold Change)")) + 
              ylab(bquote("-log"[10]~"(adj.p.value)")) + 
           #   scale_x_continuous(limits=c(xrange[1],xrange[2]), breaks=seq(xrange[1],xrange[2],xrange[3])) + 
           #   scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
              ggtitle(TITLE) +
              theme(aspect.ratio=1,
                    text=element_text(size=elementTextSize), plot.title=element_text(size=elementTextSize),
                    axis.text.x = element_text(size=elementTextSize), axis.text.y = element_text(size=elementTextSize))

   pdf(paste0(outDir, Project, "_DESeq2_DEGs_padj_", padj_cut, '_l2fc_', l2fc_cut, '_CV2_', CV2_cut, "_", TITLE, "_Volcanoplot.pdf"),width=5,height=5, onefile=T)
   par(bg=NA)
   print({ plt.vlc })
   dev.off()

return(plt.vlc)
}


outDir  <- paste0(baseDir, "/DifferentialGeneExpression/")
topN    <- 5
CV2_cut <- 1.95
res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.MA_Vcn  <- make_MA_Vcn_Plots(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_filt.tab,  Project, l2fc_cut, CV2_cut, padj_cut, topN, TITLE="Hu_CD39pKIRp_vs_Hu_CD39mKIRm",  outDir, xrange=c(-10,10,5), yrange=c(0,186,50))
res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.MA_Vcn  <- make_MA_Vcn_Plots(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_filt.tab,  Project, l2fc_cut, CV2_cut, padj_cut, topN, TITLE="Hu_CD39pKIRp_vs_Hu_CD39mKIRp",  outDir, xrange=c(-10,10,5), yrange=c(0,186,50))
res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.MA_Vcn  <- make_MA_Vcn_Plots(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_filt.tab,  Project, l2fc_cut, CV2_cut, padj_cut, topN, TITLE="Hu_CD39mKIRp_vs_Hu_CD39mKIRm",  outDir, xrange=c(-10,10,5), yrange=c(0,186,50))


make_MA_Vcn_Plots_V2 <- function(res.table.ann, Project, l2fc_cut, CV2_cut, padj_cut, reads_cut, topN, TITLE, outDir, figVer, xrange, yrange, interestingGenes) {
  
  res.table.ann.ord.u    <- res.table.ann[order(res.table.ann$log2FoldChange,decreasing=T),] 
  res.table.ann.ord.u    <- subset(res.table.ann.ord.u, log2FoldChange >= l2fc_cut & MaxMeanCV2 < CV2_cut & padj <= padj_cut)
  res.table.ann.ord.u    <- res.table.ann.ord.u[c(1:topN),]
  res.table.ann.ord.d    <- res.table.ann[order(res.table.ann$log2FoldChange,decreasing=F),] 
  res.table.ann.ord.d    <- subset(res.table.ann.ord.d, log2FoldChange <= l2fc_cut & MaxMeanCV2 < CV2_cut & padj <= padj_cut)
  res.table.ann.ord.d    <- res.table.ann.ord.d[c(1:topN),]
  res.table.ann.intr     <- res.table.ann[res.table.ann$external_gene_name %in% interestingGenes,] 
  
  colnames(res.table.ann)[8]  <- "Mean.x"
  colnames(res.table.ann)[12] <- "Mean.y"
  
  plt.vlc <- ggplot(data=res.table.ann, aes(x=log2FoldChange, y=-log(padj), label=external_gene_name)) +
                    geom_vline(xintercept = l2fc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
                    geom_vline(xintercept = -(l2fc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
                    geom_hline(yintercept = -log(padj_cut), colour="black", linetype = "dashed", alpha=0.5) +
                    geom_point(data=subset(res.table.ann, abs(log2FoldChange) < l2fc_cut | padj > padj_cut), 
                               alpha=0.75, size=1.0, colour="grey85") +
                    geom_point(data=subset(res.table.ann, padj<=padj_cut & abs(log2FoldChange) >= l2fc_cut & (MaxMeanCV2 >= CV2_cut | (Mean.x < reads_cut & Mean.y < reads_cut))),      
                               alpha=0.75, size=1.0, colour="grey40") +
                    geom_point(data=subset(res.table.ann, padj<=padj_cut & abs(log2FoldChange) >= l2fc_cut & MaxMeanCV2 < CV2_cut & (Mean.x >= reads_cut | Mean.y >= reads_cut) ),      
                               alpha=0.75, size=1.0, colour="firebrick") 
  
  if(topN < 1){ plt.vlc  +
                geom_label_repel(data=res.table.ann.ord.u, fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                                 size=3,  segment.color = 'firebrick', nudge_x = 0, nudge_y=5) +
                geom_label_repel(data=res.table.ann.ord.d, fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                                 size=3,  segment.color = 'firebrick', nudge_x = 0, nudge_y=5)  }
    
  plt.vlc <- plt.vlc  + 
             geom_label_repel(data=res.table.ann.intr, fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                              size=2.5, min.segment.length = unit(0.25, 'lines'), segment.color = 'firebrick', nudge_x = 0, nudge_y=5, force=10) +
    
             xlab(bquote("log"[2]~"(Fold Change)")) + 
             ylab(bquote("-log"[10]~"(adj.p.value)")) + 
             scale_x_continuous(limits=c(xrange[1]+3,xrange[2]-3), breaks=seq(xrange[1],xrange[2],xrange[3])) + 
             #scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
             ggtitle(TITLE) +
             theme(aspect.ratio=1,
                   text=element_text(size=elementTextSize), plot.title=element_text(size=elementTextSize),
                   axis.text.x = element_text(size=elementTextSize), axis.text.y = element_text(size=elementTextSize))
  
  pdf(paste0(outDir, Project, "_DESeq2_DEGs_padj_", padj_cut, '_l2fc_', l2fc_cut, '_CV2_', CV2_cut, "_read_", reads_cut, "_", TITLE, "_Volcanoplot_V", figVer, ".pdf"),width=5,height=5, onefile=T)
  par(bg=NA)
  print({ plt.vlc })
  dev.off()
  
  return(plt.vlc)
}

interesting_pp_mm <- c('PDGFD','LGALS9','ENTPD1', 'CXCR4','CD44','KIR2DL3','KIR2DL1','EGLN3','GPR15','IL18R1','ICOS','CDKN2A', 'GNLY', 'GZMA', 'CCL5')
interesting_pp_mp <- c('TEF', 'GFPT2', 'TMRT12', 'RAB38','NKIRAS1', 'PDGFD', 'PRKCQ', 'IL2RA', 'ZNF683', 'ACVR1B', 'INPP4B', 'DPP4', 'TNFSF14')
interesting_mp_mm <- c('CXCR6', 'KLRG1', 'ICOS', 'IL7R', 'IKZF2', 'KIR2DL3', 'KIR2DL1', 'LILRB1', 'KIR2DS4', 'CDKN2A', 'GNLY')
topN              <- 0
figVer            <- 3
l2fc_cut          <- 1.5 #3

res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.MA_Vcn  <- make_MA_Vcn_Plots_V2(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.tab,  
                                                                                 Project, l2fc_cut, CV2_cut, padj_cut, reads_cut, topN, 
                                                                                 TITLE="Hu_CD39pKIRp_vs_Hu_CD39mKIRm",  outDir, figVer,
                                                                                 xrange=c(-15,15,5), yrange=c(0,186,50), interesting_pp_mm)
res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_maxmean.MA_Vcn  <- make_MA_Vcn_Plots_V2(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_maxmean.tab,  
                                                                                 Project, l2fc_cut, CV2_cut, padj_cut, reads_cut, topN, 
                                                                                 TITLE="Hu_CD39pKIRp_vs_Hu_CD39mKIRp",  outDir, figVer,
                                                                                 xrange=c(-15,15,5), yrange=c(0,186,50), interesting_pp_mp)
res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.MA_Vcn  <- make_MA_Vcn_Plots_V2(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.tab,  
                                                                                 Project, l2fc_cut, CV2_cut, padj_cut, reads_cut, topN, 
                                                                                 TITLE="Hu_CD39mKIRp_vs_Hu_CD39mKIRm",  outDir, figVer,
                                                                                 xrange=c(-15,15,5), yrange=c(0,186,50), interesting_mp_mm)







message("+-------------------------------------------------------------------------------")
message("+ Heatmaps")
message("+-------------------------------------------------------------------------------")


l2fc_cut  <- 3     #2
padj_cut  <- 0.05  #0.01
reads_cut <- 10

annotation.df              <- as.data.frame(colData(dds.Samples)[,c("Group","Origin")])
annotation.df$Group        <- gsub("CD39pKIRp", "CD39+KIR+", annotation.df$Group)
annotation.df$Group        <- gsub("CD39mKIRm", "CD39-KIR-", annotation.df$Group)
annotation.df$Group        <- gsub("CD39mKIRp", "CD39-KIR+", annotation.df$Group)

AnnotCols                  <- list( Origin=c("Hu"="lightgrey"), Group=c("CD39+KIR+"="pink", "CD39-KIR-"="darkgreen", "CD39-KIR+"="purple"))

count.data                 <- as.data.frame(counts(dds.Samples,normalized=TRUE))
count.data$ensembl_gene_id <- rownames(count.data)

res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.ann.sub <- subset(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.tab, abs(log2FoldChange)>=l2fc_cut & MaxMeanCV2<=CV2_cut & padj<=padj_cut)
colnames(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.ann.sub) <- c("ensembl_gene_id","log2FoldChange","lfcSE","padj","Origin.G1","Group.G1","G1_CV2","G1_Mean","Origin.G2","Group.G2","G2_CV2","G2_Mean","external_gene_name", "description", "chromosome_name", "MaxMeanCV2")

res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.ann.sub <- subset(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.cv2_maxmean.tab, abs(log2FoldChange)>=l2fc_cut & MaxMeanCV2<=CV2_cut & padj<=padj_cut)
colnames(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.ann.sub) <- c("ensembl_gene_id","log2FoldChange","lfcSE","padj","Origin.G1","Group.G1","G1_CV2","G1_Mean","Origin.G2","Group.G2","G2_CV2","G2_Mean","external_gene_name", "description", "chromosome_name", "MaxMeanCV2")

res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.ann.sub <- subset(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.cv2_maxmean.tab, abs(log2FoldChange)>=l2fc_cut & MaxMeanCV2<=CV2_cut & padj<=padj_cut)
colnames(res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.ann.sub) <- c("ensembl_gene_id","log2FoldChange","lfcSE","padj","Origin.G1","Group.G1","G1_CV2","G1_Mean","Origin.G2","Group.G2","G2_CV2","G2_Mean","external_gene_name", "description", "chromosome_name", "MaxMeanCV2")


hm.hu.only.DEGs  <- rbind(res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRm.ann.sub,  
                          res.o_g.Hu_CD39pKIRp_vs_Hu_CD39mKIRp.ann.sub,  
                          res.o_g.Hu_CD39mKIRp_vs_Hu_CD39mKIRm.ann.sub)
#subset(hm.hu.only.DEGs, ensembl_gene_id=="ENSG00000138185")

res.data                   <- count.data[count.data$ensembl_gene_id %in% unique(hm.hu.only.DEGs$ensembl_gene_id),]
hm.data.filt               <- merge(res.data,ensEMBL2id , by="ensembl_gene_id" )
rownames(hm.data.filt)     <- hm.data.filt$external_gene_name
samp.list                  <- (rownames(subset(annotation.df, Origin == "Hu" & (Group == "CD39+KIR+" | Group == "CD39-KIR+" | Group == "CD39-KIR-") )))
hm.data.filt.sel           <- hm.data.filt[, samp.list]
hm.data.filt.sel           <- as.matrix( log2( hm.data.filt.sel + 1) )
numDEG                     <- nrow(hm.data.filt.sel)

annotation.df <- annotation.df[,c("Group"),drop=F]
print(annotation.df)

means.tab.1        <- hm.data.filt[, samp.list]
means.tab.1        <- means.tab.1[,c("HU12_P11_01", "HU13_P11_05", "Hu05_P11_09", "Hu7_P12_14", 
                                     "HU12_P9_02",  "HU13_P9_06",  "Hu05_P9_10",  "Hu7_P9_16",
                                     "HU12_P10_03", "HU13_P10_07", "Hu05_P10_11", "Hu7_P13_15")]
means.tab.1$MeanPP <- rowMeans( means.tab.1[,c("HU12_P11_01", "HU13_P11_05", "Hu05_P11_09", "Hu7_P12_14")])
means.tab.1$MeanMM <- rowMeans( means.tab.1[,c("HU12_P9_02",  "HU13_P9_06",  "Hu05_P9_10",  "Hu7_P9_16")])
means.tab.1$MeanMP <- rowMeans( means.tab.1[,c("HU12_P10_03", "HU13_P10_07", "Hu05_P10_11", "Hu7_P13_15")])
means.tab.1        <- means.tab.1[,c("MeanPP","MeanMM","MeanMP")]

write.csv(means.tab.1, 
          file=paste0(outDir, Project, "_DEGs_padj_", padj_cut, '_l2fc_', l2fc_cut, '_CV2_', CV2_cut, 
                      "_HeatMapData_Hu_CD39pKIRp_vs_CD39mKIRm_vs_CD39mKIRp_V1.csv"))

head(means.tab.1, 10)


CairoPDF(paste0(outDir, Project, "_DEGs_padj_", padj_cut, '_l2fc_', l2fc_cut, '_CV2_', CV2_cut, 
                "_HeatMap_Hu_CD39pKIRp_vs_CD39mKIRm_vs_CD39mKIRp_V3.pdf"),width=7.5,height=20, onefile=T)
par(bg=NA)
pheatmap(hm.data.filt.sel, cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col=annotation.df, annotation_colors=AnnotCols, 
         fontsize_col=5, fontsize_row=4, color = colorRampPalette(c("dodgerblue4", "white", "firebrick"))(25),
         cutree_cols=3, cutree_rows = 8,
         main=paste0("Hu (CD39+KIR+ Vs CD39-KIR- Vs CD39-KIR+)", "\n#Genes=", numDEG, " [l2FC>=", l2fc_cut, " MaxMeanCV2<=", CV2_cut, " padj<=", padj_cut, "]"))
dev.off()


message("+-------------------------------------------------------------------------------")
message("+ GO Term Enrichment Vs Expression")
message("+-------------------------------------------------------------------------------")

#
# To Be Aded
#
message("+-------------------------------------------------------------------------------")
message("+ END OF SCRIPT ")
message("+-------------------------------------------------------------------------------")