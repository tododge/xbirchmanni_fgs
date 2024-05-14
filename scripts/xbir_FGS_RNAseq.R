##FGS RNAseq

#install.packages(c("readr","ashr"))
#BiocManager::install(c("tximportData","tximport","GenomicFeatures","DESeq2"))
BiocManager::install(c("EnhancedVolcano"))
BiocManager::install("PCAtools")

library(tidyverse)
library(tximportData)
library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)
library(ggrepel)
library(PCAtools)
library(ggplot2)
library(ashr)

species<- "birchmanni"

## Read in sample files
dir <- "/Users/trisdodge/Desktop/Stanford/falsegravid/RNAseq"

samples <- read.table(file.path(dir, "xbir_juv_samples.txt"), header = TRUE)

# make sure that non-continuous variables are cast as factors
samples$phenotype <- factor(samples$phenotype)
samples$extraction_batch <- factor(samples$extraction_batch)

# subset data by species

focal_tissue <- "PM_EAM"
#focal_tissue <- "Flank"
focal_tissue <- "Brain"

samples <- samples[samples$tissue == focal_tissue, ]


samples

files <- file.path(dir, "birchmanni/input_files", paste(samples$prefix,"kallisto",sep="_"), "abundance.h5")
names(files) <- paste0(samples$prefix)
files

## load transcript annotations
txdb <- makeTxDbFromGFF(file=file.path(dir, species, "xbir-COAC-16-VIII-22-M_v2023.1.tx.transcriptDESeq2.gff"), format="gff3")

# create tx2gene table, match transcript name to geneid, doing this 
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# get gene-level count information from kallisto counts

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)
head(txi$counts)

#################### DGE analysis ####################

# DESeq with simple model including FGS phenotype and extraction batch to control for batch effects

dds     <- DESeqDataSetFromTximport(txi,
                                    colData = samples,
                                    design = ~ phenotype + extraction_batch)

# Set WT (non-FGS) as reference
dds$phenotype <- relevel(dds$phenotype, ref = "WT")

# Determine best fit
# local: 0.0716
dds <- DESeq(dds, fitType="local")
resultsNames(dds) #IMPORTANT! Note the results names printed here, you will use it in the step "compare expression between treatments" below

#plotDispEsts(dds, main="Dispersion plot with local fit")
#residual <- mcols(dds)$dispGeneEst - mcols(dds)$dispFit
#absres<-abs(residual)
#summary(absres) 
# local: 0.0716 better fit
#dds <- DESeq(dds, fitType="parametric")
#plotDispEsts(dds, main="Dispersion plot with parametric fit")
#residual <- mcols(dds)$dispGeneEst - mcols(dds)$dispFit
#absres<-abs(residual)
#summary(absres)  
# parametric: 0.1229 worse fit

# save dds object
saveRDS(dds, paste0(dir,"/xbir_fgs_dge_",focal_tissue,"_dds.rds"))

# save vst object
vst <- vst(dds,blind=FALSE)
saveRDS(vst, paste0(dir,"/xbir_fgs_dge_",focal_tissue,"_vst.rds"))

#load them in
dds_PMEAM <- readRDS("/Users/trisdodge/Desktop/Stanford/falsegravid/RNAseq/xbir_fgs_dge_PM_EAM_dds.rds")
vst_PMEAM <- readRDS("/Users/trisdodge/Desktop/Stanford/falsegravid/RNAseq/xbir_fgs_dge_PM_EAM_vst.rds")   
dds_Flank <- readRDS("/Users/trisdodge/Desktop/Stanford/falsegravid/RNAseq/xbir_fgs_dge_Flank_dds.rds")
vst_Flank <- readRDS("/Users/trisdodge/Desktop/Stanford/falsegravid/RNAseq/xbir_fgs_dge_Flank_vst.rds")   
dds_Brain <- readRDS("/Users/trisdodge/Desktop/Stanford/falsegravid/RNAseq/xbir_fgs_dge_Brain_dds.rds")
vst_Brain <- readRDS("/Users/trisdodge/Desktop/Stanford/falsegravid/RNAseq/xbir_fgs_dge_Brain_vst.rds")   


res.wtvsfgs_PMEAM <- lfcShrink(dds_PMEAM, coef="phenotype_FGS_vs_WT", type="ashr")
res.wtvsfgs_Flank <- lfcShrink(dds_Flank, coef="phenotype_FGS_vs_WT", type="ashr")
res.wtvsfgs_Brain <- lfcShrink(dds_Brain, coef="phenotype_FGS_vs_WT", type="ashr")

res.wtvsfgs_PMEAM <- merge(base::as.data.frame(res.wtvsfgs_PMEAM), base::as.data.frame(counts(dds_PMEAM, normalized=TRUE)), by="row.names", sort=FALSE)
names(res.wtvsfgs_PMEAM)[1] <- "Gene"
res.wtvsfgs_PMEAM$col <- ifelse(res.wtvsfgs_PMEAM$log2FoldChange>0, "#55185D", "darkgrey")

res.wtvsfgs_Flank <- merge(base::as.data.frame(res.wtvsfgs_Flank), base::as.data.frame(counts(dds_Flank, normalized=TRUE)), by="row.names", sort=FALSE)
names(res.wtvsfgs_Flank)[1] <- "Gene"
res.wtvsfgs_Flank$col <- ifelse(res.wtvsfgs_Flank$log2FoldChange>0, "#55185D", "darkgrey")

res.wtvsfgs_Brain <- merge(base::as.data.frame(res.wtvsfgs_Brain), base::as.data.frame(counts(dds_Brain, normalized=TRUE)), by="row.names", sort=FALSE)
names(res.wtvsfgs_Brain)[1] <- "Gene"
res.wtvsfgs_Brain$col <- ifelse(res.wtvsfgs_Brain$log2FoldChange>0, "#55185D", "darkgrey")

#res.wtvsfgs <- res.wtvsfgs_PMEAM
#dds <- dds_PMEAM

#res.wtvsfgs <- res.wtvsfgs_Flank
#dds <- dds_Flank

#res.wtvsfgs <- res.wtvsfgs_Brain
#dds <- dds_Brain

####Figure 3####
fig3a <- 
  ggplot() +
  geom_hline(yintercept= -log10(0.05), size = 0.5, lty = "dashed", alpha=0.4) +
  geom_point(data=subset(res.wtvsfgs_PMEAM,Gene != "g2346"), 
             aes(x=log2FoldChange, y=-log10(padj), 
                 color = col)) +
  geom_point(data=subset(res.wtvsfgs_PMEAM,Gene == "g2346"), 
             aes(x=log2FoldChange, y=-log10(padj), 
                 color = col), size=2, pch=8) +
  scale_x_continuous(limits = c(-7,7), breaks = seq(-10,10,2.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,45), breaks = seq(0,50,10)) +
  scale_color_manual(values = c("#55185D","darkgrey")) +
  xlab(expression(paste("log"[2],"(expression fold change)",sep="")))+
  ylab(expression(paste("-log"[10],"(",italic(""*P*""),"-value)",sep="")))+
  theme_bw(base_size = 7) + theme(legend.position = "none",
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())

fig3b <- 
  ggplot() +
  geom_hline(yintercept= -log10(0.05), size = 0.5, lty = "dashed", alpha=0.4) +
  geom_point(data=subset(res.wtvsfgs_Brain,Gene != "g2346"), 
             aes(x=log2FoldChange, y=-log10(padj), 
                 color = col)) +
  geom_point(data=subset(res.wtvsfgs_Brain,Gene == "g2346"), 
             aes(x=log2FoldChange, y=-log10(padj), 
                 color = col), size=2, pch=8) +
  scale_x_continuous(limits = c(-7,7), breaks = seq(-10,10,2.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,45), breaks = seq(0,50,10)) +
  scale_color_manual(values = c("#55185D","darkgrey")) +
  xlab(expression(paste("log"[2],"(expression fold change)",sep="")))+
  ylab(expression(paste("-log"[10],"(",italic(""*P*""),"-value)",sep="")))+
  theme_bw(base_size = 7) + theme(legend.position = "none",
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())


#Fig S12
  ggplot() +
  geom_hline(yintercept= -log10(0.05), size = 0.5, lty = "dashed", alpha=0.4) +
  geom_point(data=subset(res.wtvsfgs_Flank,Gene != "g2346"), 
             aes(x=log2FoldChange, y=-log10(padj), 
                 color = col)) +
  geom_point(data=subset(res.wtvsfgs_Flank,Gene == "g2346"), 
             aes(x=log2FoldChange, y=-log10(padj), 
                 color = col), size=2, pch=8) +
  #scale_x_continuous(limits = c(-7,7), breaks = seq(-10,10,2.5), expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,45), breaks = seq(0,50,10)) +
  scale_color_manual(values = c("#55185D","darkgrey")) +
  xlab(expression(paste("log"[2],"(expression fold change)",sep="")))+
  ylab(expression(paste("-log"[10],"(",italic(""*P*""),"-value)",sep="")))+
  theme_bw(base_size = 10) + theme(legend.position = "none",
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())



####PCA####

xPE <- assay(vst_PMEAM)
colnames(xPE)

rownames(samples) <- samples$prefix
rownames(samples)

all(colnames(xPE) == rownames(samples))
#remotes::install_version("matrixStats", version="1.1.0")
library(matrixStats)
p_PMEAM <- PCAtools::pca(xPE, metadata = samples, removeVar = 0.1)
#pairsplot(p_PMEAM)

biplot(p_PMEAM, x="PC1", y="PC2",
       lab = as.character(p$metadata$individual),
       colby = 'phenotype', colkey=c('FGS'="purple", 'WT'="darkgrey"),
       legendPosition = 'right', legendLabSize = 12, legendIconSize = 6.0,
       gridlines.major=FALSE, gridlines.minor=FALSE,
       drawConnectors = FALSE,
       title = 'EAM PCA',
       subtitle = 'PC1 versus PC2',
       caption = '')

xB <- assay(vst_Brain)
colnames(xB)

rownames(samples) <- samples$prefix
rownames(samples)

all(colnames(xB) == rownames(samples))

p_Brain <- PCAtools::pca(xB, metadata = samples, removeVar = 0.1)

biplot(p_Brain, x="PC1", y="PC2",
       lab = as.character(p$metadata$individual),
       colby = 'phenotype', colkey=c('FGS'="purple", 'WT'="darkgrey"),
       legendPosition = 'right', legendLabSize = 12, legendIconSize = 6.0,
       gridlines.major=FALSE, gridlines.minor=FALSE,
       drawConnectors = FALSE,
       title = 'Brain PCA',
       subtitle = 'PC1 versus PC2',
       caption = '')

#################### START PRODUCE OUT TABLES ####################

#g2346 kitlga

#res.wtvsfgs <- res.wtvsfgs[order(res.wtvsfgs$pvalue),]
#res.wtvsfgs["g2346",]

gff <- read.csv(file.path(dir, species, "xbir-COAC-16-VIII-22-M_v2023.1.tx.transcriptDESeq2.gff"), 
                header=FALSE, 
                sep="\t", quote = "", 
                row.names = NULL, 
                stringsAsFactors = FALSE)
# Extract gene_id and symbol using regex
gene_match <- gff[,c(1,9)]
gene_match$gene_id <- str_match(gene_match$V9, 'gene_id=\\s*(.*?)\\s*;')[,2]
gene_match$symbol <- str_match(gene_match$V9, 'symbol="\\s*(.*?)\\s*"')[,2]
gene_match$name <- str_match(gene_match$V9, 'name="\\s*(.*?)\\s*"')[,2]


sigdge_PMEAM <- subset(res.wtvsfgs_PMEAM, padj < 0.05) # 75 
sigdge_PMEAM[sigdge_PMEAM$Gene=="g2346",]

sigdge_Flank <- subset(res.wtvsfgs_Flank, padj < 0.05) # 85 
sigdge_Flank[sigdge_Flank$Gene=="g2346",]

sigdge_Brain <- subset(res.wtvsfgs_Brain, padj < 0.05) # 38 
sigdge_Brain[sigdge_Brain$Gene=="g2346",]

sigdge_PMEAM_names <- merge(sigdge_PMEAM, gene_match, by.x="Gene", by.y="gene_id")
sigdge_Flank_names <- merge(sigdge_Flank, gene_match, by.x="Gene", by.y="gene_id")
sigdge_Brain_names <- merge(sigdge_Brain, gene_match, by.x="Gene", by.y="gene_id")


sigdge_PMEAM_names <- sigdge_PMEAM_names[c("Gene","log2FoldChange","padj","symbol")][order(sigdge_PMEAM_names$padj), ]
sigdge_PMEAM_names$tissue <- "PM+EAM"

sigdge_Flank_names <- sigdge_Flank_names[c("Gene","log2FoldChange","padj","symbol")][order(sigdge_Flank_names$padj), ]
sigdge_Flank_names$tissue <- "Flank"

sigdge_Brain_names <- sigdge_Brain_names[c("Gene","log2FoldChange","padj","symbol")][order(sigdge_Brain_names$padj), ]
sigdge_Brain_names$tissue <- "Brain"

sigdge_all_tissues <- rbind(sigdge_PMEAM_names,sigdge_Flank_names,sigdge_Brain_names)


write.table(sigdge_all_tissues, paste0(dir,"/sigdge_alltissues_names.txt"), quote=FALSE, sep="\t", row.names=FALSE)


res.wtvsfgs_PMEAM

res.wtvsfgs_PMEAM$gene_num <- as.numeric(gsub("g","",res.wtvsfgs_PMEAM$Gene))

res.wtvsfgs_PMEAM_1mb <- subset(res.wtvsfgs_PMEAM, gene_num >= 2309 & gene_num <= 2366)
res.wtvsfgs_PMEAM_300kb <- subset(res.wtvsfgs_PMEAM, gene_num >= 2341 & gene_num <= 2352)

res.wtvsfgs_PMEAM_300kb <- merge(res.wtvsfgs_PMEAM_300kb, gene_match, by.x="Gene", by.y="gene_id")

res.wtvsfgs_PMEAM_300kb$mean_WT <- rowMeans(res.wtvsfgs_PMEAM_300kb[c("Xbir-COAC-16-VIII-23-S67-M03-WT-FGS","Xbir-COAC-3-X-23-M02-WT-FGS","Xbir-COAC-6-IX-23-S67-M01-WT-FGS","Xbir-COAC-6-IX-23-S67-M02-WT-FGS","Xbir-COAC-7-IX-23-M06-WT-FGS","Xbir-COAC-7-IX-23-M04-WT-FGS")])

res.wtvsfgs_PMEAM_300kb_long <- res.wtvsfgs_PMEAM_300kb %>%
  pivot_longer(
    cols = starts_with("Xbir"),
    names_to = "sample",
    values_to = "expression",
    values_drop_na = TRUE
  )
res.wtvsfgs_PMEAM_300kb_long$phenotype <- str_match(res.wtvsfgs_PMEAM_300kb_long$sample, ".*(.{2,3})-FGS$")[,2]
res.wtvsfgs_PMEAM_300kb_long$cor_expr <- res.wtvsfgs_PMEAM_300kb_long$expression/res.wtvsfgs_PMEAM_300kb_long$mean_WT

res.wtvsfgs_PMEAM_300kb_long$symbol <- gsub("KITLG", "kitlga",res.wtvsfgs_PMEAM_300kb_long$symbol)
res.wtvsfgs_PMEAM_300kb_long$symbol <- gsub("-", "rootletin-like",res.wtvsfgs_PMEAM_300kb_long$symbol)

res.wtvsfgs_PMEAM_300kb_long$symbol <- factor(res.wtvsfgs_PMEAM_300kb_long$symbol, levels = c("urah","cep41","poc1b","DUSP6","SLC23A2","kitlga","TMTC3","cep290","zgc:103499","Lyve1","rootletin-like","EIF4G2"))
res.wtvsfgs_PMEAM_300kb_long$phenotype <- factor(res.wtvsfgs_PMEAM_300kb_long$phenotype, levels = c("WT","GS"))

fig3c <- 
  ggplot() +
  geom_hline(yintercept=1, size = 0.5, lty = "dashed", alpha=0.4) +
  geom_boxplot(data=subset(res.wtvsfgs_PMEAM_300kb_long), 
               aes(x=symbol, y=cor_expr, color=phenotype), outlier.shape = NA) +
  geom_point(data=subset(res.wtvsfgs_PMEAM_300kb_long, Gene != "g2346"), 
             aes(x=symbol, y=cor_expr, color=phenotype), position = position_jitterdodge(), alpha=0.2) +
  geom_point(data=subset(res.wtvsfgs_PMEAM_300kb_long, Gene == "g2346"), 
             aes(x=symbol, y=cor_expr, color=phenotype),
             position = position_jitterdodge(),
             alpha=1, pch=8) +
  scale_color_manual(values=c("darkgrey", "#55185D")) +
  xlab("genes within 300kb of GWAS peak") +
  ylab("expression fold change") +
  theme_bw(base_size = 7) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 5))
