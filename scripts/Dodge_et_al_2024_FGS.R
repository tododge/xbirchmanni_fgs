##################
#####FIGURE 1#####
##################

####GWAS PLOTS AND STRUCTURAL REARANGEMENT####
library(tidyverse)
library(cowplot)

gwas_xbir_WT_FGS_329inds <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/gwas_results_final/gwas_xbir-COAC-16-VIII-22-M_v2023.1.fa_FGS_329inds.txtindivs.allindiv.allscaff.vcf.summary_LR1lessthan0.01",sep="\t", header=TRUE)
gwas_xbir_WT_FGS_142inds <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/gwas_results_final/gwas_xbir-COAC-16-VIII-22-M_v2023.1.fa_FGS_142inds.txtindivs.allindiv.allscaff.vcf.summary_LR1lessthan0.01",sep="\t", header=TRUE)
gwas_xbir_WT_FGS_329inds_chr2 <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/gwas_results_final/gwas_xbir-COAC-16-VIII-22-M_v2023.1.fa_FGS_329inds.txtindivs.allindiv.chr-02.vcf.summary",sep="\t", header=TRUE)

prelim_annot_kitlg_2 <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/kitlga_manual_annotation.gtf",sep="\t", header=FALSE)
kitlg_blast <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/kitlga_exon_blast_results_plotting.txt",sep="\t", header=TRUE)

annotation <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/gwas_results_final/xbir-COAC-16-VIII-22-M_v2023.1.gff_chr-02.genes_CDS.txt",sep="\t", header=FALSE)


gwasRes <- gwas_xbir_WT_FGS_329inds

gwasRes$AF.group1 <- format(gwasRes$AF.group1, scientific = FALSE)
gwasRes$AF.group2 <- format(gwasRes$AF.group2, scientific = FALSE)

gwasRes_chr2 <- gwas_xbir_WT_FGS_329inds_chr2
gwasRes_chr2$AF.group1 <- format(gwasRes_chr2$AF.group1, scientific = FALSE)
gwasRes_chr2$AF.group2 <- format(gwasRes_chr2$AF.group2, scientific = FALSE)

#mean_gwas_depth <- mean(gwasRes$depth)
mean_gwas_depth <- 113

gwasRes <- subset(gwasRes, group != "chr-21-Y")
gwasRes$group <- gsub("-X", "", gwasRes$group)
gwasRes$chr_num <- as.numeric(
  gsub("chr-", "", gwasRes$group))

sigThresh <- -log10(3.7E-9)

df <- subset(gwasRes, (depth < 2*mean_gwas_depth & depth > 0.5*mean_gwas_depth) & (MAF >= 0.05) & LR1 < 0.01)
gwasRes_chr2 <- subset(gwasRes_chr2, (depth < 2*mean_gwas_depth & depth > 0.5*mean_gwas_depth) & MAF >= 0.05)

range(subset(df,group == "chr-02" & -log10(LR1) >= -log10(3.7E-9))$pos)
mean(c(28341222,28356777)) #28349000

subset(subset(df,group == "chr-02" & -log10(LR1) >= -log10(3.7E-9) & pos > 28348950 & pos < 28349050))
subset(subset(df,group == "chr-02"  & pos > 28348950 & pos < 28349050))

don <- df %>% 
  
  # Compute chromosome size
  dplyr::group_by(chr_num) %>% 
  dplyr::summarise(chr_len=max(pos)+6000000) %>% 
  
  # Calculate cumulative position of each chromosome
  dplyr::mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  dplyr::left_join(df, ., by=c("chr_num"="chr_num")) %>%
  
  # Add a cumulative position of each SNP
  dplyr::arrange(chr_num, pos) %>%
  dplyr::mutate(pos2=pos+tot)



axisdf = don %>% dplyr::group_by(chr_num) %>% dplyr::summarize(center=( max(pos2) + min(pos2) ) / 2 )

gwas_FGS <- ggplot(subset(don), aes(x=pos2, y=-log10(LR1))) +
  
  geom_point( aes(color=as.factor(chr_num)), alpha=0.7, size=0.9) +
  scale_color_manual(values = rep(c("grey","#55185D"), 24)) +
  scale_x_continuous( label = subset(axisdf, axisdf$chr_num %% 2 != 1)$chr_num, breaks = subset(axisdf, axisdf$chr_num %% 2 != 1)$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(min(-log10(don$LR1)), 25)) +     # remove space between plot area and x axis
  
  geom_hline(yintercept=sigThresh, size = 0.5, lty = "dashed", alpha=0.4) +
  #ggtitle(title) +
  # Custom the theme:
  theme_classic(base_size = 7) +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent")
  ) +
  xlab(expression(paste(italic("X. birchmanni "),"chromosome",sep=""))) +
  ylab(expression(paste("-log"[10],"(",italic(""*P*""),"-value)",sep=""))) #+

gwas_FGS_chr2_300kb <- ggplot() +
  geom_point(data=subset(gwasRes_chr2, group=="chr-02"), aes(x=pos, y=-log10(LR1), color=as.factor(group)), alpha=0.7, size=0.9) +
  scale_color_manual(values = rep(c("#55185D","grey"), 25)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +
  scale_x_continuous(expand = c(0, 0), limits = c(28344924 - 3e5, 28344924 + 3e5), breaks = seq(28e6, 29e6, 1e5), labels = seq(28, 29, 0.1)) +     # remove space between plot area and x axis
  geom_hline(yintercept=sigThresh, size = 0.5, lty = "dashed", alpha=0.4) +
  theme_classic(base_size = 7) +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent")
  ) +
  xlab(expression(paste(italic(""),"Chromosome 2 position (Mb)",sep=""))) +
  ylab(expression(paste("-log"[10],"(",italic(""*P*""),"-value)",sep=""))) +
  #geom_segment(data=subset(annotation, V3 == "gene" & (V5 < 28.35e6 | V4 > 28.45e6) ), aes(x=V4, xend=V5, y=23, yend=23), col="black") +
  #geom_rect(data=subset(annotation, V3 == "CDS" & (V5 < 28.35e6 | V4 > 28.45e6) ), aes(xmin=V4, xmax=V5, ymin=22.5, ymax=23.5), fill = "#424949", col="#424949") +
  #geom_segment(data=subset(prelim_annot_kitlg_2, V3 == "transcript"), aes(x=V4, xend=V5, y=23, yend=23), col="black") +
  #geom_rect(data=subset(prelim_annot_kitlg_2, V3 == "exon"), aes(xmin=V4, xmax=V5, ymin=22.5, ymax=23.5), fill = "#424949", col="#424949")

  geom_segment(data=subset(annotation, V3 == "gene" & (V5 < 28.35e6 | V4 > 28.45e6) ), aes(x=V4, xend=V5, y=24, yend=24), col="darkgrey",lwd=0.5) +
  geom_segment(data=subset(annotation, V3 == "CDS" & (V5 < 28.35e6 | V4 > 28.45e6) ), aes(x=V4, xend=V5, y=24, yend=24), col="darkgrey",lwd=2) +
  geom_segment(data=subset(prelim_annot_kitlg_2, V3 == "transcript"), aes(x=V4, xend=V5, y=24, yend=24), col="black", lwd=0.5) +
  geom_segment(data=subset(prelim_annot_kitlg_2, V3 == "exon"), aes(x=V4, xend=V5, y=24, yend=24), col = "black",lwd=2)



figd_data <- subset(gwas_xbir_WT_FGS_142inds, group == "chr-02" & pos == 28349032)
figd_data$AF.group1 <- as.numeric(figd_data$AF.group1)
figd_data$AF.group2 <- as.numeric(figd_data$AF.group2)

figd_data <- data.frame(
  phenotype = c("non-FGS","FGS"),
  AF = c(figd_data$AF.group1, figd_data$AF.group2),    # Values for the first column
  n = c(86, 56)  # Values for the second column
)
figd_data$phenotype <- factor(figd_data$phenotype, levels = c("non-FGS", "FGS"))

figd_data$se <- sqrt(figd_data$AF * (1 - figd_data$AF) / figd_data$n)

fig1d <- ggplot(figd_data) +
  geom_errorbar(aes(x=phenotype, ymin = AF - 2 * se, ymax = AF + 2 * se), col=c("darkgrey","#55185D"), width = 0.2) +
  geom_point(aes(x = phenotype, y = AF), col=c("darkgrey","#55185D"), size=2) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete() + theme_bw(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Phenotype") + ylab(paste0("Alt. AF at central SNP"))
        
fig1a <- ggplot() + theme_void()
fig1c <- ggplot() + theme_void()

fig1ac <- plot_grid(fig1a,fig1c, ncol = 2, labels = c("A","C"), rel_widths = c(3,1),label_size = 10)

fig1ef <- plot_grid(gwas_FGS_chr2_300kb,fig1d, ncol = 2, labels = c("E","F"), rel_widths = c(13,5),label_size = 10)
fig1def <- plot_grid(gwas_FGS, fig1ef, nrow = 2, labels = c("D",""),label_size = 10)

fig1 <- plot_grid(fig1ac, fig1def, nrow = 2, labels = c("A",""), rel_heights = c(0.8,1), label_size = 10)

#ggsave("/Users/trisdodge/Desktop/Stanford/falsegravid/figs/fig1.pdf",
#       fig1, width = 114, height = 140, units = "mm")





##################
#####FIGURE 2#####
##################

####figure 2####
gwas_xbir_WT_FGS_329inds_chr2 <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/gwas_results_final/gwas_xbir-COAC-16-VIII-22-M_v2023.1.fa_FGS_329inds.txtindivs.allindiv.chr-02.vcf.summary",sep="\t", header=TRUE)
plink_LD <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/gwas_results_final/COAC_2016-2018_23_allhighcovinds_plink-LD_chr-02_28349032.ld.trimmed.ld",sep="\t", skip=1, header=FALSE)

kitlg_blast <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/kitlga_exon_blast_results_plotting.txt",sep="\t", header=TRUE)

f = "/Users/trisdodge/Desktop/Stanford/falsegravid/alignment/xbir_WT_2_xbir_FGS_mummer.m.coords"


f.column.names <- c("start1","end1","start2","end2","length1","length2","perc_similar","name1","name2")
infile.coords<-read.csv(file=f, sep="\t", skip=4, head=FALSE, col.names = f.column.names)
infile.coords$segment_ID <- row(infile.coords)
infile.coords$for_or_rev <- ifelse((infile.coords$start2 - infile.coords$end2 > 0) == TRUE, "reverse", "forward")
infile.coords <- subset(infile.coords, infile.coords$name1 == "chr_02")

infile.coords <- subset(infile.coords, infile.coords$name2 == "chr_02_FGS")

subset(infile.coords, start1 > 28.285e6 & end1 < 28.415e6 & length1)
subset(infile.coords, start1 > 28.34e6 & end1 < 28.36e6 & length1)

min_length <- 0

xbir_WT <- c(28.285e6,28.415e6)

subset(infile.coords, start1 > xbir_WT[1] & end1 < xbir_WT[2])
#xbir_FGS <- c(4043680,4194237)
sigThresh <- -log10(3.7E-9)



WT_peak_region <- c(min(subset(gwas_xbir_WT_FGS_329inds_chr2, -log10(LR1) >=sigThresh)$pos), max(subset(gwas_xbir_WT_FGS_329inds_chr2, -log10(LR1) >=sigThresh)$pos))
WT_peak_region[2] - WT_peak_region[1]

head(plink_LD)

test <- merge(gwas_xbir_WT_FGS_329inds_chr2, plink_LD[,c("V5","V7")], by.x="pos", by.y="V5")

#barcoding <- subset(clean, hap == "chr-02" | hap == "chr-02-FGS")

palette2 <- c("grey","lightblue","darkgreen","red","orange")

fig2a <- 
  ggplot() + 
  geom_hline(yintercept=sigThresh*1800 + 28060607, size = 0.5, lty = "dashed", alpha=0.4) +
  geom_segment(data=subset(infile.coords, length1 > min_length & length2 > min_length), aes(x = start1, y = start2, xend = end1, yend = end2), colour = "grey", lineend = "butt",lwd = 1) +
  #labs(color = paste0("LD (",expression (~R^2),")")) +
  labs(color = expression (~R^2)) +
  #scale_fill_manual(values=palette2) +
  
  
  ylab(expression(paste(italic("X. birchmanni"), " FGS haplotype position (Mb)",sep=" "))) +
  xlab(expression(paste(italic("X. birchmanni"), " non-FGS haplotype position (Mb)",sep=" "))) +
  theme_bw(base_size = 7) +
  theme(
        axis.text = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.8,0.4),
        legend.position = c(0.2,0.8),
        legend.key.size=unit(3, 'mm')
        ) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(27.875e6,28.35e6), breaks = seq(27.81e6,28.51e6,by=.2e5), labels = seq(27.81,28.51,by=0.02),
                     sec.axis = sec_axis(~(.- 28060607)/1800, breaks = c(5,10,15,20,25), 
                                         name = expression(paste("-log"[10],"(",italic(""*P*""),"-value)",sep="")))) +
  scale_x_continuous(expand = c(0, 0), limits = c(28.075e6,28.56e6), breaks = seq(28.01e6,28.71e6,by=.2e5), labels = seq(28.01,28.71,by=0.02))+#,
  coord_fixed(xlim=c(28288640,28411220), ylim=c(28060607,28205283)) +
  
  geom_segment(data=subset(kitlg_blast, feature == "transcript"), aes(y=FGS_start, yend=FGS_end, x=WT_start, xend=WT_end), col="black", lwd = .75) +
  geom_segment(data=subset(kitlg_blast, feature == "exon"), aes(y=FGS_start, yend=FGS_end, x=WT_start, xend=WT_end), lwd = 3, col="#424949") +
  geom_point(data = subset(test, (depth < 2*113 & depth > 0.5*113) & (MAF >= 0.05)), aes(x = pos, y =-log10(LR1) * 1800 + 28060607, color=V7), size=0.9) + 
  scale_colour_distiller(palette="RdYlBu")
  

####fig2b####

library(dplyr)
library(ggplot2)
library(gggenes)
library(ggh4x)
library(forcats)
palette <- c("grey","lightblue","darkgreen","red","orange4","orange")

clean <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/pangenome/barcoding/barcoding.fa_all_coords_clean", sep="\t", header = TRUE)
info <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/pangenome/barcoding/barcode_information.txt", sep="\t", header = TRUE)
clean <- unique(clean)

clean <- subset(clean, abs(start - stop) > 1000)

clean_group <- clean %>%
  dplyr::filter(barcode == "syntenic2") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(syntenic2_end = max(start, stop)) %>%
  dplyr::select(hap, syntenic2_end)

clean <- merge(clean, clean_group, by="hap")

clean$pos <- clean$syntenic2_end - 60000
clean$start2 <- clean$start - clean$pos
clean$stop2 <- clean$stop - clean$pos

clean <- merge(clean, info, by="individual")
clean$pop_ID <- paste(clean$pop, clean$ID, sep="_")

clean <- clean %>% 
  mutate(barcode = recode(barcode, syntenic1 = 'co-linear', syntenic2 = 'co-linear', insertion = 'INS (FGS)', wt_insertion = 'INS (non-FGS)', inversion = 'INV', dup1 = "SD1", dup2 = "SD2"))

clean$barcode <- fct_relevel(clean$barcode, "co-linear", "SD1", "SD2", "INS (FGS)", "INS (non-FGS)")

clean$pop <- fct_relevel(clean$pop, "ref", "COAC", "BEJU", "IZAP", "xmal")
clean$ID <- fct_relevel(clean$ID, "C1","F1")
clean$haplotype <- fct_relevel(clean$haplotype, "FGS (ref)", "non-FGS (ref)", "h1", "h2")

#write.table(clean, "/Users/trisdodge/Desktop/Stanford/falsegravid/pangenome/barcoding/barcoding_alldata.txt", quote=FALSE, sep="\t", row.names=FALSE)

clean_subset <- subset(clean, 
       ID=="F1" & allele == "FGS" | 
         ID=="C2" & allele == "FGS" | 
         ID=="C3" & allele == "FGS" | 
         ID=="B1" & allele == "FGS" | 
         ID=="I1" & allele == "FGS" |
         ID=="C1" & reference == "ref" | 
         ID=="C4" & allele == "non-FGS"
)


clean_subset <- clean_subset %>% 
  mutate(allele_class = recode(allele_class, 'FGS (1)' = 'FGS_ref', 'FGS (5)' = 'FGS_3', 'FGS (2)' = 'FGS_2', 'FGS (4)' = 'FGS_4', 'FGS (3)' = 'FGS_5', 'non-FGS (1)' = 'non-FGS_ref','non-FGS (2)' = 'non-FGS_2'))
clean_subset$allele_class <- fct_relevel(clean_subset$allele_class, "FGS_ref", "FGS_2", "FGS_3", "FGS_4", "FGS_5", "non-FGS_ref", "non-FGS_2")
clean_subset <- clean_subset %>% 
  mutate(n = recode(allele_class, 'FGS_ref' = '(2)', 'FGS_2' = '(2)', 'FGS_3' = '(2)', 'FGS_4' = '(1)', 'FGS_5' = '(1)', 'non-FGS_ref' = '(17)','non-FGS_2' = '(1)'))

fig2c <- 
  ggplot(clean_subset, aes(xmin = start2, xmax = stop2, y = allele_class, fill = barcode)) +
  geom_gene_arrow(data = clean_subset, arrow_body_height = unit (1.25, "mm"), arrowhead_height = unit(2.5, "mm"), arrowhead_width = unit(1.25, "mm")) +
  geom_text(data = subset(clean_subset, reference == "ref"), aes(x=-4000, label = allele_class), size=2, hjust=0) +
  geom_text(aes(x=61300, label = n), size=2, hjust=0.5) +
  #geom_point(aes(x=-1000,y=ID, color=phenotype)) +
  #facet_wrap(~allele) +
  facet_nested(
    rows = vars(allele),
    nest_line = TRUE,
    scales = "free", 
    space = 'free',
    switch = "y") +
  scale_fill_manual(values=palette) +
  scale_x_continuous(breaks = seq(0e4,6e4,1e4), labels=seq(60,0,-10)) +
  scale_y_discrete(limits=rev) +
  coord_cartesian(xlim=c(-2000,60000),clip = 'off') +
  theme_minimal(base_size = 7) +
  labs(x="position (kb)") +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        #strip.text.y.left = element_text(angle = 0),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size=unit(3, 'mm')) +
  guides(fill = guide_legend(nrow = 1))



#tree
library(ape)
library(treeio)
library(tidytree)
library(ggtree)

BiocManager::install("treeio")
inv <- ape::read.tree("/Users/trisdodge/Desktop/Stanford/falsegravid/phylo/RAxML_bestTree.clustal_phylohaps_inversion_combined_SNPsites.fa")
#p<-
  ggtree(inv) + 
    geom_tree() + geom_tiplab(offset = .6, hjust = .5)

#fig2

library(cowplot)
fig2b <- ggplot() + theme_void()
fig2d <- ggplot() + theme_void()
fig2e <- ggplot() + theme_void()
fig2bc <- plot_grid(fig2b, fig2c, nrow = 2, labels = c("B","C"), rel_heights= c(1.5,1),label_size = 10)
fig2abc <- plot_grid(fig2a, fig2bc, ncol = 2, labels = c("A",""), rel_widths= c(1,1.2),label_size = 10)
fig2de <- plot_grid(fig2d, fig2e, ncol = 2, labels = c("D","E"), rel_widths= c(1,1.5),label_size = 10)
fig2 <- plot_grid(fig2abc, fig2de, nrow = 2, labels = c("",""), rel_heights = c(1.5,1),label_size = 10)

#fig2
ggsave("/Users/trisdodge/Desktop/Stanford/falsegravid/figs/fig2.pdf",
       fig2, width = 174, height = 150, units = "mm")

  
##################
#####FIGURE 3#####
##################

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

ASE <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/2021.04.29_xbir_xmal_FGS_ASE.csv", header=TRUE)
ASE <- subset(ASE, ind != "oligo_validation" & quality == "passed")
ASE <- subset(ASE, locus == "kitlg759" & tissue != "testes") #FLAG

ASE$rel_exp <- ASE$xbir_allele/ASE$xmal_allele

ASE$tissue <- gsub("flank", "PM + EAM", ASE$tissue )
ASE$tissue <- as.factor(ASE$tissue)
ASE$locus <- as.factor(ASE$locus)
ASE$tissue <- factor(ASE$tissue, levels = c("PM + EAM", "brain", "gill"))


ASE_means <- ASE  %>%
  group_by(tissue) %>%
  summarise(
    rel_exp_mean = mean(rel_exp),
    rel_exp_se = sd(rel_exp) / sqrt(n())
  )


fig3d <- ggplot() +
  geom_point(data=ASE, aes(x=tissue, y=rel_exp), color="#55185D", alpha=0.2,
             position = position_jitter(width=0.1)) +
  geom_errorbar(data=ASE_means, aes(x=tissue, ymin = rel_exp_mean - 2 * rel_exp_se, ymax = rel_exp_mean + 2 * rel_exp_se), col="#55185D", width = 0.2) +
  geom_point(data=ASE_means, aes(x = tissue, y = rel_exp_mean), col="#55185D", size=2) +
  geom_hline(yintercept=1, size = 0.5, lty = "dashed", alpha=0.4) +
  scale_y_continuous(limits=c(0,3)) +
  xlab("tissue") +
  ylab("ratio FGS:non-FGS allele") +
  theme_bw(base_size = 7) + theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())


####fig 3e####

male_dev_prelim_follow[male_dev_prelim_follow$focal_1_5 != "focal",c("Tag")]
male_dev_prelim_follow[male_dev_prelim_follow$Tag=="BBG",c()]

male_dev_prelim_follow <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/male_dev/plotting_data_unique_final.csv", header=TRUE)
male_dev_prelim_follow$gon_stage1_round_date <- ymd(male_dev_prelim_follow$gon_stage1_round_date)
male_dev_prelim_follow$fgs_second_minus1week_date <- ymd(male_dev_prelim_follow$fgs_second_minus1week_date)
male_dev_prelim_follow$gon_stage5_round_date <- ymd(male_dev_prelim_follow$gon_stage5_round_date)

male_dev_prelim_follow$days_1tofgs2 <- as.numeric(days(male_dev_prelim_follow$fgs_second_minus1week_date - male_dev_prelim_follow$gon_stage1_round_date))/86400
male_dev_prelim_follow$days_1to5round <- as.numeric(days(male_dev_prelim_follow$gon_stage5_round_date-male_dev_prelim_follow$gon_stage1_round_date))/86400
male_dev_prelim_follow$days_fgs2to5round <- as.numeric(days(male_dev_prelim_follow$gon_stage5_round_date - male_dev_prelim_follow$fgs_second_minus1week_date))/86400

male_dev_prelim_follow_focal <- subset(male_dev_prelim_follow, focal_1_5 == "focal")

male_dev_prelim_follow_focal_means <- male_dev_prelim_follow_focal %>%
  group_by(FGS_adult_phenotype) %>%
  summarise(
    days_1to5round_mean = mean(days_1to5round),
    days_1tofgs2_mean = mean(days_1tofgs2),
    days_fgs2to5round_mean= mean(days_fgs2to5round)
  )

fig3e <- ggplot(data = male_dev_prelim_follow_focal) + 
  geom_segment(data=subset(male_dev_prelim_follow_focal, FGS_adult_phenotype == "FGS"), aes(x =0, xend = days_1tofgs2, y=Sort_focal+2, yend=Sort_focal+2), color = "#55185D", alpha=0.4) +
  geom_segment(data=subset(male_dev_prelim_follow_focal, FGS_adult_phenotype == "FGS"), aes(x =days_1tofgs2, xend = days_1tofgs2 + days_fgs2to5round, y=Sort_focal+2, yend=Sort_focal+2), color = "#55185D", alpha=1) +
  geom_segment(data=subset(male_dev_prelim_follow_focal, FGS_adult_phenotype == "non-FGS"), aes(x =0, xend = days_1to5round, y=Sort_focal+2, yend=Sort_focal+2), color = "grey", alpha=1) +
  geom_point(data=male_dev_prelim_follow_focal, aes(x=days_1tofgs2, y=Sort_focal+2), color="#55185D", alpha=1) +
  geom_segment(data=subset(male_dev_prelim_follow_focal_means, FGS_adult_phenotype == "FGS"), aes(x=0, xend = days_1tofgs2_mean, y=1, yend=1), color = "#55185D", alpha=0.4, size=2) +
  geom_segment(data=subset(male_dev_prelim_follow_focal_means, FGS_adult_phenotype == "FGS"), aes(x =days_1tofgs2_mean, xend = days_1tofgs2_mean + days_fgs2to5round_mean, y=1, yend=1), color = "#55185D", alpha=1, size=2) +
  geom_segment(data=subset(male_dev_prelim_follow_focal_means, FGS_adult_phenotype == "non-FGS"), aes(x =0, xend = days_1to5round_mean, y=2, yend=2), color = "grey", alpha=1, size=2) +
  geom_point(data=male_dev_prelim_follow_focal_means, aes(x=days_1tofgs2_mean, y=1), color="#55185D", alpha=1, size=3) +
  xlab("Days after puberty onset") +
  ylab("Individual") +
  theme_bw(base_size = 7) + theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())




####fig3f####
#write.table(male_dev_prelim_subset, "/Users/trisdodge/Desktop/Stanford/falsegravid/male_dev/male_dev_final_subset.txt", quote=FALSE, sep="\t", row.names=FALSE)
male_dev_prelim <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/male_dev/male_dev_final_fulldataset.txt", header=TRUE, sep="\t")

length(unique(male_dev_prelim$Tag))
male_dev_prelim_subset <- base::subset(male_dev_prelim, Dev_stage == "gon_stage1" | Dev_stage == "gon_stage5" | Dev_stage == "fgs_second_minus1week" | Dev_stage == "date_done")
male_dev_prelim_subset$Dev_stage <- factor(male_dev_prelim$Dev_stage, levels = c("gon_stage1", "fgs_second_minus1week", "gon_stage5", "date_done"))

fig3f <- 
  ggplot() + #paper
  geom_point(data=male_dev_prelim_subset, aes(x=Dev_stage, y=ldr_sl, color = FGS_adult_phenotype, group = Tag), size=2, alpha=0.2) + 
  geom_line(data=male_dev_prelim_subset, aes(x=Dev_stage, y=ldr_sl, color = FGS_adult_phenotype, group = Tag), alpha=0.2) +
  geom_boxplot(data=male_dev_prelim_subset, aes(x=Dev_stage, y=ldr_sl, color=FGS_adult_phenotype), position = position_dodge(width = 0.4), width=0.2, fill = NA, outlier.shape=NA) +
  theme_bw(base_size = 7) + theme(legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Developmental stage") +
  ylab("Standardized dorsal length") +
  scale_x_discrete(labels=c("puberty onset", "FGS present", "mature", "last measurement")) +
  scale_color_manual(values = c("#55185D","darkgrey"))



subset(male_dev_prelim_subset, Dev_stage == "gon_stage5")[, 7:16]
dev_numeric <- subset(male_dev_prelim, Dev_stage == "gon_stage5" & Tag != "YYG")[, c(7:13,15:16)]
scaled_dev <- scale(dev_numeric)
pca_result <- prcomp(scaled_dev)
pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_data <- cbind(pca_scores, subset(male_dev_prelim_subset, Dev_stage == "gon_stage5" & Tag != "YYG")[,1:2])
summary(pca_result)

fig3g <- ggplot(pca_data, aes(y = PC1, x = PC2, color = FGS_adult_phenotype)) +
  geom_point(size = 2, alpha=0.9) +
  ylab("PC1 (53.77%)") +
  xlab("PC2 (19.21%)") +
  theme_bw(base_size = 7) + theme(legend.position = "none",
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(-4,4,2)) +
  scale_y_continuous(breaks = seq(-4,4,2)) +

  scale_color_manual(values = c("#55185D","darkgrey")) +
  coord_fixed()

#new

#new
fig3ab <- plot_grid(fig3a,fig3b, ncol = 2, labels = c("A","B"), rel_widths = c(1,1),label_size = 10)
fig3abc <- plot_grid(fig3ab, fig3c, ncol = 2, labels = c("","C"), rel_widths = c(1,0.85),label_size = 10)
fig3abcd <- plot_grid(fig3abc, fig3d, ncol = 2, labels = c("","D"), rel_widths = c(4.3,1),label_size = 10)
fig3ef <- plot_grid(fig3e, fig3f, ncol = 2, labels = c("E","F"), rel_widths = c(1,1.4),label_size = 10)
fig3efg <- plot_grid(fig3ef, fig3g, ncol = 2, labels = c("","G"), rel_widths = c(3.6,1),label_size = 10)

fig3 <- plot_grid(fig3abcd,fig3efg, nrow = 2, labels = c("",""), rel_heights = c(1,1),label_size = 10)

#fig3
ggsave("/Users/trisdodge/Desktop/Stanford/falsegravid/figs/fig3.pdf",
       fig3, width = 174, height = 100, units = "mm")


##################
#####FIGURE 5#####
##################

####figure 5####
#map, PSMC, sims, tajD, AF change + selection, COAC freq

####map
#install.packages("elevatr")
library(elevatr)
#install.packages("rgeoboundaries")
library(rgeoboundaries)
library(sf)
library(raster)
library(ggplot2)
library(viridis)
#install.packages("scatterpie")
library(scatterpie)
library(ggrepel)
library(ggspatial)
#install.packages("devtools")
library(ggnewscale)

#devtools::install_github('oswaldosantos/ggsn')

library(rnaturalearth)
library(remotes)
#remotes::install_github("ropensci/rnaturalearthhires")
#library(rnaturalearthhires)

freq <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/FGS_frequency_all_populations.txt", sep="\t",header = TRUE)
#subset(freq, sample >= 100)$FGS_freq
sf_mexico <- ne_states(country = "mexico", returnclass = "sf")
sf_mexico <- sf_mexico[sf_mexico$postal %in% c("HG","SL","VE","PU","TM","QT"), ]

elevation_1 <- elevatr::get_elev_raster(locations = sf_mexico, z = 9, clip = "locations") #change to z=8 for final plot
cropped_elev <- crop(elevation_1, sf_mexico)
elevate <- as.data.frame(cropped_elev, xy = TRUE)

colnames(elevate)[3] <- "elevation_value"
elevate <- elevate[complete.cases(elevate), ]

xbir_range <- sf::st_read("/Users/trisdodge/Desktop/xbir_range.kml")
rios <- sf::st_read("/Users/trisdodge/Desktop/RIOS_FGS.kml")

fig5a <- 
  ggplot() +
  geom_raster(data = elevate, aes(x = x, y = y, fill = elevation_value)) +
  #scale_fill_gradientn(colours = terrain.colors(10)) +
  scale_fill_gradient(low = "white", high = "#333333", limits=c(-200,3100)) +
  labs(x = "Longitude", y = "Latitude",fill = "Elevation (m)") +
  new_scale("fill") +
  geom_sf(data = sf_mexico, color = "white", fill = NA) +
  geom_sf(data=xbir_range, alpha=0.4, fill="darkred", color="darkred") +
  geom_sf(data=rios, color="lightblue3") +
  geom_label_repel(data=freq,aes(x=long, y=lat, label=population_n, fill=nanopore),
                   size = 2,
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.size  = 0.4,
                   segment.color = "black") +
  scale_fill_manual(values=c("#00A087B2", "#7E6148B2","#8491B4B2","#FFFFFF")) + 
  new_scale("fill") +
  geom_scatterpie(data=freq,aes(x=long, y=lat, r=0.026), cols = c("FGS","nonFGS")) +
  coord_sf(xlim = c(-97.61,-99.29), ylim = c(20.54, 21.44)) +
  scale_fill_manual(values=c("#55185D","grey")) +
  annotation_scale(location="tr", height = unit(1, "mm")) +
  labs(fill = "Phenotype") +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.y = element_text(angle = 90, hjust=0.5))
  
####psmc
PSMC_data_PASS <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/nanopore/psmc_analysis/xbir_allpops_psmc_combined_primary_bootstrap100.txt", header=FALSE, sep = "\t")
#PSMC_data_PASS$V6 <- gsub(" ", "", PSMC_data_PASS$V6)
#PSMC_data_PASS$V7 <- str_extract(PSMC_data_PASS$V6, "(?<=xbir-)\\w{4}")
PSMC_data_PASS$unique <- paste(PSMC_data_PASS$V9, PSMC_data_PASS$V7, PSMC_data_PASS$V8,sep="_")


fig5b <- ggplot() +
  #geom_step(data=subset(PSMC_data_PASS, V7=="primary" & ((V7=="IZAP" & V1 > 2e2) | (V1 > 5e2))), aes(x=V1*2, y=V2*10, group=V9, col=V6)) +
  #geom_step(data=subset(PSMC_data_PASS, V7=="bootstrap"& ((V7=="IZAP" & V1 > 2e2) | (V1 > 5e2))), aes(x=V1*2, y=V2*10, group=unique, col=V6), alpha=0.05) +
  geom_step(data=subset(PSMC_data_PASS, V7=="bootstrap"), aes(x=V1*2, y=V2*10, group=unique, col=V6), alpha=0.05) +
  geom_step(data=subset(PSMC_data_PASS, V7=="primary"), aes(x=V1*2, y=V2*10, group=V9, col=V6)) +
  
  scale_x_log10( #set limits, breaks, labels however you please
    expand=c(0,0),
    breaks = c(1e3,1e4,1e5),
    labels = c("1,000", "10,000","100,000")) +
  
  scale_y_continuous(
    expand=c(0,0), breaks = c(seq(0,150,25)) #set limits, breaks, labels however you please
  ) +
  #geom_vline(xintercept =c(5e2,2e2))+
  theme_classic(base_size = 7) +
  theme(legend.key.size=unit(2.5, 'mm'),
        legend.position = c(.45,.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Generations in past") +
  ylab(expression(paste("N"[e], "*", 10^{3},sep=" "))) + #can change scaling so 10^4, 10^6 etc, just also modify the scaling factor in geom_step(). psmc plot automatically outputs 10^4
  labs(color='Population') +
  annotation_logticks(side = "b") +
  coord_cartesian(xlim=c(5e2, 1.2e5), ylim = c(0,75)) +
  scale_color_manual(values=c("#00A087B2", "#7E6148B2","#8491B4B2"))

fig5c <- ggplot() +
  geom_step(data=subset(PSMC_data_PASS, V6=="IZAP" & V7=="primary"), aes(x=V1*2, y=V2*1e4, group=V9, col=V6)) +
  geom_step(data=subset(PSMC_data_PASS, V6=="IZAP" & V7=="bootstrap"), aes(x=V1*2, y=V2*1e4, group=unique, col=V6), alpha=0.05) +
  
  scale_x_continuous( #set limits, breaks, labels however you please
    breaks = seq(0,5000,1000),
    labels = seq(0,5000,1000)) +
  
  scale_y_continuous(
    expand=c(0,0), breaks = c(seq(0,2500,500)) #set limits, breaks, labels however you please
  ) +
  #geom_vline(xintercept =c(5e2,2e2))+
  theme_classic(base_size = 7) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+#,
  xlab("Generations in past") +
  ylab(expression(paste("N"[e],sep=" "))) + #can change scaling so 10^4, 10^6 etc, just also modify the scaling factor in geom_step(). psmc plot automatically outputs 10^4
  labs(color='Population') +
  scale_color_manual(values=c("#8491B4B2")) + 
  coord_cartesian(xlim=c(0, 3.5e3), ylim = c(0,2500) )

library(cowplot)

####sims figure####
izap_sim <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/sims/IZAPM03_results.txt", header=FALSE, sep=" ")
izap_sim_end <- subset(izap_sim, V3 == 64725)
izap_sim_end$pop <- "IZAP"
coac_sim <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/sims/COACM03_short_results.txt", header=FALSE, sep=" ")
coac_sim_end <- subset(coac_sim, V3 == 63049)
coac_sim_end$pop <- "COAC"
beju_sim <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/sims/BEJUM01_short_results.txt", header=FALSE, sep=" ")
beju_sim_end <- subset(beju_sim, V3 == 72729)
beju_sim_end$pop <- "BEJU"


sims <- rbind(rbind(izap_sim_end, coac_sim_end),beju_sim_end)


sims$cut <- cut(sims$V2, breaks = c(-Inf, 0, 1-(1e-10), Inf), labels = c("lost", "polymorphic", "fixed"), include.lowest = TRUE)
sims$dummy <- 1
sims$cut <- factor(sims$cut, levels = c("fixed", "polymorphic", "lost"))

sims$cut2 <- cut(sims$V2, breaks = c(-Inf, 0, 1-(1e-10), Inf), labels = c("mono.", "poly.", "mono."), include.lowest = TRUE)
sims$cut2 <- factor(sims$cut2, levels = c("mono.", "poly."))

sims_summary <- sims%>%
  dplyr::group_by(pop) %>%
  dplyr::summarise(poly_prop = as.numeric(round(mean(cut2 == 'poly.'), 3)))

fig5d <- 
  ggplot() +
  geom_bar(data= subset(sims, pop != "coac_long"), aes(x = as.factor(pop), y = dummy/10000, fill = cut2), position="stack", stat = "identity", alpha=1) +
  labs(x = "Population", y = "Prop. neutral sims", fill = ""
       ) +
  geom_text(data=sims_summary, aes(x = as.factor(pop), y = 1.05, label = as.character(poly_prop)), size=2) +
  scale_y_continuous(breaks = c(seq(0,1,.25)))+
  geom_hline(yintercept=0.05, size = 0.5, lty = "dashed") +
  theme_classic(base_size = 7) +
  scale_fill_manual(values=c("grey","#55185D")) +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size=unit(1, 'mm'))


pi_8kb <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/balancing/COAC_2016-2018_23_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf_pixy_pi_8794_pi.txt", sep="\t", head=TRUE)

#mean(na.omit(pi_8kb$avg_pi))
fig5e <- ggplot(pi_8kb) +
  geom_density(aes(x=avg_pi), fill="grey", alpha=0.2, color="#7E6148B2") +
  labs(x = "Diversity (π)", y = "Density") +
  geom_vline(xintercept=0.007264126380375, size = 1, alpha=0.9, color="#55185D") +
  annotate("text", x = 0.007264126380375, y = 900, hjust=-0.1, size = 2, label = paste0(round(nrow(subset(pi_8kb, avg_pi > 0.007264126380375)) / nrow(pi_8kb)*100,2),"%")) +
  theme_classic(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0,0.01), breaks=seq(0,0.012,0.003))

nrow(subset(pi_8kb, avg_pi > 0.007264126380375)) / nrow(pi_8kb)

taj8.794kb <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/balancing/chr-02.tajD.Tajima.D", sep="\t", head=TRUE)


fig5f <- 
  ggplot(taj8.794kb) +
  geom_density(aes(x=TajimaD), fill="grey", alpha=0.2, color="#7E6148B2") +
  labs(x = "Tajima's D", y = "Density") +
  geom_vline(xintercept=2.03781, size = 1, alpha=0.9, color="#55185D") +
  annotate("text", x = 2.03781, y = 0.3, hjust=-0.1, size = 2, label = paste0(round(nrow(subset(taj8.794kb, TajimaD > 2.03781)) / nrow(na.omit(taj8.794kb))*100,2),"%")) +
  theme_classic(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

nrow(subset(taj8.794kb, TajimaD > 2.03781)) / nrow(na.omit(taj8.794kb))



####ncd
# install.packages("devtools")
library(devtools)
#devtools::install_github("bitarellolab/balselr")
library(balselr)
library(stringr)
library(ggplot2)

#install.packages("tictoc")

#infile_sys <- read_vcf(x=system.file(package="balselr", "example.vcf"))

ncd1_input_chr2_focal <- parse_vcf(infile="/Users/trisdodge/Desktop/Stanford/falsegravid/balancing/ncd1_preprocessed_focal.vcf", n0=23, type="ncd1")
ncd1_input_chr2_focal_static <- ncd1_input_chr2_focal

ncd1_chr2_8kb_inversion <- ncd1(x = ncd1_input_chr2_focal_static,
                      tf = 0.3,
                      w = 8794,
                      ncores = 2,
                      minIS = 8)
ncd_dataframe_inversion <- ncd1_chr2_8kb_inversion

ncd_dataframe_inversion$start <- as.integer(str_extract(ncd_dataframe_inversion$Win.ID, "(?<=_)\\d+(?=_)"))
ncd_dataframe_inversion$stop <- as.integer(str_extract(ncd_dataframe_inversion$Win.ID, "(?<=_)\\d+$"))

ncd_dataframe_inversion

ncd1_input_chr2 <- parse_vcf(infile="/Users/trisdodge/Desktop/Stanford/falsegravid/balancing/ncd1_preprocessed.vcf", n0=23, type="ncd1")
ncd1_input_chr2_static <- ncd1_input_chr2

ncd1_chr2_8kb <- ncd1(x = ncd1_input_chr2_static,
     tf = 0.3,
     w = 8794,
     ncores = 2,
     minIS = 8)
ncd_dataframe <- ncd1_chr2_8kb

ncd_dataframe$start <- as.integer(str_extract(ncd_dataframe$Win.ID, "(?<=_)\\d+(?=_)"))
ncd_dataframe$stop <- as.integer(str_extract(ncd_dataframe$Win.ID, "(?<=_)\\d+$"))


fig5g <- ggplot(ncd_dataframe[seq(1, nrow(ncd_dataframe), 2), ]) +
  geom_density(aes(x=ncd1), fill="grey", alpha=0.2, color="#7E6148B2") +
  labs(x = "NCD1", y = "Density") +
  geom_vline(xintercept=0.09552293, size = 1, alpha=0.9, color="#55185D") +
  annotate("text", x = 0.09552293, y = 10, hjust=-0.1, size = 2, label = paste0(round(nrow(subset(ncd_dataframe, ncd1 < 0.09552293)) / nrow(na.omit(ncd_dataframe))*100,2),"%")) +
  
  theme_classic(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

nrow(subset(ncd_dataframe, ncd1 < 0.09552293)) / nrow(na.omit(ncd_dataframe))

fig5g <- ggplot + theme_void()
####fig 5####
library(cowplot)
fig5bc <- plot_grid(fig5b, fig5c, ncol = 2, labels = c("B","C"), rel_widths = c(1,1),label_size = 10)
fig5bcd <- plot_grid(fig5bc, fig5d, ncol = 2, labels = c("","D"), rel_widths = c(2,1),label_size = 10)
fig5ef <- plot_grid(fig5e, fig5f, ncol = 2, labels = c("E","F"), rel_widths = c(1,1),label_size = 10)
fig5efg <- plot_grid(fig5ef, fig5g, ncol = 2, labels = c("","G"), rel_widths = c(2,1),label_size = 10)

fig5bcdefg <- plot_grid(fig5bcd, fig5efg, nrow = 2, labels = c("",""), rel_heights = c(1,1),label_size = 10)

fig5 <- plot_grid(fig5a, fig5bcdefg, nrow = 2, labels = c("A",""), rel_heights = c(1,1.1),label_size = 10)

#fig5


#fig5
ggsave("/Users/trisdodge/Desktop/Stanford/falsegravid/figs/fig5.pdf",
       fig5, width = 114, height = 150, units = "mm")



##################
#####FIGURE 6#####
##################

####fig6####


fig6b <- ggplot() +
  geom_point(data=subset(male_agg_means_mean, target != "female"), aes(x=target, y=chase_tot_d, col=target), size=3) + 
  geom_errorbar(data=subset(male_agg_means_mean, target != "female"), aes(x=target, ymin=chase_tot_d-2*se, ymax=chase_tot_d+2*se, col=target), width=.1) + 
  geom_jitter(data=subset(male_agg_means, target != "female"), aes(x=target,y=chase_tot_d, col=target), height=0, width = 0.1, alpha=0.2) +
  ylab("Dominant male chases (sec.)") +
  xlab("Stimulus male phenotype") +
  theme(legend.position = "NA") +
  scale_color_manual(values = c("#55185D","darkgrey")) +
  theme_bw(base_size = 7) +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))


fig6c <- 
  ggplot() +
  geom_point(data=subset(male_court_means_mean, target != "female"), aes(x=target, y=court_d, col=target), size=3) + 
  geom_errorbar(data=subset(male_court_means_mean, target != "female"), aes(x=target, ymin=court_d-2*se, ymax=court_d+2*se, col=target), width=.1) + 
  #geom_line(data=subset(male_agg_means, target != "female"), aes(x=target, y=court_d, group=fm_id), alpha=0.1) + 
  
  geom_jitter(data=subset(male_agg_means, target != "female"), aes(x=target,y=court_d, col=target), height=0, width = 0.1, alpha=0.2) +
    ylab("Dominant male courtship (sec.)") +
  xlab("Stimulus male phenotype") +
  theme(legend.position = "NA") +
  scale_color_manual(values = c("#55185D","darkgrey")) +
  theme_bw(base_size=7) +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))

#FGSvsWT_mean_long_summary <- FGSvsWT_mean_long_summary %>% mutate(phenotype = case_when(phenotype == "FGS animation" ~ "FGS", phenotype == "non-FGS animation" ~ "non-FGS", TRUE ~ as.character(phenotype)))
#FGSvsWT_mean_long <- FGSvsWT_mean_long %>% mutate(phenotype = case_when(phenotype == "FGS animation" ~ "FGS", phenotype == "non-FGS animation" ~ "non-FGS", TRUE ~ as.character(phenotype)))

field_obs <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/behavior/birchmanni_field_obs_forSPSS.csv", head=TRUE)

males_assoc_means <- Rmisc::summarySE(na.omit(field_obs[,c("X0","Phenotype","Males_associated")]), measurevar="Males_associated", groupvars=c("Phenotype"))

fig6g <- 
  ggplot() +
  geom_point(data=males_assoc_means, aes(x=Phenotype, y=Males_associated, col=Phenotype), size=3) + 
  geom_errorbar(data=males_assoc_means, aes(x=Phenotype, ymin=Males_associated-2*se, ymax=Males_associated+2*se, col=Phenotype), width=.1) + 
  geom_jitter(data=subset(field_obs), aes(x=Phenotype,y=Males_associated, col=Phenotype), width = 0.1,height=0, alpha=0.1) +
  ylab("Males associated (in wild)") +
  xlab("Phenotype") +
  theme(legend.position = "NA") +
  scale_color_manual(values = c("#55185D","darkgrey")) +
  theme_bw(base_size=7) +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))

SLDR <- read.csv("~/Downloads/xbir_COAC_dorsal_measurements_FGS_added.csv", header=TRUE)

SLDR <- SLDR[,c("label","dorsallength","raylength","mean_SL","phenotype","confidence")]
SLDR <- subset(SLDR, confidence == "include")
#SLDR <- na.omit(SLDR)
SLDR$phenotype <- as.factor(SLDR$phenotype)
SLDR$raylength_SL <- SLDR$raylength/SLDR$mean_SL
nrow(SLDR)

SLDR$year <- as.factor(str_extract(SLDR$label, "(?<=COAC-\\w{1,2}-)(17|18)"))

SLDR_lm <- lm(data=SLDR, raylength_SL ~ phenotype + year)
t.test(raylength_SL ~ phenotype, data=SLDR)
#SLDR_lm <- lm(data=SLDR, mean_SL ~ FGS + year)
#SLDR_lm <- lm(data=SLDR, raylength_SL ~ FGS + year)

summary(SLDR_lm)

SLDR_means <- summarySE(na.omit(SLDR), measurevar="raylength_SL", groupvars=c("phenotype"))

fig6h <- 
  ggplot() +
  geom_jitter(data=subset(SLDR), aes(x=phenotype,y=raylength_SL, col=phenotype), width = 0.1, height=0, alpha=0.1) +
  geom_point(data=SLDR_means, aes(x=phenotype, y=raylength_SL, col=phenotype), size=3) + 
  geom_errorbar(data=SLDR_means, aes(x=phenotype, ymin=raylength_SL-2*se, ymax=raylength_SL+2*se, col=phenotype), width=.1) + 
  xlab("Phenotype") +
  ylab("Standardized dorsal length (in wild)") +
  scale_color_manual(values = c("#55185D","darkgrey")) +
  theme_bw(base_size=7) +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_)) #+
#stat_compare_means(data=subset(SLDR), aes(x=phenotype, y= raylength_SL, group=phenotype), label.x = 1.5, label = "p.format", method="t.test")


###fig 6
library(cowplot)
fig6a <- ggplot() + theme_void()
fig6d <- ggplot() + theme_void()
fig6f <- ggplot() + theme_void()

fig6bc <- plot_grid(fig6b, fig6c, ncol = 2, labels = c("B","C"), rel_widths = c(1,1),label_size = 10)
fig6abc <- plot_grid(fig6a, fig6bc, ncol = 2, labels = c("A",""), rel_widths = c(1,1.3),label_size = 10)
fig6de <- plot_grid(fig6d, fig6e, ncol = 2, labels = c("D","E"), rel_widths = c(1,1.15),label_size = 10)

fig6gh <- plot_grid(fig6g, fig6h, ncol = 2, labels = c("G","H"), rel_widths = c(1,1),label_size = 10)
fig6fgh <- plot_grid(fig6f, fig6gh, ncol = 2, labels = c("F",""), rel_widths = c(1,1),label_size = 10)

fig6abcde <- plot_grid(fig6abc, fig6de, nrow = 2, labels = c("",""), rel_heights = c(1,1.1),label_size = 10)
fig6 <- plot_grid(fig6abcde, fig6fgh, nrow = 2, labels = c("",""), rel_heights = c(2,1),label_size = 10)


#fig6
ggsave("/Users/trisdodge/Desktop/Stanford/falsegravid/figs/fig6.pdf",
       fig6, width = 174, height = 174, units = "mm")







############
####SUPP####
############

####sup pi across genome

prelim_annot_kitlg_2 <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/kitlga_manual_annotation.gtf",sep="\t", header=FALSE)
kitlg_blast <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/kitlga_exon_blast_results_plotting.txt",sep="\t", header=TRUE)

annotation <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/gwas_results_final/xbir-COAC-16-VIII-22-M_v2023.1.gff_chr-02.genes_CDS.txt",sep="\t", header=FALSE)
WT_peak <- c(28341258,28356027)
pi3kb <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/balancing/COAC_2016-2018_27_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf_pixy_pi_3000_pi.txt", sep="\t", head=TRUE)
ggplot(pi3kb) +
  geom_step(aes(x=window_pos_1, y=avg_pi), lwd=0.4, color="#7E6148B2") +
  theme_classic(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(28344924 - 3e5, 28344924 + 3e5), breaks = seq(28e6, 29e6, 1e5), labels = seq(28, 29, 0.1))   +  # remove space between plot area and x axis
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.012)) +
  xlab(expression(paste(italic(""),"Chromosome 2 position (Mb)",sep=""))) +
  ylab("pi") +
  geom_segment(data=subset(annotation, V3 == "gene" & (V5 < 28.35e6 | V4 > 28.45e6) ), aes(x=V4, xend=V5, y=0.011, yend=0.011), col="darkgrey",lwd=0.5) +
  geom_segment(data=subset(annotation, V3 == "CDS" & (V5 < 28.35e6 | V4 > 28.45e6) ), aes(x=V4, xend=V5, y=0.011, yend=0.011), col="darkgrey",lwd=2) +
  geom_segment(data=subset(prelim_annot_kitlg_2, V3 == "transcript"), aes(x=V4, xend=V5, y=0.011, yend=0.011), col="black", lwd=0.5) +
  geom_segment(data=subset(prelim_annot_kitlg_2, V3 == "exon"), aes(x=V4, xend=V5, y=0.011, yend=0.011), col = "black",lwd=2) +
  annotate("rect", xmin=WT_peak[1], xmax=WT_peak[2], ymin=0, ymax=0.012, alpha=0.2, fill="#55185D") +
  coord_cartesian(xlim=c(28344924 - 3e5,28344924 + 3e5))


####ncd stat across chrom2
ggplot(ncd_dataframe) +
  #geom_segment(aes(x=start, xend=stop, y=ncd1, yend=ncd1), lwd=0.4) +
  geom_line(aes(x=start, y=ncd1), lwd=0.4, color="#7E6148B2") +
  
  theme_classic(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(28344924 - 3e5, 28344924 + 3e5), breaks = seq(28e6, 29e6, 1e5), labels = seq(28, 29, 0.1))   +  # remove space between plot area and x axis
  scale_y_continuous(expand = c(0, 0), limits = c(0.05, 0.27)) +
  theme_classic(base_size = 7) +
  xlab(expression(paste(italic(""),"Chromosome 2 position (Mb)",sep=""))) +
  ylab("ncd1") +
  annotate("rect", xmin=WT_peak[1], xmax=WT_peak[2], ymin=0.05, ymax=0.27, alpha=0.2, fill="#55185D") +
  coord_cartesian(xlim=c(28344924 - 3e5,28344924 + 3e5))


#female pref data
ggplot() +
  geom_point(data=FGSvsWT_mean_long_summary, aes(x=phenotype, y=zone_time, col=phenotype), size=3) + 
  geom_errorbar(data=FGSvsWT_mean_long_summary, aes(x=phenotype, ymin=zone_time-2*se, ymax=zone_time+2*se,col=phenotype), width=.1) + 
  geom_jitter(data=subset(FGSvsWT_mean_long), aes(x=phenotype,y=zone_time, col=phenotype), width = 0.1, height=0, alpha=0.5) +
  ylab("Female association time (prop)") +
  xlab("Stimulus male animation phenotype") +
  scale_y_continuous(limits=c(0,0.6),breaks = c(0,0.2,0.4,0.6)) +
  theme(legend.position = "NA") +
  scale_color_manual(values = c("#55185D","darkgrey")) +
  theme_bw(base_size = 7) +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))



####SUPPORTING FIGURES#####

#Fig S7
#GWAS to chr-02-FGS hap

gwas_xbir_FGS_FGS_329inds_chr2 <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/gwas_FGS/gwas_xbir-COAC-16-VIII-22-M_v2023.1_FGS_hap.fa_FGS_329inds.txtindivs.allindiv.chr-02-FGS.vcf.summary",sep="\t", header=TRUE)
kitlg_blast <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/kitlga_exon_blast_results_plotting.txt",sep="\t", header=TRUE)

mean_gwas_depth <- 113
sigThresh <- -log10(3.7E-9)


ggplot() +
  geom_point(data=subset(gwas_xbir_FGS_FGS_329inds_chr2, group=="chr-02-FGS" & depth <= 2*mean_gwas_depth & depth >= 0.5*mean_gwas_depth & MAF >= 0.05), aes(x=pos, y=-log10(LR1), color=as.factor(group)), alpha=0.7, size=0.9) +
  scale_color_manual(values = rep(c("#55185D","grey"), 25)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +
  scale_x_continuous(expand = c(0, 0), limits = c(28153355 - 3e5, 28153355 + 3e5), breaks = seq(27e6, 29e6, 1e5), labels = seq(27, 29, 0.1)) +     # remove space between plot area and x axis
  geom_hline(yintercept=sigThresh, size = 0.5, lty = "dashed", alpha=0.4) +
  theme_classic(base_size = 10) +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent")
  ) +
  xlab(expression(paste(italic(""),"Chromosome 2 position (Mb)",sep=""))) +
  ylab(expression(paste("-log"[10],"(",italic(""*P*""),"-value)",sep=""))) +
  geom_segment(data=subset(kitlg_blast, species == "birchmanni" & feature == "transcript"), aes(x=FGS_start, xend=FGS_end, y=24, yend=24), col="black", lwd=0.5) +
  geom_segment(data=subset(kitlg_blast, species == "birchmanni" & feature == "exon"), aes(x=FGS_start, xend=FGS_end, y=24, yend=24), col = "black",lwd=2)



#FIG S8

coord <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/mummer/xmal-CHIC-XI-20-M.asm.bp.hap1.p_ctg_h1tg000087l.fa_2_xbir-COAC-16-VIII-22-M_v2023.1_chr-02.fa.delta.m.coords", sep="\t", skip=4, header=FALSE, 
                  col.names = c("start1","end1","start2","end2","length1","length2","perc_similar","name1","name2"))

M1 <- ggplot(data = coord, aes(x = start1, y = start2, xend = end1, yend = end2)) +
  geom_segment(lineend = "butt", lwd = 1) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(color="Percent identity") +
  ylab(expression(paste(italic("X. malinche"), " haplotype 1 position (Mb)",sep=" "))) +
  xlab(expression(paste(italic("X. birchmanni"), " non-FGS haplotype position (Mb)",sep=" "))) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 35e6, by = .1e5), labels = seq(0, 35, by = 0.01)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 35e6, by = .1e5), labels = seq(0, 35, by = 0.01)) +
  coord_fixed(xlim = c(28332530, 28366966), 
              ylim = c(0+7700, 34436+7700))

coord <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/mummer/xmal-CHIC-XI-20-M.asm.bp.hap2.p_ctg_h2tg000026l.fa_2_xbir-COAC-16-VIII-22-M_v2023.1_chr-02.fa.delta.m.coords", sep="\t", skip=4, header=FALSE, 
                  col.names = c("start1","end1","start2","end2","length1","length2","perc_similar","name1","name2"))

M2 <- ggplot(data = coord, aes(x = start1, y = start2, xend = end1, yend = end2)) +
  geom_segment(lineend = "butt", lwd = 1) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(color="Percent identity") +
    ylab(expression(paste(italic("X. malinche"), " haplotype 2 position (Mb)",sep=" "))) +
    xlab(expression(paste(italic("X. birchmanni"), " non-FGS haplotype position (Mb)",sep=" "))) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 35e6, by = .1e5), labels = seq(0, 35, by = 0.01)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 35e6, by = .1e5), labels = seq(0, 35, by = 0.01)) +
  coord_fixed(xlim = c(28332530, 28366966), 
              ylim = c(0+28774250, 34436+28774250)
              )
  
  plot_grid(M1, M2,
    ncol = 2)

#FIG S9
  coord <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/mummer/xbir-COAC-16-VIII-22-M_v2023.1_FGS_hap_chr-02-FGS.fa_2_xbir-COAC-16-VIII-22-M_v2023.1_chr-02.fa.delta.m.coords", sep="\t", skip=4, header=FALSE, 
                    col.names = c("start1","end1","start2","end2","length1","length2","perc_similar","name1","name2"))
  
  ggplot(data = coord, aes(x = start1, y = start2, xend = end1, yend = end2, colour = perc_similar)) +
    geom_segment(lineend = "butt", lwd = 1) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(color="Percent identity") +
    ylab(expression(paste(italic("X. birchmanni"), " FGS haplotype position (Mb)",sep=" "))) +
    xlab(expression(paste(italic("X. birchmanni"), " non-FGS haplotype position (Mb)",sep=" "))) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 35e6, by = .1e5), labels = seq(0, 35, by = 0.01)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 35e6, by = .1e5), labels = seq(0, 35, by = 0.01)) +
    coord_fixed(xlim = c(28332530, 28366966), ylim = c(28104475, 28161091)) +
    scale_color_distiller(palette = "Spectral", 
    limits = c(93, 100))

  #Fig S10
  
  library(dplyr)
  library(ggplot2)
  library(gggenes)
  library(ggh4x)
  library(forcats)
  palette <- c("grey","lightblue","darkgreen","red","orange4","orange")
  
  clean <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/pangenome/barcoding/barcoding.fa_all_coords_clean", sep="\t", header = TRUE)
  info <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/pangenome/barcoding/barcode_information.txt", sep="\t", header = TRUE)
  clean <- unique(clean)
  
  clean <- subset(clean, abs(start - stop) > 1000)
  
  clean_group <- clean %>%
    dplyr::filter(barcode == "syntenic2") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(syntenic2_end = max(start, stop)) %>%
    dplyr::select(hap, syntenic2_end)
  
  clean <- merge(clean, clean_group, by="hap")
  
  clean$pos <- clean$syntenic2_end - 60000
  clean$start2 <- clean$start - clean$pos
  clean$stop2 <- clean$stop - clean$pos
  
  clean <- merge(clean, info, by="individual")
  clean$pop_ID <- paste(clean$pop, clean$ID, sep="_")
  
  clean <- clean %>% 
    mutate(barcode = recode(barcode, syntenic1 = 'co-linear', syntenic2 = 'co-linear', insertion = 'INS (FGS)', wt_insertion = 'INS (non-FGS)', inversion = 'INV', dup1 = "SD1", dup2 = "SD2"))
  
  clean$barcode <- fct_relevel(clean$barcode, "co-linear", "SD1", "SD2", "INS (FGS)", "INS (non-FGS)")
  
  clean$pop <- fct_relevel(clean$pop, "ref", "COAC", "BEJU", "IZAP", "xmal")
  clean$ID <- fct_relevel(clean$ID, "C1","F1")
  clean$haplotype <- fct_relevel(clean$haplotype, "FGS (ref)", "non-FGS (ref)", "h1", "h2")
  
  
  ggplot(clean, aes(xmin = start2, xmax = stop2, y = haplotype, fill = barcode)) +
    geom_gene_arrow(arrow_body_height = unit (1, "mm"), arrowhead_height = unit(2, "mm"), arrowhead_width = unit(1, "mm")) +
    geom_text(aes(x=61000, label = allele_class), size=1.5, hjust=0) +
    geom_text(aes(x=-3000, label = haplotype), size=1.5, hjust=0) +
    #geom_point(aes(x=-1000,y=ID, color=phenotype)) +
    facet_nested(
      rows = vars(pop,ID),
      nest_line = TRUE,
      scales = "free", 
      space = 'free',
      switch = "y") +
    scale_fill_manual(values=palette) +
    scale_x_continuous(breaks = seq(0e4,6e4,1e4), labels=seq(60,0,-10)) +
    scale_y_discrete(limits=rev) +
    coord_cartesian(xlim=c(0,65000),clip = 'off') +
    theme_minimal(base_size = 7) +
    labs(x="position (kb)") +
    theme(axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          strip.text.y.left = element_text(angle = 0),
          legend.position = "top",
          legend.title = element_blank(),
          legend.key.size=unit(3, 'mm')) +
    guides(fill = guide_legend(nrow = 1)) +
    #annotate("point", y = 1.5, x = -4000, size =3, color="purple")
    annotate("segment", y = 0.5, yend = 2.5, x = -4200, xend = -4200, size =2, color="black")
  
  #Fig S17
  
  COAC_freq <- read.csv("/Users/trisdodge/Downloads/FGS_COAC_population_frequency.csv", head=TRUE)
  COAC_freq$n <- COAC_freq$fgs + COAC_freq$wt
  
  COAC_freq$se <- sqrt(COAC_freq$freq * (1 - COAC_freq$freq) / COAC_freq$n)
  
  
  ggplot(COAC_freq) +
    geom_line(aes(x = Date, y = freq), alpha=0.5)+
    geom_errorbar(aes(x=Date, ymin = freq - 2 * se, ymax = freq + 2 * se), alpha=0.4, width = 0.1) +
    geom_point(aes(x = Date, y = freq, size=n, alpha=n), col="#55185D") +
    scale_y_continuous(limits=c(0,1)) +
    #scale_x_discrete() + 
    theme_bw(base_size = 10) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("Sampling date") + ylab("FGS phenotype frequency") + labs(size="sample (n)", alpha="sample (n)") +
    scale_alpha(range = c(0.3, 1))
  
  
  
  ####IR
  IR <- read.csv("/Users/trisdodge/Desktop/Stanford/falsegravid/pangenome/barcoding/non-beta/gfa_xbir_wt_FGS_chr-02.positions.txt", sep="\t", header = FALSE)
  
  ggplot() +
    geom_rect(data=subset(IR, V1=="chr-02-FGS" & V2 =="Inverted_Repeat"), aes(xmin=V3, xmax=V4, ymin=1.1, ymax=1.2, fill=V2)) +
    geom_gene_arrow(data=subset(clean, ID=="F1" & haplotype == "FGS (ref)"), aes(xmin = start, xmax = stop, y = 1, fill = barcode), arrow_body_height = unit (1.25, "mm"), arrowhead_height = unit(2.5, "mm"), arrowhead_width = unit(1.25, "mm")) +
    theme_minimal(base_size = 7) +
    coord_cartesian(xlim=c(28105000,28160000))
  
  ggplot() +
    geom_rect(data=subset(IR, V1=="chr-02" & V2 =="Inverted_Repeat"), aes(xmin=V3, xmax=V4, ymin=1.1, ymax=1.2, fill=V2)) +
    geom_gene_arrow(data=subset(clean, ID=="C1" & haplotype == "non-FGS (ref)"), aes(xmin = start, xmax = stop, y = 1, fill = barcode), arrow_body_height = unit (1.25, "mm"), arrowhead_height = unit(2.5, "mm"), arrowhead_width = unit(1.25, "mm")) +
    theme_minimal(base_size = 7) +
    coord_cartesian(xlim=c(28332000,28368000))
  
  
  #Fig S17
  
  library(pheatmap)
  
  file_names <- Sys.glob("/Users/trisdodge/Desktop/Stanford/falsegravid/nanopore/popgen/*_combined_polydiv")
  
  summary_tibble <- tibble(poly1=numeric(), poly2=numeric(),
                           dxy=numeric(), homo=numeric(), ind1=character(), ind2=character())
  
  
  for (file in file_names){
    polydiv_data <- NULL
    polydiv_data <- read.csv(file, header=FALSE, sep = "\t")
    
    polydiv_data$poly1_scaled <- polydiv_data$V4/polydiv_data$V8
    polydiv_data$poly2_scaled <- polydiv_data$V5/polydiv_data$V8
    polydiv_data$dxy_scaled <- polydiv_data$V6/polydiv_data$V8
    polydiv_data$homo_scaled <- polydiv_data$V7/polydiv_data$V8
    
    polydiv_poly1_wm <- weighted.mean(polydiv_data$poly1_scaled, polydiv_data$V8)
    polydiv_poly2_wm <- weighted.mean(polydiv_data$poly2_scaled, polydiv_data$V8)
    polydiv_dxy_wm <- weighted.mean(polydiv_data$dxy_scaled, polydiv_data$V8)
    polydiv_homo_wm <- weighted.mean(polydiv_data$homo_scaled, polydiv_data$V8)
    
    matches <- str_extract_all(file, "xbir-(?!COAC-16-VIII-22-M_v2023\\.[0-9]+)[^.]+")
    
    match1 <- unlist(matches)[1]
    match2 <- unlist(matches)[2]
    
    summary_tibble <- add_row(summary_tibble, poly1=polydiv_poly1_wm, poly2=polydiv_poly2_wm,
                              dxy=polydiv_dxy_wm, homo=polydiv_homo_wm, ind1=match1, ind2=match2)
  }

  tmp <- as.data.frame(cbind(summary_tibble$ind1, summary_tibble$ind2, as.numeric(summary_tibble$dxy)))

  tmp$V1 <- as.character(tmp$V1)
  tmp$V2 <- as.character(tmp$V2)
  tmp$V3 <- as.numeric(tmp$V3)
  
  wide_matrix <- as.matrix(pivot_wider(tmp, names_from = V2, values_from = V3)[, -1])
  rownames(wide_matrix) <- NULL
  
  pheatmap(wide_matrix, display_numbers=TRUE, number_format="%.5f")

  
  #Fig S21
  
  ###################
  ####power curve####
  ###################
  
  
  data<-read.csv(file="~/Desktop/Stanford/falsegravid/behavior/female_preference_proportions_power_curve.csv") #need dataframe with all data (experienced, not experienced)
  
  data <- data %>%
    mutate(ID = paste(str_extract(origin, "[A-Za-z]+"), str_extract(comparision, "[A-Za-z]+"),sep="_"))
  
  result_df <- data.frame(comparison = character(),
                          effect_size = numeric(),
                          p_value = numeric())
  prop_diff_df <- data.frame(comparison = character(),
                             effect_size = numeric())
  effect_sizes <- seq(0,0.4,0.025)
  
  for (entry in unique(data$ID)) {
    focal <- na.omit(subset(data, ID == entry)) #subset dataframe
    sample_size <- length(focal[, 1]) #count samples
    mean_prop_difference <- mean(focal$prop_diff)
    sd_prop_difference <- sd(focal$prop_diff) #compute SD
    
    for (effect_size in effect_sizes) {
      pvals_prop_trial <- numeric() #initialize p-values
      
      #simulation
      for (x in 1:10000) {
        diff_rep_prop <- rnorm(sample_size, mean = effect_size, sd = sd_prop_difference)
        pvals_prop_trial <- c(pvals_prop_trial, wilcox.test(diff_rep_prop)$p.value)
      }
      
      # Combine p-values with comparison group and add to result dataframe
      temp_df <- data.frame(comparison = rep(entry, 10000),
                            effect_size = rep(effect_size, 10000), 
                            p_value = pvals_prop_trial)
      result_df <- rbind(result_df, temp_df)
      
    }
    temp_pd_df <- data.frame(comparison = entry,
                             effect_size = mean_prop_difference)
    prop_diff_df <- rbind(prop_diff_df, temp_pd_df)
  }
  
  data_summarized <- result_df %>%
    dplyr::group_by(comparison, effect_size) %>%
    dplyr::summarize(power = length(p_value[p_value<0.05])/length(p_value)) %>%
    dplyr::mutate(comparison = recode(comparison, 'lab_ornamented'='lab born, ornamented', 'lab_unornamented'='lab born, unornamented', 'wild_ornamented'='wild-caught, ornamented', 'wild_unornamented'='wild-caught, unornamented'))
  
  prop_diff_df <- prop_diff_df %>%
    dplyr::mutate(comparison = recode(comparison, 'lab_ornamented'='lab born, ornamented', 'lab_unornamented'='lab born, unornamented', 'wild_ornamented'='wild-caught, ornamented', 'wild_unornamented'='wild-caught, unornamented'))
  
  ggplot() +
    geom_line(data=data_summarized, aes(x=effect_size,y=power*100,col=comparison), lty="dashed") +
    geom_vline(data=prop_diff_df, aes(xintercept=abs(effect_size), col=comparison)) +
    xlab("Effect size (absolute value)") +
    ylab("Power (%)") +
    theme_bw(base_size = 10)
  
  