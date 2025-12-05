

# use R v4.4.4.1 (not earlier ones, which will generate errors) 
# 
# 1 - Setup, grouping & data import
# 2 - PCA plots
# 3 - get biomart data
# 4 - limma exploration
# 5 - Sleuth across genes initial

########## 1 Setup, grouping & data import ###################

if (!require("BiocManager", quietly = T)) { 
    install.packages("BiocManager") }

#install.packages("aggregation") # for genes from tx
library(aggregation)
#install.packages("devtools")
library(devtools)
#BiocManager::install("rhdf5")
library(rhdf5)
#devtools::install_github("pachterlab/sleuth")
library(sleuth)
#vignette('intro', package = 'sleuth')
library(ggplot2) # v3.5.2
library(ggrepel)
library(dplyr)
library(ggpubr)
library(grid) # v4.4.1
library(gridExtra)
#BiocManager::install("limma")
library(limma) # v3.60.6
#BiocManager::install("edgeR")
library(edgeR)
#BiocManager::install("tximport")
library(tximport)
library(readr)
library(ggpubr)

# [1] Load a table describing our sample, conditions
# and the source directories such that the 1st column
# contains the sample names, the middle column(s)
# contain the conditions, and the last column has the
# folder containing the kallisto output
s2c <- read.csv("table.csv", header=T, sep="\t")
str(s2c) # check table 

# reallocate table into groups: 
s2c_c <- subset(s2c, condition=="C") # 
s2c_wt <- subset(s2c, condition=="WT") # 
s2c_4x <- subset(s2c, condition=="4x") # 

# Now, we will use sleuth_prep() to make an object 
# with our experiment info, model and groups.
colnames(s2c)[1] = "sample"
so <- sleuth_prep(s2c, extra_bootstrap_summary=T, read_bootstrap_tpm=T)
# 35,793 targets pass filter

######## 2 PCA plots ##################################

# Get TPM values, filter, make PCA plots

tpm_data <- so$obs_norm %>% # extract TPM values
  dplyr::select(target_id, sample, tpm) %>%
  tidyr::pivot_wider(names_from = sample, values_from = tpm)
# Convert to matrix format for PCA
tpm_matrix <- as.matrix(tpm_data[,-1])
rownames(tpm_matrix) <- tpm_data$target_id
var_per_gene <- apply(tpm_matrix, 1, var) # filter for non-zero values
tpm_matrix_filtered <- tpm_matrix[var_per_gene > 0, ]
pca_result <- prcomp(t(tpm_matrix_filtered), scale = T)
pca_df <- as.data.frame(pca_result$x)
pca_df$sample <- rownames(pca_df)
pca_df <- left_join(pca_df, so$sample_to_covariates, by = "sample")
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
unique_combinations <- unique(pca_df$condition)

# define the base colors
base_colors <- c("WT" = "darkred", "4x"= "navy", "C" = "black" )
color_palette <- base_colors[unique_combinations]

library(stringr)
pdf("PC1.PC2.pdf", width=6.7, height=3.5)
ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) +
  geom_point(size =5, alpha=0.5) +
  geom_label(
  aes(label = str_replace_all(sample, c("_S\\d" = "", "VEEV_" = "", "_RNA_" = ""))),
  size =4, alpha =.6,  position = position_jitter(width =5, height =15)) + 
    scale_color_manual(values = color_palette, name = "") +
  labs( x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")  ) +
  theme_bw() + theme( axis.text = element_text(size = 12),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5)  )
dev.off()

pdf("PC3.PC4.pdf", width=6.7, height=3.5)
ggplot(pca_df, aes(x = PC3, y = PC4, color = condition)) +
  geom_point(size =5, alpha=0.5) +
  geom_label(
  aes(label = str_replace_all(sample, c("_S\\d" = "", "VEEV_" = "", "_RNA_" = ""))),
  size =4, alpha =.6,  position = position_jitter(width =5, height =15)) + 
  scale_color_manual(values = color_palette, name = "") +
  labs( x = paste0("PC3 (", var_explained[3], "%)"),
    y = paste0("PC4 (", var_explained[4], "%)")  ) +
  theme_bw() + theme( axis.text = element_text(size = 12),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5)  )
dev.off()

#### visualise samples as heatmap

pdf("plot_sample_heatmap.pdf", width=6, height=6)
plot_sample_heatmap(so ,use_filtered = T, color_high = "white",
  color_low = "dodgerblue", x_axis_angle = 50,
  annotation_cols = setdiff(colnames(so$sample_to_covariates), "sample"),
  cluster_bool = T)
dev.off()

####### 3 get biomart data ##########################

# Let's test across genes  
# We need to map the isoforms to genes with BiomaRt
# Note we need to have annotation matching our reference cDNAs

# BiocManager::install("biomaRt")
library(biomaRt)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = 'https://jan2024.archive.ensembl.org')
str(mart) # check it worked # 398
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
      "transcript_version", "ensembl_gene_id", "external_gene_name",
            "description", "transcript_biotype"), mart = mart)
str(t2g) # 149132 rows of:
        # ensembl_transcript_id
        # ensembl_gene_id
        # external_gene_name
# now we'll rename the transcripts
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
t2g$target_id <- paste(t2g$target_id, t2g$transcript_version,sep=".") # Re-name
subset(t2g, ext_gene=="Npw")#  check for example

############ 4 limma exploration ########################

# transcripts first

#  Examine the TPM data using tximport  
files = paste(s2c$path, "abundance.h5", sep="/")

# Import Kallisto abundance.h5 files with tximport
txi.kallisto <- tximport(files, type = "kallisto", txOut = T)
str(txi.kallisto)
y <- DGEList(txi.kallisto$counts)
dim(y) # 255407

full <- read.csv("table.csv", header=T, sep="\t") # metadata
# filtering using the design information:
full$condition <- factor(full$condition, levels = c("C", "WT", "4x"))
design <- model.matrix(~condition, data = full)
colnames(design)
contrasts <- makeContrasts(
  WT_vs_C  = conditionWT,  FourX_vs_C = condition4x,
  FourX_vs_WT = condition4x - conditionWT, levels = design )

keep <- filterByExpr(y, design)
y2 <- y[keep, ]
str(y2) # 36,038 
y2 <- calcNormFactors(y2) # normalize and run voom transformation
v <- voom(y2, design) # v is now ready for lmFit()  
fit <- lmFit(v, design) # eBayes stands for empirical Bayes
fit <-  contrasts.fit(fit, contrasts)
fitm <- eBayes(fit, trend=T)
str(fitm)
wt_c <- topTable(fitm, coef = "WT_vs_C", number = Inf)
c_4x <- topTable(fitm, coef = "FourX_vs_C", number = Inf)
wt_4x <- topTable(fitm, coef = "FourX_vs_WT", number = Inf)
str(wt_c)  
str(c_4x)  
str(wt_4x)  # 36038
write.csv(wt_c, "limma.wt_c.csv") # 37207 transcripts 
write.csv(c_4x, "limma.c_4x.csv") # 37207 transcripts 
write.csv(wt_4x, "limma.wt_4x.csv") # 37207 transcripts 

pdf("plotSA.limma.pdf")
plotSA(fitm) # we have a variance trend, so keep trend=T in eBayes()
dev.off()

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    Cor <- abs(cor(x, y)) # Remove abs function if desired
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
    if(missing(cex.cor)) {  cex.cor <- 1 + 0.4 / strwidth(txt)  }
    text(0.5, 0.5, txt, cex = 1 + cex.cor * Cor) 
    } # Resize the text by level of correlation

pdf("pairs.limma.pdf", width=6, height=6)
pairs( out1[,1:4], cex=0.1, upper.panel = panel.cor,
       lower.panel = panel.smooth, cex.labels = 1.5)
dev.off()

### make a volcano plot  for wt_c

results <- wt_c
transcript_ids <- rownames(results)
str(results) # 37027

if(grepl("\\.", transcript_ids[1])) {
  results$target_id <- rownames(results)
  results_with_genes <- left_join(
    data.frame(target_id = rownames(results), results, row.names = NULL),
    t2g %>% dplyr::select(target_id, ens_gene, ext_gene),
    by = "target_id" ) } else {
  results$ensembl_transcript_id <- rownames(results)
  results_with_genes <- left_join(
    data.frame(ensembl_transcript_id=rownames(results),results, row.names=NULL),
    t2g %>% dplyr::select(ensembl_transcript_id = target_id, ens_gene,ext_gene),
    by = "ensembl_transcript_id" ) }

results_with_genes$gene_label <- ifelse(
  !is.na(results_with_genes$ext_gene) & results_with_genes$ext_gene != "", 
  results_with_genes$ext_gene,
  ifelse(!is.na(results_with_genes$ens_gene) & results_with_genes$ens_gene !="",
    results_with_genes$ens_gene,
    ifelse(grepl("\\.", transcript_ids[1]), results_with_genes$target_id,
           results_with_genes$ensembl_transcript_id) ) )

results_with_genes$logFC = -results_with_genes$logFC 

results_with_genes$significant <- ifelse(
     (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC > +3.18), "Higher in WT",
     ifelse(
     (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC < -3.18), "Higher in Control", "Not" ) )

top_genes <- results_with_genes %>% filter(
       (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC >3.18)| 
  (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC< -3.18))  
write.csv(top_genes, "limma.DE.wt_c.csv",quote=T,row.names=F) 
length(unique(sort(top_genes$target_id))) # transcripts 2244
length(unique(sort(top_genes$ext_gene))) # total genes 1592
length(unique(sort(subset(top_genes, significant=="Higher in WT")$target_id))) ## 1227 transcripts, higher WT
length(unique(sort(subset(top_genes, significant=="Higher in WT")$ext_gene))) # 1024 genes, higher WT
length(unique(sort(subset(top_genes, significant=="Higher in Control")$target_id))) ## 1017 transcripts, higher control
length(unique(sort(subset(top_genes, significant=="Higher in Control")$ext_gene))) # 601 genes, higher control

p1 <- ggplot(results_with_genes,aes(x=logFC,y=-log10(adj.P.Val),color=significant)) + 
  geom_point(shape =1, size =1.2,  alpha=.5) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +
  geom_vline(xintercept = c(-3.18,3.18), linetype = "dashed", color = "darkgray") +
  geom_text_repel(  data = top_genes, aes(label = gene_label), box.padding =0.2,
   max.overlaps =10,   force=1, size=4, force_pull = 0.01,  point.padding = 0.3,
   		nudge_y=0.1) + xlim(-15.1,15.1) + 
  scale_color_manual(values = c("Higher in WT" = "red", "Higher in Control"="blue",
  			    "Not" = "black"), guide="none") + ylim(0,7) + 
  labs(    x = "log2(FC)", y = "-log10(p)",      title = "WT vs Control") + 
  theme_bw() + theme( legend.position = "none",
  	     plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14), axis.text=element_text(size = 14),
    legend.text = element_text(size=14), legend.title=element_text(size=14))   +
    geom_text(  x =9, y = max(-log10(wt_c$adj.P.Val)) + .2,
      label = "Higher in WT", color = "darkred", size = 5, hjust = 0) +
    geom_text(  x =-7, y = max(-log10(wt_c$adj.P.Val)) + .2,
      label = "Higher in Control", color = "navy", size = 5, hjust = 1)

pdf("volcano.wt_c.pdf", height=5, width=6.5)
print(p1)
dev.off()

### make a volcano plot  for c_4x

results <- c_4x
transcript_ids <- rownames(results)

if(grepl("\\.", transcript_ids[1])) {
  results$target_id <- rownames(results)
  results_with_genes <- left_join(
    data.frame(target_id = rownames(results), results, row.names = NULL),
    t2g %>% dplyr::select(target_id, ens_gene, ext_gene),
    by = "target_id" ) } else {
  results$ensembl_transcript_id <- rownames(results)
  results_with_genes <- left_join(
    data.frame(ensembl_transcript_id=rownames(results),results, row.names=NULL),
    t2g %>% dplyr::select(ensembl_transcript_id = target_id, ens_gene,ext_gene),
    by = "ensembl_transcript_id" ) }

results_with_genes$gene_label <- ifelse(
  !is.na(results_with_genes$ext_gene) & results_with_genes$ext_gene != "", 
  results_with_genes$ext_gene,
  ifelse(!is.na(results_with_genes$ens_gene) & results_with_genes$ens_gene !="",
    results_with_genes$ens_gene,
    ifelse(grepl("\\.", transcript_ids[1]), results_with_genes$target_id,
           results_with_genes$ensembl_transcript_id) ) )

results_with_genes$logFC = -results_with_genes$logFC # so infection ND is +ve 

results_with_genes$significant <- ifelse(
     (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC > +3.18), "Higher in 4x",
     ifelse(
     (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC < -3.18), "Higher in Control", "Not" ) )

top_genes <- results_with_genes %>% filter(
       (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC >3.18)| 
  (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC< -3.18))  
write.csv(top_genes, "limma.DE.c_4x.csv",quote=T,row.names=F) # c1

results_with_genes$direction <- ifelse(results_with_genes$logFC > 0, "Higher in 4x", "Higher in Control")

length(unique(sort(top_genes$target_id))) # transcripts 42
length(unique(sort(top_genes$ext_gene))) # total genes 37
length(unique(sort(subset(top_genes, significant=="Higher in 4x")$target_id))) ## 2 transcripts, higher WT
length(unique(sort(subset(top_genes, significant=="Higher in 4x")$ext_gene))) # 2 genes, higher WT
length(unique(sort(subset(top_genes, significant=="Higher in Control")$target_id))) ## 40 transcripts, higher control
length(unique(sort(subset(top_genes, significant=="Higher in Control")$ext_gene))) # 36 genes, higher control


p2 <- ggplot(results_with_genes,aes(x=logFC,y=-log10(adj.P.Val),color=significant)) + 
  geom_point(shape =1, size =1.2,  alpha=.5) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +
  geom_vline(xintercept = c(-3.18,3.18), linetype = "dashed", color = "darkgray") +
  geom_text_repel(  data = top_genes, aes(label = gene_label), box.padding =0.12,
   max.overlaps =10,   force=1, size=4, force_pull = 0.01,  point.padding = 0.25,
   		nudge_y=0.1) + xlim(-15.1,15.1) + 
  scale_color_manual(values = c("Higher in 4x" = "darkgreen", "Higher in Control"="blue",
  			    "Not" = "black"), guide="none") + ylim(0,7) + 
  labs(    x = "log2(FC)", y = "-log10(p)",      title = "4x vs Control") + 
  theme_bw() + theme( legend.position = "none",
  	     plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14), axis.text=element_text(size = 14),
    legend.text = element_text(size=14), legend.title=element_text(size=14))   +
    geom_text(  x =9, y = max(-log10(wt_c$adj.P.Val)) - .01,
      label = "Higher in 4x", color = "darkgreen", size = 5, hjust = 0) +
    geom_text(  x =-7, y = max(-log10(wt_c$adj.P.Val)) - .01,
      label = "Higher in Control", color = "navy", size = 5, hjust = 1)

pdf("volcano.c_4x.pdf", height=5, width=6.5)
print(p2)
dev.off()

###
results <- wt_4x
transcript_ids <- rownames(results)

if(grepl("\\.", transcript_ids[1])) {
  results$target_id <- rownames(results)
  results_with_genes <- left_join(
    data.frame(target_id = rownames(results), results, row.names = NULL),
    t2g %>% dplyr::select(target_id, ens_gene, ext_gene),
    by = "target_id" ) } else {
  results$ensembl_transcript_id <- rownames(results)
  results_with_genes <- left_join(
    data.frame(ensembl_transcript_id=rownames(results),results, row.names=NULL),
    t2g %>% dplyr::select(ensembl_transcript_id = target_id, ens_gene,ext_gene),
    by = "ensembl_transcript_id" ) }

results_with_genes$gene_label <- ifelse(
  !is.na(results_with_genes$ext_gene) & results_with_genes$ext_gene != "",
  results_with_genes$ext_gene,
  ifelse(!is.na(results_with_genes$ens_gene) & results_with_genes$ens_gene !="",
    results_with_genes$ens_gene,
    ifelse(grepl("\\.", transcript_ids[1]), results_with_genes$target_id,
           results_with_genes$ensembl_transcript_id) ) )

results_with_genes$logFC = -results_with_genes$logFC # so infection ND is +ve

results_with_genes$significant <- ifelse(
     (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC > +3.18), "Higher in 4x",
     ifelse(
     (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC < -3.18), "Higher in WT", "Not" ) )

top_genes <- results_with_genes %>% filter(
       (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC >3.18)|
  (results_with_genes$adj.P.Val < 0.05  & results_with_genes$logFC< -3.18))
write.csv(top_genes, "limma.DE.wt_4x.csv",quote=T,row.names=F) # c1

length(unique(sort(top_genes$target_id))) # transcripts 1911
length(unique(sort(top_genes$ext_gene))) # total genes 1533
length(unique(sort(subset(top_genes, significant=="Higher in 4x")$target_id))) ## 1506 transcripts, higher WT
length(unique(sort(subset(top_genes, significant=="Higher in 4x")$ext_gene))) # 1275 genes, higher WT
length(unique(sort(subset(top_genes, significant=="Higher in WT")$target_id))) ## 405 transcripts, higher control
length(unique(sort(subset(top_genes, significant=="Higher in WT")$ext_gene))) # 287 genes, higher control

p3 <- ggplot(results_with_genes,aes(x=logFC,y=-log10(adj.P.Val),color=significant)) +
  geom_point(shape =1, size =1.2,  alpha=.5) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +
  geom_vline(xintercept = c(-3.18,3.18), linetype = "dashed", color = "darkgray") +
  geom_text_repel(  data = top_genes, aes(label = gene_label), box.padding =0.2,
   max.overlaps =10,   force=1, size=4, force_pull = 0.01,  point.padding = 0.3,
                nudge_y=0.15) + xlim(-15.1,15.1) +
  scale_color_manual(values = c("Higher in 4x" = "darkgreen", "Higher in WT"="red",
                            "Not" = "black"), guide="none") + ylim(0,7) +
  labs(    x = "log2(FC)", y = "-log10(p)",      title = "4x vs WT") +
  theme_bw() + theme( legend.position = "none",
             plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14), axis.text=element_text(size = 14),
    legend.text = element_text(size=14), legend.title=element_text(size=14))   +
    geom_text(  x =9.5, y = max(-log10(wt_c$adj.P.Val)) - .005,
      label = "Higher in 4x", color = "darkgreen", size = 5, hjust = 0) +
    geom_text(  x =-7, y = max(-log10(wt_c$adj.P.Val)) - .01,
      label = "Higher in WT", color = "darkred", size = 5, hjust = 1)

pdf("volcano.wt_4x.pdf", height=5, width=6.5)
print(p3)
dev.off()

# combined
pdf("volcano_combined.pdf", height=5, width=18)
ggarrange(  p1, p2, p3, ncol = 3, nrow = 1,
  labels = c("A", "B", "C"), common.legend = F)
  dev.off()
#pdf("ILR1L1_2.pdf", height=6, width=8)
#plot_bootstrap(so2, "ENSMUST00010033183.1", x_axis_angle = 90,
#               units="est_counts", color_by="condition")
#dev.off()

#### get pariwse correlation plots of logFCs ##

wt_c2 <- read.csv("limma.wt_c.csv", header=T)
str(wt_c2)
wt_4x2 <- read.csv("limma.wt_4x.csv", header=T)
str(wt_4x2)
c_4x2 <- read.csv("limma.c_4x.csv", header=T)
str(c_4x2)

# Rename columns for clarity
wt_c2 <- wt_c2 %>% rename(transcript_id = X, logFC_wt_c = logFC)
wt_4x2 <- wt_4x2 %>% rename(transcript_id = X, logFC_wt_4x = logFC)
c_4x2 <- c_4x2 %>% rename(transcript_id = X, logFC_c_4x = logFC)

# Merge dataframes by transcript_id
merged_data <- wt_c2 %>%
  dplyr::select(transcript_id, logFC_wt_c) %>%
  inner_join(wt_4x2 %>% dplyr::select(transcript_id, logFC_wt_4x), by = "transcript_id") %>%
  inner_join(c_4x2 %>% dplyr::select(transcript_id, logFC_c_4x), by = "transcript_id")

# Invert c_4x logFC so positive means higher in Control
merged_data$logFC_c_inv <- -merged_data$logFC_c_4x

# Calculate correlation coefficients (using inverted c_4x)
cor_wt_c <- -cor(merged_data$logFC_wt_c, merged_data$logFC_c_inv, use = "complete.obs")
cor_4x_c <- -cor(merged_data$logFC_wt_4x, merged_data$logFC_c_inv, use = "complete.obs")
cor_wt_4x <- -cor(merged_data$logFC_wt_c, merged_data$logFC_wt_4x, use = "complete.obs")

# Get axis limits for quadrant labels
x_range_a <- range(merged_data$logFC_c_inv, na.rm = TRUE)
y_range_a <- range(merged_data$logFC_wt_c, na.rm = TRUE)
x_range_b <- range(merged_data$logFC_c_inv, na.rm = TRUE)
y_range_b <- range(merged_data$logFC_wt_4x, na.rm = TRUE)
x_range_c <- range(merged_data$logFC_wt_4x, na.rm = TRUE)
y_range_c <- range(merged_data$logFC_wt_c, na.rm = TRUE)

# Create plot A: WT vs Control (both positive = higher in WT and Control)
plot_a <- ggplot(merged_data, aes(x = logFC_c_inv, y = logFC_wt_c)) + # c_4x inverted
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  annotate("text", x = -Inf, y = Inf, 
           label = paste0("r = ", round(cor_wt_c, 2)), 
           hjust = -0.1, vjust = 1.5, size = 5, fontface = "bold") +
  annotate("text", x = x_range_a[2] * 0.7, y = y_range_a[2] * 0.85,
           label = "Higher in WT\nHigher in Control", size = 3, color = "gray30") +
  annotate("text", x = x_range_a[1] * 0.7, y = y_range_a[2] * 0.85,
           label = "Higher in WT\nHigher in 4x", size = 3, color = "gray30") +
  annotate("text", x = x_range_a[1] * 0.7, y = y_range_a[1] * 0.85,
           label = "Higher in Control\nHigher in 4x", size = 3, color = "gray30") +
  annotate("text", x = x_range_a[2] * 0.7, y = y_range_a[1] * 0.85,
           label = "Higher in Control\nHigher in Control", size = 3, color = "gray30") +
  labs(x = "logFC2(Control/4x)", # positive = higher in Control
       y = "logFC2(WT/Control)", # positive = higher in WT
       title = "A") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

# Create plot B: 4x vs Control (both positive = higher in 4x and Control)
plot_b <- ggplot(merged_data, aes(x = logFC_c_inv, y = logFC_wt_4x)) + # c_4x inverted
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  annotate("text", x = -Inf, y = Inf, 
           label = paste0("r = ", round(cor_4x_c, 2)), 
           hjust = -0.1, vjust = 1.5, size = 5, fontface = "bold") +
  annotate("text", x = x_range_b[2] * 0.7, y = y_range_b[2] * 0.85,
           label = "Higher in 4x\nHigher in Control", size = 3, color = "gray30") +
  annotate("text", x = x_range_b[1] * 0.7, y = y_range_b[2] * 0.85,
           label = "Higher in 4x\nHigher in 4x", size = 3, color = "gray30") +
  annotate("text", x = x_range_b[1] * 0.7, y = y_range_b[1] * 0.85,
           label = "Higher in WT\nHigher in 4x", size = 3, color = "gray30") +
  annotate("text", x = x_range_b[2] * 0.7, y = y_range_b[1] * 0.85,
           label = "Higher in WT\nHigher in Control", size = 3, color = "gray30") +
  labs(x = "logFC2(Control/4x)", 
       y = "logFC2(4x/WT)",
       title = "B") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

# Create plot C: WT vs 4x
plot_c <- ggplot(merged_data, aes(x = logFC_wt_4x, y = logFC_wt_c)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  annotate("text", x = -Inf, y = Inf, 
           label = paste0("r = ", round(cor_wt_4x, 2)), 
           hjust = -0.1, vjust = 1.5, size = 5, fontface = "bold") +
  annotate("text", x = x_range_c[2] * 0.7, y = y_range_c[2] * 0.85,
           label = "Higher in WT\nHigher in 4x", size = 3, color = "gray30") +
  annotate("text", x = x_range_c[1] * 0.7, y = y_range_c[2] * 0.85,
           label = "Higher in WT\nHigher in WT", size = 3, color = "gray30") +
  annotate("text", x = x_range_c[1] * 0.7, y = y_range_c[1] * 0.85,
           label = "Higher in Control\nHigher in WT", size = 3, color = "gray30") +
  annotate("text", x = x_range_c[2] * 0.7, y = y_range_c[1] * 0.85,
           label = "Higher in Control\nHigher in 4x", size = 3, color = "gray30") +
  labs(x = "logFC2(4x/WT)", 
       y = "logFC2(WT/Control)",
       title = "C") +  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

combined_plot <- grid.arrange(plot_a, plot_b, plot_c, nrow = 1)
ggsave("logFC_correlation_plots.pdf", combined_plot, width = 15, height = 5)
ggsave("logFC_correlation_plots.png", combined_plot, width = 15, height = 5, dpi = 300)

# get correlation 
library(ppcor)  # For partial correlations
# Calculate partial correlations # c_4x
pcor_a <- pcor.test(merged_data$logFC_wt_c, merged_data$logFC_c_4x, 
                    merged_data$logFC_wt_4x, method = "pearson")
                    # r for WT/C vs 4x/C controlling for 4x/WT = -0.73
pcor_b <- pcor.test(merged_data$logFC_wt_4x, merged_data$logFC_c_inv,
                    merged_data$logFC_wt_c, method = "pearson")
                    # r for WT/C vs 4x/C controlling for 4x/WT = -0.73
pcor_c <- pcor.test(merged_data$logFC_wt_c, merged_data$logFC_wt_4x,
                    merged_data$logFC_c_inv, method = "pearson")
# Create a summary table of correlations
correlation_summary <- data.frame(
  Comparison = c("A: WT vs Control", "B: 4x vs Control", "C: WT vs 4x"),
  Pearson_r = round(c(cor_wt_c, cor_4x_c, cor_wt_4x), 3),
  Partial_r = round(c(pcor_a$estimate, pcor_b$estimate, pcor_c$estimate), 3),
  Controlled_for = c("4x", "WT", "Control"),
  Partial_p_value = format.pval(c(pcor_a$p.value, pcor_b$p.value, pcor_c$p.value), digits = 3) )

print(correlation_summary) 

####### 5 Sleuth across genes ###########

so2 <- sleuth_prep(s2c, target_mapping=t2g, aggregation_column='ens_gene',
                  extra_bootstrap_summary=T, read_bootstrap_tpm=T)

# Next, we will smooth the tpm per sample using a parameter
# based on our model - so here we estimate parameters for
# response error measurement (full) model
# this is our alternative model with DE
so2 <- sleuth_fit(so2, ~condition, 'full')  # XXX passed filter
str(so2)
# get our null model r where the isoform levels are equal
so2 <- sleuth_fit(so2, ~1  , 'r')
# we compare our null and alternative models
# differential analysis using a likelihood ratio test (LRT)
so2 <- sleuth_lrt(so2, 'r', 'full')
sr <- sleuth_results(so2, 'r:full', 'lrt', show_all=F)
out11 <- dplyr::filter(sr, qval <= 0.05)  # 35793 transcripts from x genes
str(unique(out11$target_id)) # 7101 unique genes
str( (out11$target_id)) # 28,648 transcripts
write.csv(out11, "sleuth.DE.genes.csv", quote=T, row.names=T)
write.csv(dplyr::filter(sr, qval <=1), "sleuth.all.genes.csv", quote=T, row.names=T)
models(so2) # check model details
##########