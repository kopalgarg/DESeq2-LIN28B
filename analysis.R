setwd('~/Desktop/DESeq2_LIN28B')
require(DESeq2)
require(tximport)
require(edgeR)
require(limma)
require(dplyr)
require(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(BuenColors)
library(gridExtra)

dir = '/Volumes/broad_sankaranlab/kgarg/cMyb/out'

samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples$condition <- factor(rep(c("PB", "CB"), each = 24))
samples$time = factor(rep(c(1, 2, 3, 4, 5, 6, 7, 8), each = 3))
files <- file.path(dir, paste0(samples$sample, '_',samples$run, ".genes.results"))

names(files) <- paste0(samples$sample, '_',samples$run)

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

# DESeqDataSet object 
dds = DESeqDataSetFromMatrix(countData = round(txi.rsem$counts),
                       colData = samples,
                       design= ~ condition + time)

# remove rows w/ no counts or only 1 count across all samples

keep <- rowSums(counts(dds)) > 1 

dds <- dds[keep,] # from 60662 to 33303

# VST
vsd <- vst(dds, blind = FALSE)

# rlog-transformation
rld <- rlog(dds, blind = FALSE)

# estimate size factors for log2(count + 1)
dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 2:3]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)


# 4.2 scatterplot of transformed counts from two samples (L:log2 of normalized counts, 
#                                                     M:using the VST, 
#                                                     R:using the rlog)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

# 4.3 sample distances: calculate the Euclidean distance between sample

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# heatmap of Euclidean distances using the variance stabilizing transformed values.

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# using Poisson Distance
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$dex, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# PCA: using VST
plotPCA(vsd, intgroup = c("condition", "sample"))

pcaData <- plotPCA.san(vsd, intgroup = c( "condition", "sample"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# specify cell line (symbol) and condition (color)
p3 = ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC3: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_manual(values = jdb_palette("corona")) + pretty_plot()
# MDS
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

# 9. Time Course

ddsTC <- DESeq(dds)

# resTC = readRDS('resTC.rds')
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol

tc <- plotCounts(ddsTC, which.min(resTC$padj), returnData = TRUE)
tc$stage = factor(rep(c('Pro', 'Poly', 'Ortho', 'lBaso', 'eBaso', 'CFUE', 'CD34', 'BFUE'), each = 3))
tc$stage <- factor(tc$stage, levels=c('CD34', 'BFUE', 'CFUE', 'Pro', 
                                      'eBaso', 'lBaso', 'Poly', 'Ortho'))
p1 = ggplot(pcaData, aes(x = PC1, y = PC2, color = stage, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC3: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_manual(values = jdb_palette("corona")) + pretty_plot()

p2 = ggplot(tc,
       aes(x = stage, y = count, color = stage, group = condition, shape = condition)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10() + xlab('stage') + ylab('count') + pretty_plot() +
  scale_color_manual(values = jdb_palette("corona")) 

ggsave(grid.arrange(pc12, pc13,p2), file = "out.pdf", width = 7, height = 10)

pcaData <- plotPCA.san(vsd, intgroup = c( "condition", "sample"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$stage = factor(rep(c('Pro', 'Poly', 'Ortho', 'lBaso', 'eBaso', 'CFUE', 'CD34', 'BFUE'), each = 3))
ggplot(pcaData, aes(x = PC1, y = PC2, color = stage, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_manual(values = jdb_palette("corona")) + pretty_plot()

plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                  intgroup.df, name = colData(rld)[,1])
  if (returnData) {
    attr(d, "percentVar") <- cbind(percentVar[1], percentVar[2])
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}


# Differential Expression analysis
countData = as.data.frame(round(txi.rsem$counts))
countData$genes = rownames(countData)
RNA.counts <- countData[,1:48]
meta <- stringr::str_split_fixed(colnames(RNA.counts), "_", 4)
RNA.counts.df <- as.data.frame(RNA.counts)
RNA.condition = c('PB_proerythroblast_1', 'PB_proerythroblast_2','PB_proerythroblast_3',
                  'PB_polychromatic_1', 'PB_polychromatic_2', 'PB_polychromatic_3', 
                  'PB_orthochromatic_1', 'PB_orthochromatic_2', 'PB_orthochromatic_3', 
                  'PB_late_basophilic_1', 'PB_late_basophilic_2', 'PB_late_basophilic_3', 
                  'PB_early_basophilic_1', 'PB_early_basophilic_2', 'PB_early_basophilic_3', 
                  'PB_CFU_1', 'PB_CFU_2', 'PB_CFU_3', 
                  'PB_CD34_1', 'PB_CD34_2', 'PB_CD34_3', 
                  'PB_BFU_1', 'PB_BFU_2', 'PB_BFU_3', 
                  'CB_proerythroblast_1', 'CB_proerythroblast_2', 'CB_proerythroblast_3', 
                  'CB_polychromatic_1', 'CB_polychromatic_2', 'CB_polychromatic_3', 
                  'CB_orthochromatic_1', 'CB_orthochromatic_2', 'CB_orthochromatic_3', 
                  'CB_late_basophilic_1',  'CB_late_basophilic_2',  'CB_late_basophilic_3', 
                  'CB_early_basophilic_1',   'CB_early_basophilic_2',   'CB_early_basophilic_3', 
                  'CB_CFU_1',  'CB_CFU_2',  'CB_CFU_3', 
                  'CB_CD34_1','CB_CD34_2', 'CB_CD34_3',  
                  'CB_BFU_1', 'CB_BFU_2', 'CB_BFU_3' )
colData <- as.data.frame(RNA.condition)
colnames(RNA.counts.df) = RNA.condition
row.names(colData) <- colnames(RNA.counts.df)

conditionPB=as.data.frame(factor(rep(c('Pro_PB', 'Poly_PB', 'Ortho_PB', 'lBaso_PB', 'eBaso_PB', 'CFUE_PB', 'CD34_PB', 'BFUE_PB'), each = 3)))
colnames(conditionPB)='condition'
conditionPB %>% mutate_if(is.factor, as.character) -> conditionPB
conditionCB=as.data.frame(factor(rep(c('Pro_CB', 'Poly_CB', 'Ortho_CB', 'lBaso_CB', 'eBaso_CB', 'CFUE_CB', 'CD34_CB', 'BFUE_CB'), each = 3)))
colnames(conditionCB)='condition'
conditionCB %>% mutate_if(is.factor, as.character) -> conditionCB
conditionPB[25:48,]=conditionCB[1:24,]
condition=conditionPB
condition=as.data.frame(condition)
RNA.dds <- DESeqDataSetFromMatrix(countData = RNA.counts.df, colData = condition, design = ~condition)
RNA.dds <- DESeq(RNA.dds)

comp <- expand.grid(unique(condition$condition),unique(condition$condition), stringsAsFactors = FALSE)
comp <- subset(comp, Var1!=Var2)
comp$cond <- "condition"
comp <- comp %>% 
  mutate(key = paste0(pmin(Var1, Var2), "_" ,pmax(Var1, Var2), sep = "")) %>%
  dplyr::distinct(key, .keep_all=TRUE)

sig.pairwise <- function(x) {
  print(x)
  # Do the DESeq contrast
  resSig.RNA <- data.frame(results(RNA.dds, contrast=c(comp[x,3], comp[x,2], comp[x,1]), filterFun=ihw, alpha = 0.01))
  resSig.RNA$gene <- countData$genes
  resSig.RNA %>% filter(complete.cases(.)) %>% arrange(padj) %>% filter(padj < 0.01) -> out
  # Do some rounding
  out$baseMean <- round( out$baseMean, 1)
  out$log2FoldChange <- round( out$log2FoldChange, 1)
  out$pvalue <-sprintf("%.3e", out$pvalue)
  out$padj <-sprintf("%.3e", out$padj)
  # Output table
  out_table = out[,c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")]
  write.table(out[,c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")],
              file = paste0("", comp[x,4], "_Padj01.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  comp[x,4]
}
all.pairwise.P1 <- function(x) {
  print(x)
  # Do the DESeq contrast
  resSig.RNA <- data.frame(results(RNA.dds, contrast=c(comp[x,3], 'lBaso', 'Pro'), filterFun=ihw, alpha = 0.01))
  resSig.RNA$gene <- countData$genes
  resSig.RNA %>% filter(complete.cases(.)) %>% arrange(padj) -> out
  # Do some rounding
  out$baseMean <- round( out$baseMean, 1)
  out$log2FoldChange <- round( out$log2FoldChange, 1)
  out$pvalue <-sprintf("%.3e", out$pvalue)
  out$padj <-sprintf("%.3e", out$padj)
  # Output table
  out_table = out[,c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")]
  comp[x,4]
  write.table(out[,c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")],
              file = paste0("", comp[x,4], "_Padj01.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

dumb <- lapply(1:dim(comp)[1], sig.pairwise)

#strongest upregulation
resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ]


require(annotables)

grch38 = grch38 %>% 
  dplyr::select(ensgene, symbol, chr, start, end, description) %>%
  subset(., select = c(ensgene, symbol, chr))

countData = countData %>% mutate(symbol=factor(unlist(lapply(strsplit(genes, split = "_"), function(x) x[2])))) %>%
  mutate(ensgene=factor(unlist(lapply(strsplit(genes, split = "_"), function(x) x[1])))) %>%
  mutate(ensgene=factor(unlist(lapply(strsplit(ensgene, split = "\\."), function(x) x[1]))))

countData = inner_join(countData, grch38)

# make long  so that you have one column that has replicate (poly_cb ...) then has gene, chr, count etc
countData_long <- gather(countData, replicate, count, hs_PB_proerythroblast_1:CB_BFU_3, factor_key=TRUE)

# group_by on the replicate and chr

countData_long = countData_long %>% group_by(replicate, chr) %>% summarise(count = sum(count))
countData_long = countData_long %>% subset(chr %in% c('X', 'Y', 'MT'))

countData_wide = spread(countData_long, chr, count)

# volcano plot

volcano_plot <- function(data, cell_type){
  data = data %>% mutate(symbol=factor(unlist(lapply(strsplit(gene, split = "_"), function(x) x[2]))))
threshold_OE <- data$padj < 0.05
data$threshold <- threshold_OE 
data <- data[order(data$padj), ] 
data$genelabels <- ""
data$genelabels[1:10] <-T

ggplot(data) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text_repel(aes(x = log2FoldChange, 
                      y = -log10(padj), 
                      label = ifelse(genelabels == T, as.character(symbol),""))) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none",
        legend.title = element_blank()) +
  labs(colour = NULL) + 
  ggtitle(cell_type) +
  pretty_plot()
}


