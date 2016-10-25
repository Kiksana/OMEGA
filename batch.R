# Set a logfile and working directory
current <- getwd()
workdir <- setwd("OMEGA_R")
logfile <- file("log_R.txt", open="wt")
# Load dependencies
library("DESeq2")
# # Open matrix files
read_mat <- read.delim("counts_mat.txt", row.names=1)
cond_mat <- read.delim("condition.txt", row.names=1)
divers_mat_A <- read.delim("alpha_diversity_A.txt", row.names=1)
divers_mat_B <- read.delim("alpha_diversity_B.txt", row.names=1)
# Creating density plot from read count files. Saving as “.jpg”
d_a <- density(t(divers_mat_A))
d_b <- density(t(divers_mat_B))
jpeg(width = 835, height = 577, "dens_plot_cond.jpg")
plot(d_a, col = "blue", main = "Density plot of condition A and B", xlab = "number of clades")
lines(d_b, col = "red")
legend("topright",levels(cond_mat$condition), fill = c("blue", "red"))
dev.off()
sink(logfile, type="message")
# Generating DESeq2 objects, applying custom geometries for low read counts
dds <- DESeqDataSetFromMatrix(countData = read_mat, colData = cond_mat, design = ~ condition)
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds$condition <- factor(dds$condition, levels=c("A","B"))
cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
# run DESeq2 and save results
dds <- DESeq(dds, quiet = TRUE, fitType = "mean", test = "Wald")
resMLE <- results(dds, alpha = 0.05, addMLE = TRUE)
res <- results(dds, alpha = 0.05)
# Save images and result files
jpeg(width = 835, height = 577, "MAplot shrunken vs unshrunken.jpg")
par(mfrow=c(1,2), pch=20, oma = c(0,0,2,0))
plotMA(resMLE, main="MAplot of unshrunken DESeq2 estimates", ylim=c(-7,7), MLE = TRUE, alpha = 0.05, ylab = "MLE log2 fold change", xlab = "mean normalized counts for each clade", cex=0.7, type = "p")
plotMA(res, main="MAplot of shrunken DESeq2 estimates", ylim=c(-7,7), MLE = FALSE, alpha = 0.05, ylab = "MAP log2 fold change", xlab = "mean normalized counts for each clade", cex=0.7, type = "p")
dev.off()
jpeg(width = 835, height = 577, "min_padj_normcounts.jpg")
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()
jpeg(width = 835, height = 577, "min_padj_log2counts.jpg")
plotCounts(dds, gene=which.min(res$padj), intgroup="condition", transform = TRUE)
dev.off()
res2Ordered <- res[order(res$padj),]
write.csv(as.data.frame(res2Ordered), file="results.csv")
sink()