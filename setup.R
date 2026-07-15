
# ---- SCA2/PD Project Setup ----
library(Biobase)
library(GEOquery)
library(limma)
library(annotate)

# load saved data
load("data/processed/eset_final.RData")
load("data/processed/fit2.RData")

# rebuild groups and results
groups <- pData(eset_final)$group
results_SCA2 <- topTable(fit2, coef = "SCA2_vs_Control", n = Inf)
results_PD <- topTable(fit2, coef = "PD_vs_Control", n = Inf)

cat("Setup complete! eset_final, fit2, results_SCA2, and results_PD are ready.
")


library(ggplot2)
library(pheatmap)

results_PD$significant <- ifelse(results_PD$adj.P.Val < 0.05, TRUE, FALSE)
results_SCA2$significant <- ifelse(results_SCA2$adj.P.Val < 0.05, TRUE, FALSE)

sig_genes <- unique(c(
  results_PD$GENE_SYMBOL[results_PD$significant == TRUE],
  results_SCA2$GENE_SYMBOL[results_SCA2$significant == TRUE]
))

expr_mat <- exprs(eset_final)[fData(eset_final)$GENE_SYMBOL %in% sig_genes, ]
rownames(expr_mat) <- fData(eset_final)$GENE_SYMBOL[fData(eset_final)$GENE_SYMBOL %in% sig_genes]

annotation <- data.frame(Group = groups)
rownames(annotation) <- colnames(expr_mat)

library(Biobase)
library(ggplot2)
load("data/processed/progress.RData")
