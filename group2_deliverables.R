library(affy)
library(affyPLM)
library(hgu133plus2.db)
library(AnnotationDbi)
library(simpleaffy)
library(gcrma)

gse <- ReadAffy(celfile.path="GSE19804_RAW")

dataset_removed = gse[-c(1,11,34:44)]
normalized = rma(dataset_removed)
exprs_data= as.data.frame(exprs(normalized))


x <- hgu133plus2.db 
library(dplyr)
exprs_data$PROBEID <- as.character(rownames(exprs_data))
exprs_annotation <- AnnotationDbi::select(x,keys=exprs_data$PROBEID, columns = "SYMBOL")

exprs_data$SYMBOL = exprs_annotation[!duplicated(exprs_annotation$PROBEID),]$SYMBOL

exprs_data = exprs_data[!duplicated(exprs_data$SYMBOL),]
exprs_data = na.omit(exprs_data)
row.names(exprs_data) = exprs_data$SYMBOL
exprs_data = exprs_data[-c(32,33)]


library(EnhancedVolcano)
library(limma)

interest <- factor(paste(rep(c("C", "N"), each = 60), sep = ""))
design.mat <- model.matrix(~ 0 + interest)


fit <- limma::lmFit(exprs_data, design.mat)

make.names(c("interestlung cancer", "interestpaired normal adjacent"), unique = TRUE)
contrast.matrix <- makeContrasts(
  interestC-interestN, 
  levels = design.mat)

fit.contrast = contrasts.fit(fit, contrast.matrix)

efit <- eBayes(fit.contrast)

genes=geneNames(gse)
limma_output <- topTable(efit, n = 54670, p.value=0.05, adjust.method="fdr")

EnhancedVolcano( toptable = limma_output, 
                 lab =rownames(exprs_data), 
                 x = "logFC", 
                 y = "P.Value")
                 
#still fixing

