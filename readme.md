
```
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
tab <- getGeneKEGGLinks(species="hsa")
tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID,column="SYMBOL", keytype="ENTREZID")
head(tab)
kegg<-read.table("")
```
