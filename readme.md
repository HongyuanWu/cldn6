
```
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
tab <- getGeneKEGGLinks(species="hsa")
tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID,column="SYMBOL", keytype="ENTREZID")
head(tab)
kegg<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/cldn6/master/cldn6.kegg.txt",as.is=T)[,1]
out<-unique(tab[tab$PathwayID %in% paste("path:",kegg,sep=""),])
write.table(out,file="cldn6.kegg.list.txt",quote=F,col.names=T,row.names=F,sep="\t")
write.csv(out,file="cldn6.kegg.list.csv",quote=F)
```
