##### Figure 2b ----------

do_venn_AvC <- attr(venn(list(ballgown=as.character(res_sig_FC1$id), DESeq=as.character(rownames(res_sig_DESeq_FC1)))), "intersection")$'A:B'

res_sig_DESeq_FC1_compilation <- res_sig_DESeq_FC1[rownames(res_sig_DESeq_FC1) %in% do_venn_AvC,] 
res_sig_DESeq_FC1_compilation_trans <- merge(res_sig_DESeq_FC1_compilation, gtf_stringtie_zator_transcript, by.x=0, by.y="gene_id")
res_sig_DESeq_FC1_compilation_trans <- res_sig_DESeq_FC1_compilation_trans[order(res_sig_DESeq_FC1_compilation_trans$Row.names, res_sig_DESeq_FC1_compilation_trans$ref_gene_id, decreasing = F),]
res_sig_DESeq_FC1_compilation_trans <- res_sig_DESeq_FC1_compilation_trans[!duplicated(res_sig_DESeq_FC1_compilation_trans$Row.names),]
res_sig_DESeq_FC1_compilation_trans_biomart <- merge(res_sig_DESeq_FC1_compilation_trans, allgenes.Ensembl, by.x="ref_gene_id" , by.y="ensembl_gene_id", all.x=T)

listaGenow_A_C <- res_sig_DESeq_FC1_compilation_trans_biomart$log2FoldChange
names(listaGenow_A_C) <- res_sig_DESeq_FC1_compilation_trans_biomart$entrezgene_id
path <- pathview(gene.data = listaGenow_A_C, pathway.id = "ssc04668", species = "ssc")
path <- pathview(gene.data = listaGenow_A_C, pathway.id = "ssc04064", species = "ssc")
path <- pathview(gene.data = listaGenow_A_C, pathway.id = "ssc04657", species = "ssc")
path <- pathview(gene.data = listaGenow_A_C, pathway.id = "ssc04210", species = "ssc")
path <- pathview(gene.data = listaGenow_A_C, pathway.id = "ssc04115", species = "ssc")
path <- pathview(gene.data = listaGenow_A_C, pathway.id = "ssc04630", species = "ssc")

write.csv2(res_sig_DESeq_FC1_compilation_trans_biomart, "AvsC_DEGs_2methods.csv")
write.table(unique(res_sig_DESeq_FC1_compilation_trans$gene_name), "gene_name_AvsC_2methods.txt", col.names = F, quote = F, row.names = F)


res$color <- "black"
res[is.na(res$padj),]$padj <- 0.99

res[(rownames(res) %in% do_venn_AvC) & (res$log2FoldChange >= 0.5),]$color <- "red"
head(res)

res[(rownames(res) %in% do_venn_AvC) & res$log2FoldChange <= -0.5,]$color <- "green"
table(res$color)


##### Figure 2d ----------

tiff("Volcano.tiff", width = 140, height = 140, units = 'mm', res = 500, compression = "lzw")
myvolcano3(res$log2FoldChange, -log10(res$padj),-0.5, 0.5,-log10(0.05), main="", kolory = res$color, maksY=16, limitX = 10.8)
points(res[res$color %in% "red",]$log2FoldChange, -log10(res[res$color %in% "red",]$padj), col="red", cex = 0.6, pch=19)
points(res[res$color %in% "green",]$log2FoldChange, -log10(res[res$color %in% "green",]$padj), col="green3", cex = 0.6, pch=19)
dev.off()

##### Figure 2c ----------

tiff("MA.tiff", width = 140, height = 140, units = 'mm', res = 500, compression = "lzw")
myplotMA4(myx = log10(res$baseMean+1), myy = res$log2FoldChange, limitY = 10.8, limitX = c(0,7), abline = 0.5 , kolory = res$color, qval = res$padj)
points(log10(res[res$color %in% "red",]$baseMean +1), res[res$color %in% "red",]$log2FoldChange, col="red", cex = 0.6, pch=19)
points(log10(res[res$color %in% "green",]$baseMean +1), res[res$color %in% "green",]$log2FoldChange, col="green3", cex = 0.6, pch=19)
dev.off()
