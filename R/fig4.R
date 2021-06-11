
##########  Figure 4a -----------

GO_results <- read.csv2("gProfiler_sscrofa_A_vs_C_2metods.txt", sep = '\t')
head(GO_results)

GO_results <- GO_results[,c(1,3,2,10,4)]
head(GO_results)
colnames(GO_results) <- c("Category","ID","Term","Genes","adj_pval")

res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$gene_name <- ifelse(is.na(res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$gene_name), res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$ref_gene_id, res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$gene_name)


res_GO <- res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein[,c("gene_name","log2FoldChange","baseMean","stat","pvalue","padj")]
colnames(res_GO) <- c("ID","logFC","AveExpr","t","P.Value","adj.P.Val")

circ <- circle_dat(GO_results, res_GO)
circ$adj_pval <- as.numeric(as.character(circ$adj_pval))

circ_GO <- circ[circ$category %in% c("GO:MF","GO:BP","GO:CC"),]

IDs <- c("GO:GO:0006954","GO:0001816","GO:0006915","GO:0012501","GO:0080134","GO:0071345","GO:0001525","GO:0001944","GO:0001568","GO:0007249")
IDs2 <- c("GO:0070851","GO:0001525","GO:1902494")


tiff("GO2.tiff", width = 250, height = 250, units = 'mm', res = 600, compression = "lzw")
GOBubble(circ_GO, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 40)  
dev.off()

KEGG_circos <- c("KEGG:04668","KEGG:04064","KEGG:04657","KEGG:04210")
circ_KEGG <- circ[circ$ID %in% KEGG_circos,]

KEGG_genes <- res_GO[res_GO$ID %in% circ_KEGG$genes,]
KEGG_process <- unique(circ_KEGG$term)

chord <- chord_dat(circ_KEGG, KEGG_genes, KEGG_process)


##########  Figure 4b -----------
tiff("GO.tiff", width = 250, height = 250, units = 'mm', res = 600, compression = "lzw")
GOCircle(circ_GO, nsub = IDs)
dev.off()
