######circos creation

FPKM_sig_DEGs <- gene_FPKM[rownames(gene_FPKM) %in% res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$Row.names,]
FPKM_sig_lncRNA <- FPKM_sig_all[rownames(FPKM_sig_all) %in% res_sig_table_wspolne_biomart_uniq_lncRNA$Row.names,]

FPKM_sig_DEGs <- as.data.frame(FPKM_sig_DEGs)
FPKM_sig_DEGs$nr <- seq(1:nrow(FPKM_sig_DEGs))
head(FPKM_sig_DEGs)
FPKM_sig_DEGs_null <- as.matrix(FPKM_sig_DEGs[,1:11])
rownames(FPKM_sig_DEGs_null) <- NULL 
rownames(FPKM_sig_DEGs_null) <- seq(1:nrow(FPKM_sig_DEGs_null))
head(FPKM_sig_DEGs_null)
colnames(FPKM_sig_DEGs_null) <- NULL


DE_all_FPKM <- FPKM_sig_DEGs_null
DE_all_FPKM <- as.matrix(DE_all_FPKM)
head(DE_all_FPKM)

DE_all_FPKM_null = log2(DE_all_FPKM+1)
DE_all_FPKM_null = t(scale(t(DE_all_FPKM_null), scale=F))

DE_all_FPKM_null<- DE_all_FPKM_null[,c(6,7,8,9,10,11,1,2,3,4,5)] 


split = rep("DEGs", 1474)
split = factor(split, levels = "DEGs")
summary(DE_all_FPKM_null)
col_fun1 = colorRamp2(c(-3, 0, 3), c("green", "black", "red"))


names_circos_DEGs <- merge(FPKM_sig_DEGs_to_links,res_sig_table_wspolne_biomart_uniq[,c("Row.names","gene_name")], by.x=0, by.y="Row.names" )
names_circos_lncRNA <- FPKM_sig_lncRNA_to_links
names_circos_lncRNA$name <- rownames(names_circos_lncRNA)

rownames(DE_all_FPKM_null) <- c(as.character(names_circos_DEGs$gene_name), as.character(rownames(names_circos_lncRNA)) )

##########  Figure 3 (cicros heatmap) -----------

dev.off()
circos.clear()
tiff("circos.tiff", width = 250, height = 250, units = 'mm', res = 600, compression = "lzw")
{
  circos.par(gap.after = 8,track.margin=c(0.02, 0.01))
  circos.heatmap(DE_all_FPKM_null, split = split, col = col_fun1, dend.side = NULL, track.height = 0.3)
  
  res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein = res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein[match(rownames(FPKM_sig_DEGs),res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$Row.names),]
  row_mean = res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$log2FoldChange
  
  names(row_mean) <- 1:1474
  circos.track(ylim = range(row_mean), track.height = 0.15, panel.fun = function(x, y) {
    y = row_mean[CELL_META$subset]
    y = y[CELL_META$row_order]
    circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
    circos.points(seq_along(y) - 0.5, y, col = ifelse(y > 0, "red", "blue"))
  }, cell.padding = c(0.02, 0, 0.02, 0))
  
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == 1) { 
      circos.text(CELL_META$cell.xlim[2] + convert_x(5.5, "mm"), 3,
                  "controls", cex = 0.6, facing = "inside")
      circos.text(CELL_META$cell.xlim[2] + convert_x(6.1, "mm"), 8,
                  "thrombosis", cex = 0.6, facing = "inside")
    }
  }, bg.border = NA)
  FPKM_sig_DEGs_to_links <- FPKM_sig_DEGs
  
  sim_matrix5_unique_links <- merge(sim_matrix5_unique, FPKM_sig_DEGs_to_links, by.x="Var1" , by.y=0)
  sim_matrix5_unique_links <- merge(sim_matrix5_unique_links, FPKM_sig_DEGs_to_links, by.x="Var2" , by.y=0)
  head(sim_matrix5_unique_links)
  sim_matrix5_unique_links <- sim_matrix5_unique_links[,c("Var1", "Var2", "value", "nr.x", "nr.y")]
  
  df_link = data.frame(
    from_index = sim_matrix5_unique_links$nr.x,
    to_index = sim_matrix5_unique_links$nr.y,
    value = sim_matrix5_unique_links$value
  )
  
  head(df_link)
  for(i in seq_len(nrow(df_link))) {
    circos.heatmap.link(df_link$from_index[i],
                        df_link$to_index[i],
                        col = ifelse(df_link$value[i] > 0.7, "#00AFBB", "#E7B800"))
  }
  
  
}

dev.off()
