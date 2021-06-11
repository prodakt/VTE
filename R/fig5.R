
##########  Figure 5 -----------
tiff("circos_KEGG.tiff", width = 250, height = 250, units = 'mm', res = 600, compression = "lzw")
  GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5,ribbon.col = c("red","green","blue","grey"))
dev.off()
