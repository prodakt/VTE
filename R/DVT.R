
source("funs.R")
source("libs.R")


ensemblSs = useMart(host="ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="sscrofa_gene_ensembl") 

allgenes.Ensembl = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype", "entrezgene_id", "description", "entrezgene_accession"),mart=ensemblSs)
gtf_stringtie_zator <- import.gff("stringtie_merged.gtf")
gtf_stringtie_zator <- as.data.frame(gtf_stringtie_zator)


pheno_data <- read.csv2("pheno_data.csv" , header = T, sep = ";") 
pheno_data_AvsC <- pheno_data[pheno_data$group %in% c("A","C"),]
pheno_data_AvsC <- pheno_data_AvsC[order(dir("ballgown_AvsC/")),]
colnames(pheno_data_AvsC)[2] <- "treatment"
pheno_data_AvsC$treatment <- as.character(pheno_data_AvsC$treatment)

# ballgown -----

rownames(pheno_data_AvsC) <- dir("ballgown_AvsC/")
bg <- ballgown(dataDir = "ballgown_AvsC/", samplePattern=".", pData=pheno_data_AvsC)
bg_filt <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset=TRUE)
gene_FPKM = as.data.frame(gexpr(bg_filt))
colnames(gene_FPKM) <- pheno_data_AvsC$ids

results_genes <-  stattest(bg_filt, feature='gene', covariate='treatment', 
                           getFC=TRUE, meas='FPKM')

res_sig <- subset(results_genes, results_genes$qval <= 0.05)
res_sig$log2FC <- log2(res_sig$fc)
res_sig_FC1 <- subset(res_sig, abs(res_sig$log2FC) > 1)
nrow(res_sig_FC1)


gtf_stringtie_zator_transcript <- gtf_stringtie_zator[gtf_stringtie_zator$type %in% "transcript",]
gtf_stringtie_zator_transcript_DEGs_AvsC <- gtf_stringtie_zator_transcript[gtf_stringtie_zator_transcript$gene_id %in% res_sig_FC1$id,]


gtf_stringtie_zator_novel_trans <- gtf_stringtie_zator[gtf_stringtie_zator$gene_id %in% ENSEMBL_annot,]
gtf_stringtie_zator_novel_trans <- gtf_stringtie_zator_novel_trans[gtf_stringtie_zator_novel_trans$type %in% "transcript",]
gtf_stringtie_zator_novel_trans <- gtf_stringtie_zator_novel_trans[is.na(gtf_stringtie_zator_novel_trans$ref_gene_id),]


# DESeq2 -------

count_table <- read.csv2("gene_count_matrix.csv", sep=",", row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = count_table[,colnames(count_table) %in% pheno_data_AvsC$ids], colData = pheno_data_AvsC, design = ~ treatment)
dds$treatment <- relevel(dds$treatment, ref = "ctr")
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- DESeq(dds)
res <- results(dds)
res_sig_DESeq <- subset(res, padj < 0.05)
res_sig_DESeq <- as.data.frame(res_sig_DESeq)
res_sig_DESeq_FC1 <- subset(res_sig_DESeq, abs(res_sig_DESeq$log2FoldChange) > 1)




#validation 

validation_genes_up <- c("IL1R1", "TNFRSF1A", "TNFRSF1B", "CD40", "TRAF2", "TRAF3", "NFKB1", "NFKB2","BIRC3","PLAU","PLAT","VEGFC","CASP3","CASP7","CASP10", "IL6","CCL2","CCL19", "CXCL2")
validation_genes_down <- c("TAB1","MAPK14","VEGFD", "VEGFB", "TNFSF10","SPTAN1", "PARP1", "PARP2", "DFFA" )

res_sig_DESeq_FC1_compilation_trans_biomart_up_validation <- res_sig_DESeq_FC1_compilation_trans_biomart[res_sig_DESeq_FC1_compilation_trans_biomart$entrezgene_accession %in% validation_genes_up,]
res_sig_DESeq_FC1_compilation_trans_biomart_down_validation <- res_sig_DESeq_FC1_compilation_trans_biomart[res_sig_DESeq_FC1_compilation_trans_biomart$entrezgene_accession %in% validation_genes_down,]

gene_FPKM_valid_up <- gene_FPKM[rownames(gene_FPKM) %in% res_sig_DESeq_FC1_compilation_trans_biomart_up_validation$Row.names,]
gene_FPKM_valid_down <- gene_FPKM[rownames(gene_FPKM) %in% res_sig_DESeq_FC1_compilation_trans_biomart_down_validation$Row.names,]

gene_FPKM_valid_up <- merge(gene_FPKM_valid_up, res_sig_DESeq_FC1_compilation_trans_biomart_up_validation[,c("Row.names","entrezgene_accession")], by.x=0, by.y="Row.names")
gene_FPKM_valid_down <- merge(gene_FPKM_valid_down, res_sig_DESeq_FC1_compilation_trans_biomart_down_validation[,c("Row.names","entrezgene_accession")], by.x=0, by.y="Row.names")

gene_FPKM_valid_down_final <- gene_FPKM_valid_down[gene_FPKM_valid_down$entrezgene_accession %!in% "VEGFB",]
gene_FPKM_valid_up_final<- gene_FPKM_valid_up[gene_FPKM_valid_up$entrezgene_accession %in% c("IL1R1", "CXCL2", "NFKB1", "BIRC3", "IL6", "CCL2", "PLAU", "NFKB2", "CASP7",  "PLAT"),]




gtf_stringtie_zator_transcript_up_final <- gtf_stringtie_zator_transcript[gtf_stringtie_zator_transcript$gene_id %in% gene_FPKM_valid_up_final$Row.names, ]
gtf_stringtie_zator_transcript_down_final <- gtf_stringtie_zator_transcript[gtf_stringtie_zator_transcript$gene_id %in% gene_FPKM_valid_down_final$Row.names, ]
transcript_count_matrix <- read.csv2("transcript_count_matrix.csv", sep = ",", row.names = 1)
head(transcript_count_matrix)

transcript_count_matrix_up_final_valid <- transcript_count_matrix[rownames(transcript_count_matrix) %in% gtf_stringtie_zator_transcript_up_final$transcript_id, ]
transcript_count_matrix_up_final_valid_trim <- transcript_count_matrix_up_final_valid[,colnames(transcript_count_matrix_up_final_valid) %in% colnames(gene_FPKM_valid_up_final)]
transcript_count_matrix_up_final_valid_trim <- merge(transcript_count_matrix_up_final_valid_trim, gtf_stringtie_zator_transcript_up_final[,c("transcript_id","gene_name","gene_id"),], by.x=0, by.y="transcript_id" )

transcript_count_matrix_down_final_valid <- transcript_count_matrix[rownames(transcript_count_matrix) %in% gtf_stringtie_zator_transcript_down_final$transcript_id, ]
transcript_count_matrix_down_final_valid_trim <- transcript_count_matrix_down_final_valid[,colnames(transcript_count_matrix_down_final_valid) %in% colnames(gene_FPKM_valid_down_final)]
transcript_count_matrix_down_final_valid_trim <- merge(transcript_count_matrix_down_final_valid_trim, gtf_stringtie_zator_transcript_down_final[,c("transcript_id","gene_name", "gene_id"),], by.x=0, by.y="transcript_id" )

#####trans corelation
gene_FPKM_all = as.data.frame(gexpr(bg))
log_FPKM <- log2(gene_FPKM_all + 1)

log_FPKM_sig_protein <- log_FPKM[rownames(log_FPKM) %in% res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$Row.names,]
sim_matrix <- cordist(log_FPKM_sig_protein)

sim_matrix2<- sim_matrix[rownames(sim_matrix) %in% res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$Row.names,] #for genes

sim_matrix2 <- sim_matrix2[,colnames(sim_matrix2) %in% res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$Row.names]
sim_matrix3 = melt(as.matrix(sim_matrix2))

sim_matrix3<-sim_matrix3[order(sim_matrix3$Var1),]

sim_matrix5<- sim_matrix3[sim_matrix3$value > 0.9,]
sim_matrix5_unique <- sim_matrix5[sim_matrix5$value %!in% 1.0,]

res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein <- res_sig_DESeq_FC1_compilation_trans_biomart[res_sig_DESeq_FC1_compilation_trans_biomart$gene_biotype %in% "protein_coding",]
res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein <- res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein[!duplicated(res_sig_DESeq_FC1_compilation_trans_biomart_circos_protein$Row.names),]





#------------------------------------------


