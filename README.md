# VTE
This repository contains the scripts for the "venous thromboembolism and transcriptomic analysis of deep vein thrombosis in the femoral veins" project ([ENA: PRJEB43020](https://www.ebi.ac.uk/ena/browser/view/PRJEB43020)).

<img align="middle" src="img/abs.jpg" alt="graphical abstract" wigth="300px" />

### Abstract 1
Venous thromboembolism (VTE), including deep vein thrombosis (DVT) and pulmonary embolism (PE), is a severe disease affecting the human venous system, accompanied by high morbidity and mortality rates. The aim of the study was to establish a new porcine VTE model based on the formation of the thrombus in vivo. The study was performed on 10 castrated male pigs: thrombus was formed in each closed femoral vein and then successfully released from the right femoral vein into the circulation of animals. In six pigs PE was confirmed via both computed tomography pulmonary angiography and an autopsy. Our research presents a novel experimental porcine model of VTE that involves inducing DVT and PE in the same animal in vivo, making it suitable for advanced clinical research and testing of future therapies.
* Gromadziński, L.; Skowrońska, A.; Holak, P.; Smoliński, M.; Lepiarczyk, E.; Żurada, A.; Majewski, M.K.; Skowroński, M.T.; Majewska, M. A New Experimental Porcine Model of Venous Thromboembolism. J. Clin. Med. 2021, 10, 1862. https://doi.org/10.3390/jcm10091862

---------------------------------------------------------------

## Script

The script consist of a set of files, where the main file is named DVT.R (Deep Vein Thrombosis). This file contains:
### reading the data
```
# setting up samples - pheno data
head(pheno_data)
           IDs Condition Day Object Group
1 PBS_1d_W2_52       PBS  1d     W2     A
2 PBS_1d_W5_55       PBS  1d     W5     A
3 PBS_1d_W8_58       PBS  1d     W8     A
4 PBS_1d_W9_59       PBS  1d     W9     A
5 MC2_1d_W3_43       MC2  1d     W3     B
6 MC2_1d_W6_46       MC2  1d     W6     B 
```
(...)
```
# expression data
count_table <- read.csv2("gene_count_matrix.csv", sep=",", row.names = 1)
```
(...)
```
# reading the annotation fles
gtf_stringtie <- import.gff("stringtie/stringtie_merged.gtf")

# extend annotation using BIOMART
ensemblSs = useMart(host="ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="sscrofa_gene_ensembl") # robię bazę GI i ENSEMBL i Symboli genów
allgenes.Ensembl = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype", "entrezgene_id", "description", "entrezgene_accession"),mart=ensemblSs)
```

### DE analyses
```
# ballgown
bg <- ballgown(dataDir = "ballgown_AvsC/", samplePattern=".", pData=pheno_data_AvsC)
bg_filt <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset=TRUE)
res_sig <- subset(results_genes, results_genes$qval <= 0.05)
res_sig_FC1 <- subset(res_sig, abs(res_sig$log2FC) > 1)

# DESeq2
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
```
 ### other analyses
 * validation
 * correlation

### other files
* libs.R - used external libraries
* funs.R - used own functions
* fig2.R - preparing figure 2
* fig3.R - preparing figure 3
* fig4.R - preparing figure 4
* fig5.R - preparing figure 5

