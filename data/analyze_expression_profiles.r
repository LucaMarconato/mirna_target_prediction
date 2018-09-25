library(hashmap)

mirnas <- read.table("patients/TCGA-CJ-4642/tumor_mirna_expression_profile.tsv", header = T, colClasses = c("numeric", "numeric"))
top10_mirnas <- mirnas[order(-mirnas$rpm),][1:10,]

genes <- read.table("patients/TCGA-CJ-4642/tumor_gene_expression_profile.tsv", header = T, colClasses = c("numeric", "numeric"))
gene_id_dictionary <- read.table("processed/gene_id_dictionary.tsv", header = T, colClasses = c("character", "numeric"))
gene_map <- hashmap(keys = gene_id_dictionary$gene_id_cpp, values = gene_id_dictionary$gene_id)

interactions <- read.table("processed/scored_interactions_processed.tsv", colClasses = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"), header = T)

gene_expression_profile <- read.table("patients/TCGA-CJ-4642/481383ce-b91b-4051-be71-0742f26cb178/d9c7a1f7-f303-4617-b46c-2003e871095b.htseq.counts", colClasses = c("character", "numeric"))
gene_expression_profile[[1]] <- unlist(lapply(gene_expression_profile[[1]], function(x) strsplit(x, "\\.")[[1]][1]))

total_reads <- sum(gene_expression_profile[[2]])

for(mirna_id in top10_mirnas$mirna_id) {
    interacting_genes <- interactions[interactions$mirna_id == mirna_id, "gene_id"]
    ensembl_ids <- gene_map[[interacting_genes]]
    targets <- gene_expression_profile[gene_expression_profile[[1]] %in% ensembl_ids, ]
    print(paste(dim(targets)[[1]], " protein coding genes interacting with the mirna with mirna_id = ", mirna_id, sep = ""))
    target_reads <- sum(targets[[2]])
    print(paste("target relative fraction = ", round(target_reads/total_reads*100,2), "%", sep = ""))
}

## the output of this code is consistent with the results of the simulation: mirnas with more available sites (considering expression profiles) get binded more
