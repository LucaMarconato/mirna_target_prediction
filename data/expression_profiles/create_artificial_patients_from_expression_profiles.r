## This is a tempory conversion that I make, in the future I will modify the C++ code so to accept expression profiles in a more general format
library(base)

mirna_id_dictionary <- read.table("../processed/mirna_id_dictionary.tsv", header = T, colClasses = c("character", "numeric"))
gene_id_dictionary <- read.table("../processed/gene_id_dictionary.tsv", header = T, colClasses = c("character", "numeric"))

expression_profiles <- read.table("expression_profiles.tsv", header = T, colClasses = c("character", "character", "character"))
for(i in 1:nrow(expression_profiles)) {
    mirna_dataset <- expression_profiles$mirna_dataset[i]
    gene_dataset <- expression_profiles$gene_dataset[i]

    patient_folder <- paste("../patients/artificial_", expression_profiles$artificial_patient_id[i], sep = "")
    if(dir.exists(patient_folder)) {
        print(paste(patient_folder, "already exists, skipping it"))
    } else {
        mirna_expression_profile <- read.table(mirna_dataset, header = T, colClasses = c("character", "numeric", "numeric"))
        gene_expression_profile <- read.table(gene_dataset, header = T, colClasses = c("character", "numeric"))
        gene_expression_profile$gene_id <- unlist(lapply(gene_expression_profile$gene_id, function(x) strsplit(x, "\\.")[[1]][1]))

        mirna_expression_profile$mirbase_id <- tolower(mirna_expression_profile$mirbase_id)
        mirna_indices <- match(mirna_expression_profile$mirbase_id, mirna_id_dictionary$mirna_family)
        gene_indices <- match(gene_expression_profile$gene_id, gene_id_dictionary$gene_id)

        mirna_expression_profile$mirna_id_cpp <- mirna_id_dictionary$mirna_id_cpp[mirna_indices]
        mirna_not_recognized_count <- sum(is.na(mirna_expression_profile$mirna_id_cpp))
        gene_expression_profile$gene_id_cpp <- gene_id_dictionary$gene_id_cpp[gene_indices]
        gene_not_recognized_count <- sum(is.na(gene_expression_profile$gene_id_cpp))

        mirna_expression_profile_path <- paste(patient_folder, "/mirna_expression_profile.tsv", sep = "")
        gene_expression_profile_path <- paste(patient_folder, "/gene_expression_profile.tsv", sep = "")

        lost_mirna_reads <- sum(mirna_expression_profile$reads[is.na(mirna_expression_profile$mirna_id_cpp)])
        total_mirna_reads <- sum(mirna_expression_profile$reads)
        lost_mirna_ratio <- lost_mirna_reads/total_mirna_reads
        lost_gene_reads <- sum(gene_expression_profile$reads[is.na(gene_expression_profile$gene_id_cpp)])
        total_gene_reads <- sum(gene_expression_profile$reads)
        lost_gene_ratio <- lost_gene_reads/total_gene_reads

        print(paste(mirna_not_recognized_count, "/", length(mirna_expression_profile$mirna_id_cpp),
                    " mirnas not recognized (", length(mirna_expression_profile[[1]]) - mirna_not_recognized_count,
                    " recognized), losing ", lost_mirna_reads, "/", total_mirna_reads,
                    " (=", round(lost_mirna_ratio, 2), ") reads",
                    sep = ""))
        print(paste(gene_not_recognized_count, "/", length(gene_expression_profile$gene_id_cpp),
                    " genes not recognized, (", length(gene_expression_profile[[1]]) - gene_not_recognized_count,
                    " recognized), losing ", lost_gene_reads, "/", total_gene_reads,
                    " (=", round(lost_gene_ratio, 2), ") reads",
                    sep = ""))

        mirna_not_recognized_path <- paste(patient_folder, "/mirna_not_recognized.tsv", sep = "")
        mirna_recognized_path <- paste(patient_folder, "/mirna_expression_profile.tsv", sep = "")
        gene_not_recognized_path <- paste(patient_folder, "/gene_not_recognized.tsv", sep = "")
        gene_recognized_path <- paste(patient_folder, "/gene_expression_profile.tsv", sep = "")

        mirna_not_recognized <- mirna_expression_profile[is.na(mirna_expression_profile$mirna_id_cpp), ]
        mirna_not_recognized <- base::subset(mirna_not_recognized, select = c("mirbase_id", "reads"))
        mirna_recognized <- mirna_expression_profile[!is.na(mirna_expression_profile$mirna_id_cpp), ]
        mirna_recognized <- base::subset(mirna_recognized, select = c("mirbase_id", "reads", "rpm"))
        gene_not_recognized <- gene_expression_profile[is.na(gene_expression_profile$gene_id_cpp), ]
        gene_not_recognized <- base::subset(gene_not_recognized, select = c("gene_id", "reads"))
        colnames(gene_not_recognized)[1] <- "gene_id"
        gene_recognized <- gene_expression_profile[!is.na(gene_expression_profile$gene_id_cpp), ]
        gene_recognized <- base::subset(gene_recognized, select = c("gene_id", "reads"))
        colnames(gene_recognized)[1] <- "gene_id"

        dir.create(patient_folder)

        write.table(mirna_not_recognized, mirna_not_recognized_path, row.names = F, quote = F, sep = "\t")
        write.table(mirna_recognized, mirna_recognized_path, row.names = F, quote = F, sep = "\t")
        write.table(gene_not_recognized, gene_not_recognized_path, row.names = F, quote = F, sep = "\t")
        write.table(gene_recognized, gene_recognized_path, row.names = F, quote = F, sep = "\t")
        ## browser()
        ## browser()
    }
}
