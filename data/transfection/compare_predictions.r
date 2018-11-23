transfection_data <- read.table("targetscan_data/transfection_data.tsv", stringsAsFactors = F, header = T, fill = T)

patients <- c("artificial_ENCFF360IHM-hela", "artificial_ENCFF495ZXC-hela", "artificial_ENCFF612ZIR-hela", "artificial_ENCFF729EQX-hela", "artificial_ENCFF806EYY-hela", "artificial_ENCFF902KUU-hela", "artificial_TCGA-CJ-4642")

gene_id_dictionary <- read.table("../processed/gene_id_dictionary.tsv", header = T, colClasses = c("character", "numeric"))

for(patient in patients) {
    print(paste("patient = ", patient, sep = ""))
    path <- paste("../patients/", patient, "/matchings_predictor_output/original_data/predicted_downregulation/p_j_downregulated_values_999999.tsv", sep = "")
    predictions <- read.table(path, header = T, colClasses = c("numeric", "numeric"))
    ## TODO: the column gene_id in the file p_j_downregulated_values_xxxxxx.tsv should be called gene_id_cpp
    ## TODO: the column gene_id in the file gene_id_dictionary.tsv should be called ensembl_gene_id
    predictions$ensembl_gene_id <- gene_id_dictionary$gene_id[match(predictions$gene_id, gene_id_dictionary$gene_id_cpp)]
    na_count <- sum(is.na(predictions$ensembl_id))
    if(na_count > 0) {
        print(paste("na_count = ", na_count, sep = ""))
        stop("abort")
    }
    predictions$fold_change <- 1 - predictions$p_j_downregulated_values
    predictions$log2_fold_change <- log2(predictions$fold_change)

    common_ensembl_gene_ids <- intersect(transfection_data$ensembl_gene_id, predictions$ensembl_gene_id)
    ids_missing_in_targetscan_data_count <- length(setdiff(transfection_data$ensembl_gene_id, predictions$ensembl_gene_id))
    ids_missing_in_patient_data_count <- length(setdiff(predictions$ensembl_gene_id, transfection_data$ensembl_gene_id))
    print(paste(ids_missing_in_targetscan_data_count, "Ensembl ids are not present in TargetScan data"))
    print(paste(ids_missing_in_patient_data_count, "Ensembl ids are not present in the patient data"))
    a <- transfection_data[match(common_ensembl_gene_ids, transfection_data$ensembl_gene_id), ]
    b <- predictions[match(common_ensembl_gene_ids, predictions$ensembl_gene_id), ]
    browser()
}
