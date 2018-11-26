transfection_data <- read.table("targetscan_data/transfection_data.tsv", stringsAsFactors = F, header = T, fill = T)

mirna_id_dictionary <- read.table("../processed/mirna_id_dictionary.tsv", header = T, colClasses = c("character", "numeric"))
gene_id_dictionary <- read.table("../processed/gene_id_dictionary.tsv", header = T, colClasses = c("character", "numeric"))

## load only if necessary since the process is slow
if(!exists("interactions")) {
    interactions <- read.table("../processed/scored_interactions_processed.tsv", colClasses = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"), header = T)
}

hela_patients <- c("artificial_ENCFF360IHM-hela", "artificial_ENCFF495ZXC-hela", "artificial_ENCFF612ZIR-hela", "artificial_ENCFF729EQX-hela", "artificial_ENCFF806EYY-hela", "artificial_ENCFF902KUU-hela")
transfected_mirnas <- c("hsa-mir-122-5p", "hsa-mir-128-3p", "hsa-mir-132-3p", "hsa-mir-142-5p")
corresponding_experiments <- c("GSM210901", "GSM210903", "GSM210904", "GSM210909")

correlations <- data.frame()

for(hela_patient in hela_patients) {
    i <- 1
    for(i in 1:length(transfected_mirnas)) {
        transfected_mirna <- transfected_mirnas[[i]]
        corresponding_experiment <- corresponding_experiments[[i]]
        parse_predictions <- function(patient) {
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
            return(predictions)
        }
        patient_not_transfected <- hela_patient
        patient_transfected <- paste(hela_patient, "_", transfected_mirna, "_transfected", sep = "")
        predictions_not_transfected <- parse_predictions(patient_not_transfected)
        predictions_transfected <- parse_predictions(patient_transfected)

        common_ensembl_gene_ids <- intersect(transfection_data$ensembl_gene_id, predictions_not_transfected$ensembl_gene_id)
        common_ensembl_gene_ids <- intersect(common_ensembl_gene_ids, predictions_transfected$ensembl_gene_id)
        ids_missing_in_targetscan_data_count <- length(setdiff(transfection_data$ensembl_gene_id, common_ensembl_gene_ids))
        ids_missing_in_the_not_transfected_patient_data_count <- length(setdiff(predictions_not_transfected$ensembl_gene_id, common_ensembl_gene_ids))
        ids_missing_in_the_transfected_patient_data_count <- length(setdiff(predictions_transfected$ensembl_gene_id, common_ensembl_gene_ids))

        print(paste(ids_missing_in_targetscan_data_count, "Ensembl ids are not present in TargetScan data"))
        print(paste(ids_missing_in_the_not_transfected_patient_data_count, "Ensembl ids are not present in the patient data (not transfected one)"))
        print(paste(ids_missing_in_the_transfected_patient_data_count, "Ensembl ids are not present in the patient data (transfected one)"))

        a <- transfection_data[match(common_ensembl_gene_ids, transfection_data$ensembl_gene_id),
                               match(c("refseq_id",
                                       "used_in_training",
                                       "gene_symbol",
                                       "predicted_gene_symbol",
                                       "ensembl_gene_id",
                                       corresponding_experiment),
                                     colnames(transfection_data))]
        b <- predictions_not_transfected[match(common_ensembl_gene_ids, predictions_not_transfected$ensembl_gene_id), ]
        d <- predictions_transfected[match(common_ensembl_gene_ids, predictions_transfected$ensembl_gene_id), ]

        a <- a[order(a$ensembl_gene_id), ]
        colnames(a)[[match(corresponding_experiment, colnames(a))]] <- "log2_fold_change"
        b <- b[order(b$ensembl_gene_id), ]
        d <- d[order(d$ensembl_gene_id), ]

        rows_without_na <- !is.na(a$log2_fold_change)
        a <- a[rows_without_na, ]
        b <- b[rows_without_na, ]
        d <- d[rows_without_na, ]
        print(paste("removed ", sum(!rows_without_na), "/", sum(rows_without_na) + sum(!rows_without_na),
                    " rows which contained NAs in the TargetScan data", sep = ""))
        if(sum(is.na(b$ensembl_gene_id)) || sum(is.na(d$ensembl_gene_id))) {
            print("error: found NAs in predictions")
            stop("abort")
        }

        mirna_id_cpp <- mirna_id_dictionary$mirna_id_cpp[mirna_id_dictionary$mirna_family == transfected_mirna]
        filename <- paste("../interactions/mirna_id", mirna_id_cpp, "_interactions.rds", sep = "")
        if(!file.exists(filename)) {
            ## I should have used a sql database here...
            print(paste("generating target file for mirna_id_cpp =", mirna_id_cpp))
            source("analyze_expression_profiles.r")
            get_targets_for_mirna(mirna_id_cpp)
            print("generated")
        }
        gene_ids_cpp_for_mirna <- readRDS(filename)
        ensembl_ids_for_mirna <- gene_id_dictionary$gene_id[match(gene_ids_cpp_for_mirna, gene_id_dictionary$gene_id_cpp)]

        a <- a[a$ensembl_gene_id %in% ensembl_ids_for_mirna, ]
        b <- b[b$ensembl_gene_id %in% ensembl_ids_for_mirna, ]
        d <- d[d$ensembl_gene_id %in% ensembl_ids_for_mirna, ]

        la <- a$log2_fold_change
        lb <- b$log2_fold_change
        ld <- d$log2_fold_change
        ld_lb <- ld - lb

        current_interactions <- interactions[interactions$mirna_id == mirna_id_cpp, c("mirna_id", "gene_id", "context_score", "weighted_context_score", "conserved")]
        current_interactions$ensembl_gene_id <- gene_id_dictionary$gene_id[match(current_interactions$gene_id, gene_id_dictionary$gene_id_cpp)]
        context_scores <- current_interactions$context_score[match(a$ensembl_gene_id, current_interactions$ensembl_gene_id)]
        weighted_context_scores <- current_interactions$weighted_context_score[match(a$ensembl_gene_id, current_interactions$ensembl_gene_id)]

        le <- context_scores
        lf <- weighted_context_scores
        all_the_predictions <- data.frame(la = la, lb = lb, ld = ld, le = le, lf = lf, ld_lb = ld - lb)

        training_set <- a$used_in_training == T
        validation_set <- a$used_in_training == F
        la_t <- la[training_set]
        la_v <- la[validation_set]

        lb_t <- lb[training_set]
        lb_v <- lb[validation_set]

        ld_t <- ld[training_set]
        ld_v <- ld[validation_set]

        ld_lb_t <- ld_lb[training_set]
        ld_lb_v <- ld_lb[validation_set]

        le_t <- le[training_set]
        le_v <- le[validation_set]

        lf_t <- lf[training_set]
        lf_v <- lf[validation_set]

        correlations <- rbind(correlations, data.frame(b = cor(la_v, lb_v, method = "spearman"),
                                                       d = cor(la_v, ld_v, method = "spearman"),
                                                       d_b = cor(la_v, ld_lb_v, method = "spearman"),
                                                       e = cor(la_v, le_v, method = "spearman"),
                                                       f = cor(la_v, lf_v, method = "spearman")))

        print(paste("--->Correlation expected to be low:", round(cor(la_v, lb_v, method = "spearman"), 2), "<---", sep = ""))
        print(paste("--->Correlation expected to be low:", round(cor(la_v, ld_v, method = "spearman"), 2), "<---", sep = ""))
        print(paste("--->TargetScan predictions correlation:", round(cor(la_v, le_v, method = "spearman"), 2), "<---", sep = ""))
        print(paste("--->New method predictions correlation:", round(cor(la_v, ld_lb_v, method = "spearman"), 2), "<---", sep = ""))
        pairs(all_the_predictions)
        ## browser()
        i <- i + 1
    }
}
