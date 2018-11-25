library(rjson)

patients <- c("artificial_ENCFF360IHM-hela", "artificial_ENCFF495ZXC-hela", "artificial_ENCFF612ZIR-hela", "artificial_ENCFF729EQX-hela", "artificial_ENCFF806EYY-hela", "artificial_ENCFF902KUU-hela", "artificial_TCGA-CJ-4642")
transfection_experiments <- fromJSON(file = "targetscan_data/info.json")

patients_folder <- "../patients/"
for(patient in patients) {
    for(dataset in transfection_experiments$datasets) {
        if(dataset$mirbase_id != "") {
            dataset$mirbase_id <- tolower(dataset$mirbase_id)
            new_patient_name <- paste(patient, "_", dataset$mirbase_id, "_transfected", sep = "")
            new_patient_folder_path <- paste(patients_folder, new_patient_name, "/", sep = "")
            old_patient_folder_path <- paste(patients_folder, patient, "/", sep = "")
            if(!dir.exists(new_patient_folder_path)) {
                dir.create(new_patient_folder_path)
            }
            file.copy(from = paste(old_patient_folder_path, "gene_expression_profile.tsv", sep = ""),
                      to = paste(new_patient_folder_path, "gene_expression_profile.tsv", sep = ""),
                      overwrite = T)
            old_mirna_expression_profile <- read.table(paste(old_patient_folder_path, "mirna_expression_profile.tsv", sep = ""),
                                                       header = T,
                                                       colClasses = c("character", "numeric", "numeric"))
            index <- match(dataset$mirbase_id, old_mirna_expression_profile$mirbase_id)
            transfected_reads <- 500000
            if(is.na(index)) {
                old_reads <- 0
                ## I will remove rpm from this file since they are useless
                new_mirna_expression_profile <- rbind(old_mirna_expression_profile,
                                                      data.frame(mirbase_id = dataset$mirbase_id, reads = transfected_reads, rpm = -10000))
            } else {
                old_reads <- old_mirna_expression_profile$reads[index]
                new_mirna_expression_profile <- old_mirna_expression_profile
                new_mirna_expression_profile$reads[index] <- old_reads + transfected_reads
            }
            print(paste("adding ", transfected_reads, " to miRNA ", dataset$mirbase_id, ", which previously had ", old_reads, sep = ""))
            write.table(new_mirna_expression_profile, file = paste(new_patient_folder_path, "mirna_expression_profile.tsv", sep = ""),
                        row.names = F, quote = F, sep = "\t")
        }
    }
}

print("now run the C++ simulations; I will create a pipeline with Snakemake for making the code reproducible")
