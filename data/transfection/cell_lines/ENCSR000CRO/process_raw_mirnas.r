library(data.table)
library(base)
library(rjson)
source("../../../my_device.r")

process_file <- function(file, mirna_threshold_rpm)
{
    a <- read.table(paste("raw/", file, sep = ""), colClasses = c("character", "numeric", "numeric", "numeric"))
    colnames(a) <- c("gene_id", "unstranded_count", "first_strand_count", "second_strand_count")

    a <- subset(a, grepl("ENSG", gene_id))
    ## TODO: consider the version
    a$gene_id <- unlist(lapply(a$gene_id, function(x) strsplit(x, "\\.")[[1]][1]))
    ensembl_ids_mirbase_ids_dictionary <- read.table("../../../processed/ensembl_id_to_mirbase_id.tsv", header = T, colClasses = c("character", "character"))

    x <- a$gene_id
    y <- ensembl_ids_mirbase_ids_dictionary$ensembl_id
    ensembl_ids_recognized <- intersect(x, y)
    ensembl_ids_not_recognized <- setdiff(union(x, y), intersect(x, y))
    ## WARNING: note that different ensembl ids can map to the same mirna
    print(paste(length(ensembl_ids_recognized), "ensembl ids recognized as mirnas"))
    print(paste(length(ensembl_ids_not_recognized), "ensembl ids not recognized as mirnas"))
    indexes_of_recognized <- match(ensembl_ids_recognized, a$gene_id)

    recognized_first_strand <- a[indexes_of_recognized, c("gene_id", "first_strand_count", "second_strand_count")]
    recognized_second_strand <- a[indexes_of_recognized, c("gene_id", "first_strand_count", "second_strand_count")]

    distinguish_isoforms <- function(mirnas, suffix)
    {
        to_return <- mirnas
        ## one could load these data only once, but the code is fast anyway
        if(!exists("mirna_isoforms_dictionary")) {
            print("loading data")
            mirna_isoforms_dictionary <<- read.table("../../../processed/mirna_isoforms_dictionary.tsv", header = T, colClasses = c("character", "character", "character"))
            mirna_isoforms_table <<- table(mirna_isoforms_dictionary$ambiguous_mirna_id)
        }

        for(i in seq_along(mirnas$mirna_id)) {
            ambiguous_mirna_id <- mirnas$mirna_id[i]
            if(mirna_isoforms_table[ambiguous_mirna_id] == 1) {
                disambiguated_mirna_id <- mirna_isoforms_dictionary[mirna_isoforms_dictionary$ambiguous_mirna_id == ambiguous_mirna_id, "disambiguated_mirna_id"]
                if(mirnas[i, "first_strand_count"] > 0 && mirnas[i, "second_strand_count"] > 0) {
                    print(paste("error: mirnas[", i, ", \"first_strand_count\"] = ",
                               mirnas[i, "first_strand_count"], ", ",
                               "mirnas[", i, ", \"second_strand_count\"] = ",
                               mirnas[i, "second_strand_count"],
                               ", do not ignore this error",
                               sep = ""))
                }
                to_return[i, "mirna_id"] <- disambiguated_mirna_id
            } else if(mirna_isoforms_table[ambiguous_mirna_id] == 2) {
                rows <- mirna_isoforms_dictionary[mirna_isoforms_dictionary$ambiguous_mirna_id == ambiguous_mirna_id, ]
                valid_disambiguations_found <- 0
                for(row_index in c(1, 2)) {
                    disambiguated_mirna_id <- rows$disambiguated_mirna_id[row_index]
                    if(base::endsWith(disambiguated_mirna_id, suffix)) {
                        to_return[i, "mirna_id"] <- disambiguated_mirna_id
                        valid_disambiguations_found = valid_disambiguations_found + 1
                    }
                }
                if(valid_disambiguations_found != 1) {
                    stop(paste("error: valid_disambiguations_found = ", valid_disambiguations_found, sep = ""))
                }
            } else {
                stop(paste("error: mirna_isoforms_table[", ambiguous_mirna_id, "] = ", mirna_isoforms_table[ambiguous_mirna_id], sep = ""))
            }
        }
        return(to_return)
    }

    recognized_first_strand$gene_id <- ensembl_ids_mirbase_ids_dictionary$mirbase_id[match(recognized_first_strand$gene_id, ensembl_ids_mirbase_ids_dictionary$ensembl_id)]
    recognized_second_strand$gene_id <- ensembl_ids_mirbase_ids_dictionary$mirbase_id[match(recognized_second_strand$gene_id, ensembl_ids_mirbase_ids_dictionary$ensembl_id)]

    colnames(recognized_first_strand)[[1]] <- "mirna_id"
    colnames(recognized_second_strand)[[1]] <- "mirna_id"

    recognized_first_strand <- distinguish_isoforms(recognized_first_strand, "3p")
    recognized_second_strand <- distinguish_isoforms(recognized_second_strand, "5p")

    recognized_first_strand <- recognized_first_strand[recognized_first_strand$mirna_id != "", ]
    recognized_second_strand <- recognized_second_strand[!is.null(recognized_second_strand$mirna_id)]

    mirna_expression_profile <- data.table(mirbase_id = c(
                                               recognized_first_strand$mirna_id,
                                               recognized_second_strand$mirna_id
                                           ), reads = c(
                                                  recognized_first_strand$first_strand_count,
                                                  recognized_second_strand$second_strand_count
                                              ), stringsAsFactors = F)

    aggregate(reads ~ mirbase_id, data = mirna_expression_profile, FUN = sum)
    mirna_expression_profile$rpm <- mirna_expression_profile$reads/sum(mirna_expression_profile$reads)*1000000
    write.table(mirna_expression_profile, file = paste("processed/", file, sep = ""), sep = "\t", row.names = F, quote = F)

    to_filter <- mirna_expression_profile$rpm < mirna_threshold_rpm
    print(paste("filtering ", sum(to_filter), "/", nrow(mirna_expression_profile),
                " mirnas, having an rpm value of less than ", mirna_threshold_rpm,
                " (", nrow(mirna_expression_profile) - sum(to_filter), " remaining)",
                sep = ""))
    rpm_filtered <- sum(mirna_expression_profile[to_filter, "rpm"])
    print(paste("total rpms filtered: ", round(rpm_filtered, 2),
                " (=", round(rpm_filtered/1000000, 4), "%)",
                sep = ""))
    mirna_expression_profile <- mirna_expression_profile[!to_filter, ]
    ## plot(sort(mirna_expression_profile$reads))
}

analyze_correlation <- function(files, mirna_threshold_rpm)
{
    dataframes <- lapply(files, function(file) read.table(paste("processed/", file, sep = ""), header = T, colClasses = c("character", "numeric", "numeric")))

    mirnas_involved <- unique(unlist(lapply(dataframes, function(x) x$mirbase_id)))
    ## TODO: find the mirnas missing for each dataset

    ## TODO: plot the correlation matrix
    for(df in dataframes) {
        
    }
    browser()
    browser()
}

info_json <- fromJSON(file = "../../../../global_parameters.json")
mirna_threshold_rpm <- info_json["mirna_threshold_rpm"]

mirna_cpp_id_dictionary <- read.table("../../../processed/mirna_id_dictionary.tsv", header = T, colClasses = c("character", "numeric"))
targetscan_mirnas <- mirna_cpp_id_dictionary[[1]]

files <- c("ENCFF495ZXC.tsv", "ENCFF902KUU.tsv")
for(file in files) {
    process_file(file, mirna_threshold_rpm)
}

analyze_correlation(files, mirna_threshold_rpm)
