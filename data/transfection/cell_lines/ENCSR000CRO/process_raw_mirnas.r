library(data.table)
library(rjson)

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

    recognized_first_strand <- a[indexes_of_recognized, c("gene_id", "first_strand_count")]
    recognized_second_strand <- a[indexes_of_recognized, c("gene_id", "first_strand_count")]

    recognized_first_strand$gene_id <- paste(ensembl_ids_mirbase_ids_dictionary$mirbase_id[match(recognized_first_strand$gene_id, ensembl_ids_mirbase_ids_dictionary$ensembl_id)], "-3p", sep = "")
    recognized_second_strand$gene_id <- paste(ensembl_ids_mirbase_ids_dictionary$mirbase_id[match(recognized_second_strand$gene_id, ensembl_ids_mirbase_ids_dictionary$ensembl_id)], "-5p", sep = "")

    mirna_expression_profile <- data.table(mirbase_id = c(
                                               recognized_first_strand$gene_id,
                                               recognized_second_strand$gene_id
                                           ), reads = c(
                                                  recognized_first_strand$first_strand_count,
                                                  recognized_second_strand$second_strand_count
                                              ), stringsAsFactors = F)

    mirna_expression_profile$rpm <- mirna_expression_profile$reads/sum(mirna_expression_profile$reads)*1000000
    to_filter <- mirna_expression_profile$rpm < mirna_threshold_rpm
    print(paste("filtering ", sum(to_filter), "/", nrow(mirna_expression_profile), " mirnas, having an rpm value of less than ", mirna_threshold_rpm, sep = ""))
    rpm_filtered <- sum(mirna_expression_profile[to_filter, "rpm"])
    print(paste("total rpms filtered: ", round(rpm_filtered, 2),
                " (=", round(rpm_filtered/1000000, 2), "%)",
                sep = ""))
    mirna_expression_profile <- mirna_expression_profile[!to_filter, ]
    plot(sort(mirna_expression_profile$reads))
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
