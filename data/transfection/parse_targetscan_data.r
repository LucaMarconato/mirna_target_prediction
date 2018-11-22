## to install biomaRt do not use install.packages() but use the following procedure
## 1) install Bioconductor by using the following code
##    if (!requireNamespace("BiocManager"))
##        install.packages("BiocManager")
##    BiocManager::install()
## 2) install biomaRt using Bioconductor
##    BiocManager::install(c("biomaRt"))
library(biomaRt)

a <- read.table("targetscan_data/elife-05005-supp1-v1.csv", sep = ";", stringsAsFactors = F)
a <- a[3:nrow(a), ]
column_names <- as.character(a[1, ])
a <- a[2:nrow(a), ]
colnames(a) <- column_names
colnames(a)[[1]] <- "refseq_id"
colnames(a)[[2]] <- "used_in_training"
colnames(a)[[3]] <- "gene_symbol"
a$gene_symbol <- toupper(a$gene_symbol)
a$used_in_training <- as.factor(a$used_in_training)

as.num = function(x, na.strings = "NA") {
    stopifnot(is.character(x))
    na <- x %in% na.strings
    ## browser()
    x[na] <- 0
    x <- as.numeric(x)
    x[na] <- NA_real_
    return(x)
}

for(i in 4:ncol(a)) {
    comma_replaced <- sub(",", ".", a[[i]], fixed = T)
    a[[i]] <-  as.num(comma_replaced, na.strings = "N/A")
}

describe_a <- function()
{
    stopifnot(dim(a)[[2]] == 77)
    print(paste("the matrix has", dim(a)[[1]], "rows and", dim(a)[[2]], "columns"))
    print("there is one row for each transcript of each human protein-coding gene (many different transcripts can be copied from the same gene)")
    print("the first 3 column names are refseq_id, used_in_training and gene_symbol")
    print("the remaining 74 columns refers to the 74 different transfection experiments")
    print("the value of a generic cell belonging to the last 74 columns is the log-fold change (log2) of the transcript expression after the transfection of the miRNA (sometimes a siRNA) specified in the column name")
    print("such a value is NA if the expression values of the transcript before and after the transfection were too low to accurately quantify the log-fold change (these expression values are not shown)")

    used_in_training_count <- sum(a$used_in_training == "yes")
    not_used_in_training_count <- sum(a$used_in_training == "no")
    print(paste(used_in_training_count, "rows have been used in training"))
    print(paste(not_used_in_training_count, "rows have not been used in training"))

    row_has_na <- apply(a, 1, function(x) sum(is.na(x)) > 0)
    rows_containing_na_count <- sum(row_has_na == T)
    rows_not_containing_na_count <- sum(row_has_na == F)
    print(paste(rows_containing_na_count, "rows have at least 1 NA"))
    print(paste(rows_not_containing_na_count, "rows do not have NAs"))
    few_na <- 5
    row_has_few_na <- apply(a, 1, function(x) sum(is.na(x)) < few_na)
    rows_with_few_na_count <- sum(row_has_few_na == T)
    print(paste(rows_with_few_na_count, "rows with less than", few_na, "NAs"))
}

describe_a()

## for the moment we will not consider different transcript, but for each gene we will choose one representative transcript. For this reason we will now find, for each transcript, the Ensembl id for the corresponding gene
if(!file.exists("refseq_embl_conversion_table.tsv")) {
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    refseq_embl_conversion_table <- getBM(attributes = c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), filters = "refseq_mrna", values = a$refseq_id, mart = ensembl)
    refseq_embl_conversion_table$hgnc_symbol <- toupper(refseq_embl_conversion_table$hgnc_symbol)
    write.table(refseq_embl_conversion_table, file = "refseq_embl_conversion_table.tsv", quote = F, sep = "\t", row.names = F)
} else {
    refseq_embl_conversion_table <- read.table("refseq_embl_conversion_table.tsv", header = T, colClasses = c("character", "character", "character"), fill = T)
}

transcript_without_gene_name_count <- sum(refseq_embl_conversion_table$hgnc_symbol == "")
## some transcript are originating from genes which do not have a gene name, only an Ensembl id. This is not a problem, I pointed out this just to get some insights of the data
print(paste("there are", transcript_without_gene_name_count, "transcript originating from a gene with only the Ensembl id and not a gene name"))

number_of_dinstinct_duplicate_genes <- sum((table(refseq_embl_conversion_table$hgnc_symbol) - 1) > 0)
number_of_non_duplicate_genes <- sum((table(refseq_embl_conversion_table$hgnc_symbol) - 1) == 0)
print(paste(number_of_non_duplicate_genes, "genes have only one transcript"))
print(paste(number_of_dinstinct_duplicate_genes, "genes have more than one transcript"))

## integrity check
matched_transcripts <- match(a$refseq_id, refseq_embl_conversion_table$refseq_mrna)
not_found_transcripts <- sum(is.na(matched_transcripts))
found_transcripts <- sum(!is.na(matched_transcripts))
print(paste(not_found_transcripts, "transcript not found in the refseq/ensembl conversion table"))
print(paste(found_transcripts, "found"))
predicted_gene_symbols <- refseq_embl_conversion_table$hgnc_symbol[matched_transcripts]
na_gene_symbols <- sum(is.na(predicted_gene_symbols))
if(na_gene_symbols != not_found_transcripts) {
    print(paste("warning:", na_gene_symbols, "gene symbols are NA,", not_found_transcripts, "expected"))
}
rows_differing_by_gene_symbol <- a$gene_symbol != predicted_gene_symbols
equal_rows_count <- sum(rows_differing_by_gene_symbol == F, na.rm = T)
differing_rows_count <- sum(rows_differing_by_gene_symbol == T, na.rm = T)
print(paste("the gene symbol coincides for", equal_rows_count, "rows"))
print(paste("the gene symbol differs for", differing_rows_count, "rows and is NA for", na_gene_symbols, "rows"))

print_the_differing_rows <- F
if(print_the_differing_rows) {
    rows_na_or_differing_by_gene_symbol <- rows_differing_by_gene_symbol | is.na(rows_differing_by_gene_symbol)
    to_print <- data.frame(refseq_id = a$refseq_id[rows_na_or_differing_by_gene_symbol],
                           targetscan_gene_symbol = a$gene_symbol[rows_na_or_differing_by_gene_symbol],
                           predicted_gene_symbols = predicted_gene_symbols[rows_na_or_differing_by_gene_symbol])
    print(to_print)
}
print("the reason behind so many transcript annotated to different genes is probably that I am using a wrong conversion table version, I will fix this. For the moment I will continue working with a subset of the data")
## for fixing the versions I will have to look at the output of these commands and at the url shown below
## listDatasets(useMart("ensembl"))
## listFilters(ensembl)
## listAttributes(ensembl)
## various databases maintained in Bioconductor: org.Hs.eg.db

## the matrix with the valid data is called b, I will only consider genes having only one transcript
b <- a
b$ensembl_gene_id <- predicted_gene_symbols
STOPPPED WORKING HEREEEEEEEEEEEE... I NEED TO SELECT THE ROWS WITH UNIQUE GENE SYMBOL
b$transcripts_with_the_same_gene <- 
b <- b[!is.na(b$ensembl_gene_id), ]
print(dim(b))
