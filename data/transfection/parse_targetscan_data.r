a <- read.table("targetscan_data/elife-05005-supp1-v1.csv", sep = ";", stringsAsFactors = F)
a <- a[3:nrow(a), ]
column_names <- as.character(a[1, ])
a <- a[2:nrow(a), ]
colnames(a) <- column_names
colnames(a)[[1]] <- "refseq_id"
colnames(a)[[2]] <- "used_in_training"
colnames(a)[[3]] <- "gene_symbol"
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

describe_a <- function(a)
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

describe_a(a)
