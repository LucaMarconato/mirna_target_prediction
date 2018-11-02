process_file <- function(file)
{
    a <- read.table(paste("raw/", file, sep = ""), colClasses = c("character", "numeric", "numeric", "numeric"))
    colnames(a) <- c("gene_id", "unstranded_count", "first_strand_count", "second_strand_count")

    a <- subset(a, grepl("ENSG", gene_id))
    print(a[1:10,])
}

files <- c("ENCFF495ZXC.tsv", "ENCFF902KUU.tsv")

for(file in files) {
    process_file(file)
}
