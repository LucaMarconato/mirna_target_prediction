a <- read.table("raw/ensembl_id_to_mirbase_id.tsv", header = T, colClasses = c("character", "character", "character"), sep = "\t")
a <- a[a[[3]] != "", ]
write.table(a, file = "processed/ensembl_id_to_mirbase_id.tsv", sep = "\t", row.names = F, quote = F)

b <- read.table("processed/mirna_id_dictionary.tsv", header = T, colClasses = c("character", "numeric"))

mirnas1 <- a[[3]]
mirnas2 <- b[[1]]
symdiff <- function( x, y) { setdiff( union(x, y), intersect(x, y))}

