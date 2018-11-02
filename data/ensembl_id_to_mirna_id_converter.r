a <- read.table("raw/ensembl_id_to_mirbase_id.tsv", header = T, colClasses = c("character", "character", "character"), sep = "\t")
a <- a[a[[3]] != "", ]
colnames(a) <- c("ensembl_id", "transcript_id", "mirbase_id")
write.table(a, file = "processed/ensembl_id_to_mirbase_id.tsv", sep = "\t", row.names = F, quote = F)

