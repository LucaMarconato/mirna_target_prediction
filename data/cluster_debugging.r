d <- read.table("../cluster_debugging", colClasses = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")); colnames(d) <- c("cluster_address", "sites", "front_mirna_id", "front_gene_id", "front_utr_start", "previous_value", "value", "original", "increment")
print(paste("max value =", max(d["previous_value"] + d["increment"])))
extracted <- d[d["previous_value"] + d["increment"] == max(d["previous_value"] + d["increment"]), "cluster_address"]
print(paste("rows with max value =", length(extracted)))
v <- extracted[1]
print(d[d$cluster_address == v,])

a <- read.table("processed/scored_interactions_processed.tsv", colClasses = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"), header = T)
