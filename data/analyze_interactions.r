library("rjson")

get_screen_resolution <- function() {
    return(c(2560, 1440))
    ## or automatically
    ## xdpyinfo_path <- "/opt/X11/bin/xdpyinfo"    
    ## cmd <- sprintf("%s | grep dimensions | perl -pe 's/^.*?([0-9]+x[0-9]+).*/$1/g' | tr 'x' ' '", xdpyinfo_path)
    ## output <- system(cmd, intern = T, ignore.stderr = T)
    ## return(as.numeric(unlist(strsplit(output, split = " "))))
}

## returns the screen size in inches
get_screen_physical_size <- function() {    
    return(c(24 + 1/4 + 1/8, 13 + 1/32))
    ## or automatically, but I preferred to do it manually via a "binary search"
    ## cmd <- sprintf("%s | grep dimensions | perl -pe 's/^.*?(\\([0-9]+x[0-9]+)\ millimeters.*/$1/g' | tr -d '(' | tr 'x' ' '", xdpyinfo_path)
    ## output <- system(cmd, intern = T, ignore.stderr = T)
    ## v <- as.numeric(unlist(strsplit(output, split = " ")))
    ## return(v/25.4)
}

new_maximized_device <- function() {
    resolution <- get_screen_resolution()
    physical_size <- get_screen_physical_size()
    ## for some reasons the ratio of the result I get is different from the ratio of the resolution so I am using only one value
    ## physical_size[2] <- physical_size[1] * resolution[2]/resolution[1]
    ## the option unit or units = "px" does not work for me, so I am using inches
    ## print(physical_size)
    dev.new(width = physical_size[1], height = physical_size[2])
}

close_all_devices <- function() {
    if(length(dev.list()) > 0) {
        for(l in dev.list()) {
            dev.off(l)   
        }
    }
}

patient_id <- "TCGA-CJ-4642"
patient_folder <- paste("patients/", patient_id, "/", sep = "")
matrix_filename <- paste(patient_folder, "mirna_gene_adjacency_matrix.mat", sep = "")
m <- read.table(matrix_filename, sep = "\t")
m <- as.matrix(m)
close_all_devices()
new_maximized_device()
layout(matrix(c(1,2,3,4,5,5,6,6,7,8),5,2,T),c(1,1))

title <- "distribution of the number of sites when varying (miRNA, gene) pairs"
barplot(table(m), main = title)
grid()

barplot(table(m[m >= 1]), main = title)
grid()

barplot(table(m[m >= 5]), main = title)
grid()

barplot(table(m[m >= 10]), main = title)
grid()

sites_per_mirna <- sapply(1:nrow(m), function(i) sum(m[i,]))
sites_per_mirna <- t(sites_per_mirna)
colnames(sites_per_mirna) <- rownames(m)
o <- order(sites_per_mirna)
barplot(sites_per_mirna[, o], main = "number of sites per miRNA")
grid()

genes_per_mirna <- apply(m, 1, function(x) sum(x != 0))
genes_per_mirna <- t(genes_per_mirna)
colnames(genes_per_mirna) <- rownames(m)
barplot(genes_per_mirna[, o], main = "number of genes per miRNA")
grid()

sites_per_gene <- sapply(1:ncol(m), function(i) sum(m[, i]))
hist(sites_per_gene, main = "number of sites per gene")
## grid()

mirnas_per_gene <- apply(m, 2, function(x) sum(x != 0))
hist(mirnas_per_gene, main = "number of miRNAs per gene")
## grid()

print("miRNA 846 escluded")
new_maximized_device()
layout(rbind(matrix(rep(1,6),1,6),matrix(2:37, 6, 6, T)), heights = c(0.5,c(rep(3,10))))
old_par <- par(mar=c(0,0,0,0), font = 2)
plot.new()
text(0.5, 0.5, "distributions of the number of sites for the various genes, one miRNA per plot")
par(old_par)
for(i in o) {
    mirna_id <- rownames(m)[i]
    if(mirna_id != "846") {
        ## hist(m[i,][m[i,] > 0], main = paste("histogram of", mirna_id))
        barplot(table(m[i,][m[i,] > 0]), main = paste("miRNA", mirna_id))
    }
}

mirna_id_dictionary <- read.table("processed/mirnas_with_scored_interactions.tsv", header = T)
mirna_id_dictionary <- cbind(mirna_id_dictionary, 1:length(mirna_id_dictionary[[1]]))
colnames(mirna_id_dictionary)[2] <- "mirna_id"

rpm_from_gdc_mirna_data <- function(filename, mirna_ids) {
    a <- read.table(filename, colClasses = c("character", "numeric", "numeric", "factor"), header = T)
    rpm <- a["reads_per_million_miRNA_mapped"]
    ## rpm_filtered <- rpm[rpm[,1] > 1,1]
    ## t <- table(rpm_filtered)
    ## barplot(t)
}

patient_info_file <- paste(patient_folder, "info.json", sep = "")
info_json <- fromJSON(file = patient_info_file)
new_maximized_device()
## normal mirna
if(info_json$normal_mirnas$uuid != "NA") {
    normal_mirna_file_unfiltered <- paste(patient_folder, info_json$normal_mirnas$uuid, "/", info_json$normal_mirnas$file)
} else {
    normal_mirna_file_unfiltered <- ""
}

## normal gene
if(info_json$normal_genes$uuid != "NA") {
    normal_gene_file_unfiltered <- paste(patient_folder, info_json$normal_genes$uuid, "/", info_json$normal_genes$file)
} else {
    normal_gene_file_unfiltered <- ""
}

## tumor mirna
if(info_json$tumor_mirnas$uuid != "NA") {
    tumor_mirna_file_unfiltered <- paste(patient_folder, info_json$tumor_mirnas$uuid, "/", info_json$tumor_mirnas$file)
} else {
    tumor_mirna_file_unfiltered <- ""
}

## tumor gene
if(info_json$tumor_genes$uuid != "NA") {
    tumor_gene_file_unfiltered <- paste(patient_folder, info_json$tumor_genes$uuid, "/", info_json$tumor_genes$file)
} else {
    tumor_gene_file_unfiltered <- ""
}

layout(matrix(c(1,2,5,6,3,4,7,8), 4, 2, T))
plot_files_filtered <- list("normal_mirna_expression_profile.tsv", "normal_gene_expression_profile.tsv", "tumor_mirna_expression_profile.tsv", "tumor_gene_expression_profile.tsv")
plot_files_filtered <- lapply(plot_files_filtered, function(x) paste(patient_folder, x, sep = ""))
plot_titles_temp <- list("log(RPM) for miRNAs in the healthy cells", "log(RPM) for genes in the healthy cells", "log(RPM) for miRNAs in the tumor cells", "log(RPM) for genes in the healthy cells")
plot_titles_filtered <- lapply(plot_titles_temp, function(x) paste(x, "(filtered)"))
plot_files_unfiltered <- list(normal_mirna_file_unfiltered, normal_gene_file_unfiltered, tumor_mirna_file_unfiltered, tumor_gene_file_unfiltered)
plot_titles_unfiltered <- lapply(plot_titles_temp, function(x) paste(x, "(unfiltered)"))

plot_files <- c(plot_files_filtered, plot_files_unfiltered)
plot_titles <- c(plot_titles_filtered, plot_files_unfiltered)

plot_data <- list()
for(i in 1:8) {
    if(i <= 4) {
        ## the file always exists, but it may contain only the header
        expression_profile <- read.table(plot_files[[i]], colClasses = c("numeric", "numeric", header = T))[[2]]
        plot_data <- c(plot_data, expression_profile)
    } else {
        ## the file could not exist
        if(plot_files[[i]] == "") {
            expression_profile = numeric()
        } else {
            if(i in c(5,7)) {
                expression_profile <- rpm_from_gdc_mirna_data(plot_files[[i]])
            } else {
                expression_profile <- rpm_from_gdc_gene_data(plot_files[[i]])
            }
        }
        plot_data <- c(plot_data, expression_profile)
    }
}

for(i in 1:8) {
    expression_profile <- plot_data[[i]]
    if(length(expression_profile) == 0) {
        plot.new()
    } else {
        plot(density(log10(expression_profile)), main = plot_titles[[i]])
    }   
}
