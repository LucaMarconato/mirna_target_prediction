library("rjson")
library("hashmap")

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

print("miRNA 845 escluded")
new_maximized_device()
layout(rbind(matrix(rep(1,6),1,6),matrix(2:37, 6, 6, T)), heights = c(0.5,c(rep(3,10))))
old_par <- par(mar=c(0,0,0,0), font = 2)
plot.new()
text(0.5, 0.5, "distributions of the number of sites for the various genes, one miRNA per plot")
par(old_par)
for(i in o) {
    mirna_id <- rownames(m)[i]
    if(mirna_id != "845") {
        ## hist(m[i,][m[i,] > 0], main = paste("histogram of", mirna_id))
        barplot(table(m[i,][m[i,] > 0]), main = paste("miRNA", mirna_id))
    }
}

## plot expression profiles
new_maximized_device()
layout(matrix(c(1,2,5,6,3,4,7,8), 4, 2, T))
plot_files_filtered <- list("normal_mirna_expression_profile.tsv", "normal_gene_expression_profile.tsv", "tumor_mirna_expression_profile.tsv", "tumor_gene_expression_profile.tsv")
plot_files_filtered <- lapply(plot_files_filtered, function(x) paste(patient_folder, x, sep = ""))
plot_titles_temp <- list("log(RPM) for miRNAs in the healthy cells", "log(RPM) for genes in the healthy cells", "log(RPM) for miRNAs in the tumor cells", "log(RPM) for genes in the healthy cells")
plot_titles_filtered <- lapply(plot_titles_temp, function(x) paste(x, "(filtered)"))
for(i in 1:4) {
    ## the file always exists, but it may contain only the header
    expression_profile <- read.table(plot_files_filtered[[i]], colClasses = c("numeric", "numeric"), header = T)
    expression_profile <- expression_profile[order(expression_profile[[1]]),]
    expression_profile <- expression_profile[[2]]
    if(length(expression_profile) == 0) {
        plot.new()
    } else {
        ## plot(density(log10(expression_profile)), main = "")
        ## hist(log10(expression_profile), main = "")
        
        my_min <- -2
        my_max <- 2
        hist_data <- log10(expression_profile)
        hist_data <- hist_data[hist_data >= my_min & hist_data <= my_max]
        epsilon <- 0.001
        my_breaks <- seq(max(-2, min(hist_data)) - epsilon, min(2, max(hist_data)) + epsilon, length.out = 50)
        hist(hist_data, main = "", breaks = my_breaks)        
        ## plot(density(expression_profile), main = "")
        ## plot(log10(expression_profile), main = "")
    }
    title(plot_titles_filtered[[i]])
}

mirna_id_dictionary <- read.table("processed/mirna_id_dictionary.tsv", colClasses = c("character", "numeric"), header = T)
mirna_id_dictionary <- hashmap(mirna_id_dictionary[[1]], mirna_id_dictionary[[2]])
gene_id_dictionary <- read.table("processed/gene_id_dictionary.tsv", colClasses = c("character", "numeric"), header = T)
gene_id_dictionary <- hashmap(gene_id_dictionary[[1]], gene_id_dictionary[[2]])

## the graphs we are going to plot may be different from the one plotted just above for two reasons:
## first, before we plotted only mirnas and genes above the selected threshold
## second, before we plotted only mirnas and genes that were found in the interaction graph, here we plot the mirnas and the genes that belongs to the intersection between the GDC list and the TargetScan list of mirnas and gens, so we plot even mirnas and genes that are not found in the interaction graph

rpm_from_gdc_mirna_data <- function(filename_unfiltered_data) {
    a <- read.table(filename_unfiltered_data, colClasses = c("character", "numeric", "numeric", "factor"), header = T)
    a <- cbind(a, mirna_id_dictionary[[a$miRNA_ID]])
    colnames(a)[[5]] <- "mirna_id_cpp"
    a <- a[order(a$mirna_id_cpp),]    
    selected_rows <- a[!is.na(a$mirna_id_cpp),]
    rpm <- selected_rows["reads_per_million_miRNA_mapped"][[1]]
    rpm <- rpm[rpm > 0]
    return(rpm)
}

rpm_from_gdc_gene_data <- function(filename_unfiltered_data) {
    a <- read.table(filename_unfiltered_data, colClasses = c("character", "numeric"))    
    a <- a[!(a[[1]] %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")),]
    colnames(a)[[2]] <- "reads"
    total_reads <- sum(a$reads)
    a$reads <- a$reads/total_reads*1000000
    colnames(a)[[2]] <- "rpm"
    colnames(a)[[1]] <- "gene_id_and_version"
    a$gene_id_and_version <- unlist(lapply(a$gene_id_and_version, function(x) strsplit(x, "\\.")[[1]][[1]]))
    colnames(a)[[1]] <- "gene_id"
    a <- cbind(a, gene_id_dictionary[[a$gene_id]])
    colnames(a)[[3]] <- "gene_id_cpp"
    a <- a[!is.na(a$gene_id_cpp),]
    a <- a[order(a$gene_id_cpp),]
    rpm <- a$rpm
    rpm <- rpm[rpm > 0]
    return(rpm)
}

patient_info_file <- paste(patient_folder, "info.json", sep = "")
info_json <- fromJSON(file = patient_info_file)

## normal mirna
if(info_json$normal_mirnas$uuid != "NA") {
    normal_mirna_file_unfiltered <- paste(patient_folder, info_json$normal_mirnas$uuid, "/", info_json$normal_mirnas$file, sep = "")
} else {
    normal_mirna_file_unfiltered <- ""
}
## normal gene
if(info_json$normal_genes$uuid != "NA") {
    normal_gene_file_unfiltered <- paste(patient_folder, info_json$normal_genes$uuid, "/", info_json$normal_genes$file, sep = "")
} else {
    normal_gene_file_unfiltered <- ""
}
## tumor mirna
if(info_json$tumor_mirnas$uuid != "NA") {
    tumor_mirna_file_unfiltered <- paste(patient_folder, info_json$tumor_mirnas$uuid, "/", info_json$tumor_mirnas$file, sep = "")
} else {
    tumor_mirna_file_unfiltered <- ""
}
## tumor gene
if(info_json$tumor_genes$uuid != "NA") {
    tumor_gene_file_unfiltered <- paste(patient_folder, info_json$tumor_genes$uuid, "/", info_json$tumor_genes$file, sep = "")
} else {
    tumor_gene_file_unfiltered <- ""
}

plot_files_unfiltered <- list(normal_mirna_file_unfiltered, normal_gene_file_unfiltered, tumor_mirna_file_unfiltered, tumor_gene_file_unfiltered)
plot_titles_unfiltered <- lapply(plot_titles_temp, function(x) paste(x, "(unfiltered)"))

threeshold_rpm = 1

for(i in 1:4) {
    if(plot_files_unfiltered[[i]] != "") {
        if(i %in% c(1,3)) {
            rpm <- rpm_from_gdc_mirna_data(plot_files_unfiltered[[i]])
        } else {
            rpm <- rpm_from_gdc_gene_data(plot_files_unfiltered[[i]])
            ## browser()
        }
    }
    if(plot_files_unfiltered[[i]] == "" || length(rpm) == 0) {
        plot.new()
    } else {
        ## plot(density(log10(rpm)), main = "")
        ## hist(log10(rpm), main = "")
        
        my_min <- -2
        my_max <- 2
        hist_data <- log10(rpm)
        hist_data <- hist_data[hist_data >= my_min & hist_data <= my_max]
        epsilon <- 0.001
        my_breaks <- seq(max(-2, min(hist_data)) - epsilon, min(2, max(hist_data)) + epsilon, length.out = 50)
        hist(hist_data, main = "", breaks = my_breaks)
        abline(v = log10(threeshold_rpm))
        
        ## plot(density(rpm), main = "")
        ## abline(v = threeshold_rpm)
        ## plot(log10(rpm), main = "")
        ## abline(h = log10(threeshold_rpm))
    }
    title(plot_titles_unfiltered[[i]])
}
