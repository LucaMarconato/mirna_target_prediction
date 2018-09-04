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
