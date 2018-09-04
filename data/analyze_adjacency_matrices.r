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
    print(physical_size)
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
layout(matrix(c(1,2,3,4,5,5,5,5),4,2, byrow = T))
barplot(table(m))
grid()
barplot(table(m[m >= 1]))
grid()
barplot(table(m[m >= 5]))
grid()
barplot(table(m[m >= 10]))
grid()
genes_per_mirna <- sapply(1:nrow(m), function(i) sum(m[i,][m[i,]>0]))
genes_per_mirna <- t(genes_per_mirna)
colnames(genes_per_mirna) <- rownames(m)
barplot(genes_per_mirna)
grid()

print("miRNA 846 escluded")
new_maximized_device()
layout(matrix(1:36, 6, 6, T))
for(i in 1:37) {
    mirna_id <- rownames(m)[i]
    if(mirna_id != "846") {
        hist(m[i,][m[i,] > 0], main = paste("Histogram of", mirna_id))
    }
}
