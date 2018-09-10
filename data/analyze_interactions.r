library("rjson")
library("hashmap")
## library("RColorBrewer")

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

plot_adjacency_matrix_insights <- function(patient_folder) {
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

    ## print("miRNA 845 escluded")
    ## new_maximized_device()
    ## layout(rbind(matrix(rep(1,6),1,6),matrix(2:37, 6, 6, T)), heights = c(0.5,c(rep(3,10))))
    ## old_par <- par(mar=c(0,0,0,0), font = 2)
    ## plot.new()
    ## text(0.5, 0.5, "distributions of the number of sites for the various genes, one miRNA per plot")
    ## par(old_par)
    ## for(i in o) {
    ##     mirna_id <- rownames(m)[i]
    ##     if(mirna_id != "845") {
    ##         ## hist(m[i,][m[i,] > 0], main = paste("histogram of", mirna_id))
    ##         barplot(table(m[i,][m[i,] > 0]), main = paste("miRNA", mirna_id))
    ##     }
    ## }    
}

rpm_from_gdc_mirna_data <- function(filename_unfiltered_data, mirna_id_dictionary)
{
    a <- read.table(filename_unfiltered_data, colClasses = c("character", "numeric"), header = T)
    total_reads <- sum(a$reads)
    a$reads <- a$reads/total_reads*1000000
    colnames(a)[[2]] <- "rpm"
    a <- cbind(a, mirna_id_dictionary[[a$mirna_id]])
    colnames(a)[[3]] <- "mirna_id_cpp"
    a <- a[order(a$mirna_id_cpp),]    
    selected_rows <- a[!is.na(a$mirna_id_cpp),]
    rpm <- selected_rows["rpm"][[1]]
    rpm <- rpm[rpm > 0]
    return(rpm)
}

rpm_from_gdc_gene_data <- function(filename_unfiltered_data, gene_id_dictionary)
{
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

plot_expression_profiles_insights <- function(patient_folder)
{
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
            
            hist(log10(expression_profile), main = "")

            ## plot(density(expression_profile), main = "")
            
            ## plot(log10(expression_profile), main = "")
            
            ## my_min <- -2
            ## my_max <- 2
            ## hist_data <- log10(expression_profile)
            ## hist_data <- hist_data[hist_data >= my_min & hist_data <= my_max]
            ## epsilon <- 0.001
            ## my_breaks <- seq(max(-2, min(hist_data)) - epsilon, min(2, max(hist_data)) + epsilon, length.out = 50)
            ## hist(hist_data, main = "", breaks = my_breaks)      
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
    
    patient_info_file <- paste(patient_folder, "info.json", sep = "")
    info_json <- fromJSON(file = patient_info_file)

    ## normal mirna
    filename <- paste(patient_folder, "mirna_expression_profile_normal.tsv", sep = "")
    if(file.exists(filename)) {
        normal_mirna_file_unfiltered <- filename
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
    filename <- paste(patient_folder, "mirna_expression_profile_tumor.tsv", sep = "")
    if(file.exists(filename)) {
        tumor_mirna_file_unfiltered <- filename
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

    global_parameters_file <- "../global_parameters.json"
    global_parameters_json <- fromJSON(file = global_parameters_file)
    mirna_threshold_rpm = global_parameters_json$mirna_threshold_rpm
    gene_threshold_rpm = global_parameters_json$gene_threshold_rpm
    print(paste("mirna_threshold_rpm =", mirna_threshold_rpm))
    print(paste("gene_threshold_rpm =", gene_threshold_rpm))

    for(i in 1:4) {
        if(plot_files_unfiltered[[i]] != "") {
            if(i %in% c(1,3)) {
                rpm <- rpm_from_gdc_mirna_data(plot_files_unfiltered[[i]], mirna_id_dictionary)
            } else {
                rpm <- rpm_from_gdc_gene_data(plot_files_unfiltered[[i]], gene_id_dictionary)
                ## browser()
            }
        }
        if(plot_files_unfiltered[[i]] == "" || length(rpm) == 0) {
            plot.new()
        } else {
            ## plot(density(log10(rpm)), main = "")
            
            hist(log10(rpm), main = "")
            
            ## my_min <- -2
            ## my_max <- 2
            ## hist_data <- log10(rpm)
            ## hist_data <- hist_data[hist_data >= my_min & hist_data <= my_max]
            ## epsilon <- 0.001
            ## my_breaks <- seq(max(-2, min(hist_data)) - epsilon, min(2, max(hist_data)) + epsilon, length.out = 50)
            ## hist(hist_data, main = "", breaks = my_breaks)
            
            if(i %in% c(1,3)) {
                abline(v = log10(mirna_threshold_rpm), col = "red")
            } else {
                abline(v = log10(gene_threshold_rpm), col = "red")
            }
            
            ## plot(density(rpm), main = "")
            ## abline(v = threeshold_rpm)
            ## plot(log10(rpm), main = "")
            ## abline(h = log10(threeshold_rpm))
        }
        title(plot_titles_unfiltered[[i]])
    }
}

plot_overlapping_sites_insights <- function(patient_folder)
{
    ## I should use data.table instead of the code that I wrote
    filename <- paste(patient_folder, "overlapping_sites.tsv", sep = "")
    a <- read.table(filename, header = T)
    aggregate_values <- unlist(lapply(unique(a[[1]]), function(x) max(a[a[[1]] == x,2])))
    ## aggregate_values <- unlist(lapply(unique(a[[1]]), function(x) length(a[a[[1]] == x,2])))
    gene_id_to_aggregate_value <- hashmap(unique(a[[1]]), aggregate_values)
    aggregate_values_columns <- gene_id_to_aggregate_value[[a[[1]]]]
    a <- cbind(a, aggregate_values_columns)
    a <- a[order(a[,4]),]

    new_maximized_device()
    rows <- 12
    layout(matrix(1:rows,rows,1))
    old_par <- par(mar = c(0,4,0,0))
    u <- unique(a[,1])
    u_split <- split(u, ceiling(seq_along(u)/(length(u)/rows)))
    for(i in seq_len(rows)){
        a_split <- a[a[[1]] %in% u_split[[i]],]
        gene_id_to_seq_along <- hashmap(u_split[[i]], seq_along(u_split[[i]]))
        x_data <- gene_id_to_seq_along[[a_split[[1]]]]
        y_data <- a_split[[2]]
        ## barplot(u_split[[i]])
        plot(x_data, y_data, pch = ".", cex = 1, xaxt = 'n', col = rgb(0.5,0.5,0.5,0.25)) #, ylim = c(0, max(a[[2]]))
    }

    aggregate_values <- unlist(lapply(unique(a[[1]]), function(x) length(unique(a[a[[1]] == x,3]))))
    gene_id_to_aggregate_value <- hashmap(unique(a[[1]]), aggregate_values)
    aggregate_values_columns <- gene_id_to_aggregate_value[[a[[1]]]]
    a <- cbind(a, aggregate_values_columns)
    a <- a[order(a[,5]),]
    ## browser()
    new_maximized_device()
    
    rows <- 10
    first_rows_to_show <- 10
    last_rows_to_show <- 0
    total_rows_to_show <- first_rows_to_show + last_rows_to_show
    layout_vector <- rep(0, 2 * total_rows_to_show - 1)
    layout_vector[seq(1,2 * total_rows_to_show - 1, 2)] <- 1:total_rows_to_show
    layout_vector[seq(2,2 * total_rows_to_show - 1, 2)] <- seq(total_rows_to_show + 1, 2 * total_rows_to_show - 1, 1)
    widths_vector <- rep(1, 2 * total_rows_to_show - 1)
    widths_vector[seq(2,2 * total_rows_to_show - 1, 2)] <- lcm(1)
    layout(matrix(layout_vector, 2 * total_rows_to_show - 1, 1), widths = widths_vector)
    par(mar = c(0,4,0,0))
    u <- unique(a[,1])
    u_split <- split(u, ceiling(seq_along(u)/(length(u)/rows)))
    rows_drawn <- 0
    for(i in seq_len(rows)){
        if(i <= first_rows_to_show || i > rows - last_rows_to_show) {
            a_split <- a[a[[1]] %in% u_split[[i]],]
            gene_id_to_seq_along <- hashmap(u_split[[i]], seq_along(u_split[[i]]))
            x_data <- gene_id_to_seq_along[[a_split[[1]]]]
            y_data <- a_split[[3]]
            unique_y <- unique(y_data)
            ## barplot(u_split[[i]])
            ## use a light color when you are visualizing thousands of genes at the same time, a dark one otherwise
            ## my_color <- rgb(0.5,0.5,0.5,0.25) 
            my_color <- "black"
            plot(x_data, y_data, pch = ".", cex = 1, xaxt = 'n', col = my_color, axes = F) #, ylim = c(0, max(a[[2]]))
            light_gray <- rgb(0.9,0.9,0.9,0.3)
            ## abline(h = unique_y, lty = 1, col = light_gray)
            axis(side = 2, at = unique_y)
            rows_drawn <- rows_drawn + 1
        }
    }
    for(i in seq_len(rows_drawn - 1)) {
        plot(1, 1, type = 'n', axes = F, ann = F)
        abline(h = 1, col = "gray")
    }

    new_maximized_device()
    layout(matrix(1:2,2,1))
    barplot(table(a[[3]]))
    grid()
    barplot(table(a[[3]][a[[3]] > 1]))
    grid()
}

study_overlaps_of_specific_gene <- function(patient_folder, gene_id = -1)
{
    filename <- paste(patient_folder, "overlapping_sites.tsv", sep = "")
    a <<- read.table(filename, header = T)
    if(gene_id == -1) {
        gene_id <- a[sample(nrow(a), 1), ][[1]]
    }
    my_palette <- palette()
    my_palette[7] <- "orange"
    new_maximized_device()
    rows_first_plot <- 60
    layout(matrix(seq(1,rows_first_plot + 2), rows_first_plot + 2, 1), heights = c(lcm(2), rep(1, rows_first_plot), rows_first_plot/4))
    old_par <- par(mar = c(0,0,2,0))
    plot.new()
    title(paste("sites for gene_id =", gene_id))
    legend_labels = c(paste(0:6, "overlaps"),">= 7 overlaps")
    legend("top", legend_labels, col = my_palette[8:1], ncol = 8, pch = 15)
    mtext("position in the 3' UTR", side=3, at=par("usr")[1]+0.01*diff(par("usr")[1:2]), adj = 0, cex = 0.5)
    
    rows <- a[a[[1]] == gene_id, ]
    max_site_position <- max(rows[[2]])
    intervals <- seq(0, max_site_position + 1, length.out = rows_first_plot + 1)
    left_values <- intervals[1:(length(intervals) - 1)]
    ## this -1 compensates for the +1 of two rows above
    right_values <- intervals[2:length(intervals)] - 1
    utr_length_so_far <- 0
    ## print(intervals)
    ## print(left_values)
    ## print(right_values)
    for(j in seq_along(left_values)) {
        left_value <- left_values[j]
        right_value <- right_values[j]        
        if(j == 1) {
            par(mar = c(0,0,0,0))
        } else {
            par(mar = c(0,0,0,0))
        }
        rows_split <- rows[rows[[2]] >= left_value & rows[[2]] <= right_value,]
        for(i in 1:8) {
            if(i < 8) {
                sites <- rows_split[rows_split[[3]] == i - 1, ]
            } else {
                sites <- rows_split[rows_split[[3]] >= 7, ]
            }
            my_xlim <- c(left_value, right_value)
            my_color <- my_palette[8 - i + 1]
            my_pch <- 1
            my_cex <- 1
            if(length(sites[[2]])) {
                if(i == 1) {
                    plot(sites[[2]], rep(1, length(sites[[2]])), xlim = my_xlim, xlab = "position in 3' UTR", yaxt = 'n', ylab = '', col = my_color, pch = my_pch, cex = my_cex, axes = F)
                } else {
                    points(sites[[2]], rep(1, length(sites[[2]])), xlim = my_xlim, xlab = "position in 3' UTR", yaxt = 'n', ylab = '', col = my_color, pch = my_pch, cex = my_cex)
                }
            } else {
                if(i == 1) {
                    plot.new()
                }
            }
            box(col = "lightgray")
            ## TODO: BUG: the text is not placed in the last row
            mtext(paste(round(left_value)),side=3, at=par("usr")[1]+0.01*diff(par("usr")[1:2]), cex = 0.5)
        }
    }

    par(old_par)
    counts <- sapply(0:6, function(x) length(rows[[3]][rows[[3]] == x]))
    greater_than7 <- length(rows[[3]][rows[[3]] >= 7])
    if(greater_than7 == 0) {
        greater_than7 <- 0
    }
    counts <- c(counts, greater_than7)
    ## print(counts); print(greater_than7)
    ## h <- hist(rows[[3]], plot = F)
    ## colors_for_breaks <- 8 - h$breaks
    ## colors_for_breaks[colors_for_breaks < 0] <- 0
    colors_for_breaks <- 8:1
    ## my_names <- h$mids - 0.5
    ## found <- match(7, my_names)
    ## if(!is.na(found)) {
    ##     my_names[found] <- ">= 7"
    ## }
    my_names <- paste(0:6)
    my_names <- c(my_names, ">= 7")
    barplot(counts, names.arg = my_names, col = my_palette[colors_for_breaks], main = "Distribution of overlapping sites", horiz = T, las = 1, axes = F)
    my_ticks <- round(seq(0,max(counts),5))
    axis(side = 1, at = my_ticks)
    abline(v = my_ticks, lty = 2, col = "lightgray")
}

patient_id <- "TCGA-CJ-4642"
patient_folder <- paste("patients/", patient_id, "/", sep = "")
## plot_adjacency_matrix_insights(patient_folder)
## plot_expression_profiles_insights(patient_folder)
## plot_overlapping_sites_insights(patient_folder)
study_overlaps_of_specific_gene(patient_folder)
study_overlaps_of_specific_gene(patient_folder, gene_id = 118)
## study_overlaps_of_specific_gene(patient_folder, gene_id = 545)
