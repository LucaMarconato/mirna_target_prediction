library(purrr)
library(raster)
library(RColorBrewer)
library(plotrix)
library(gridExtra)
library(lattice)
library(latticeExtra)
rm(list = ls())
source("my_device.r")

analyze_convergence <- function(simulation_output_path)
{
    filename <- paste(simulation_output_path, "matchings_prediction.tsv", sep = "")
    a <<- read.table(filename, header = T, colClasses = c("numeric", "numeric", "numeric", "numeric", "numeric"))
    x_values <- seq(1, length(a[[1]]))

    new_maximized_device()
    layout(matrix(c(1,2,4,1,3,5,6,8,10,7,9,11), 4, 3, byrow = T))
    plot(x_values, a$cumulative_scaling, main = "cumulative scaling", pch = 20)
    grid()
    plot(x_values, a$avg_mirna_level, main = "average mirna level", pch = 20)
    grid()
    plot(x_values, log10(a$avg_mirna_level), main = "average mirna level (log10)", pch = 20)
    grid()
    plot(x_values, a$avg_cluster_level, main = "average cluster level", pch = 20)
    grid()
    plot(x_values, log10(a$avg_cluster_level), main = "average cluster level (log10)", pch = 20)
    grid()
    plot(x_values, a$mirna_total_exchange, main = "mirna total exchange in the step", pch = 20)
    grid()
    plot(x_values, log10(a$mirna_total_exchange), main = "mirna total exchange in the step (log10)", pch = 20)
    grid()
    plot(x_values, a$cluster_total_exchange, main = "cluster total exchange in the step", pch = 20)
    grid()
    plot(x_values, log10(a$cluster_total_exchange), main = "cluster total exchange in the step (log10)", pch = 20)
    grid()
    plot(x_values, cumsum(a$mirna_total_exchange), main = "mirna cumulative total exchange", pch = 20)
    grid()
    plot(x_values, cumsum(a$cluster_total_exchange), main = "cluster cumulative total exchange", pch = 20)
    grid()
}

analyze_mirna_expression_profiles <- function(simulation_output_path)
{
    print("analyzing mirna expression profiles")
    expression_profile_folder <- paste(simulation_output_path, "mirna_expression_profiles", sep = "")
    expression_profile_copy_folder <- paste(simulation_output_path, "mirna_expression_profiles_copy", sep = "")
    system(paste("mkdir -p ", expression_profile_copy_folder, sep = ""))
    system(paste("rm ", expression_profile_copy_folder, "/mirna_expression_profile_*", sep = ""))
    system(paste("cp -r ", expression_profile_folder, "/* ", expression_profile_copy_folder, sep = ""))
    files <- list.files(path = expression_profile_copy_folder, pattern = "*.tsv", full.names = T, recursive = F)

    ## mirna_dynamics is the list of all the y_data
    mirna_dynamics <- list()
    mirna_ids_ordered <- NULL
    original_y_data <- NULL
    original_log_difference <- NULL
    i <- 0
    for(file in files) {
        t <- read.table(file, header = T, colClasses = c("numeric", "numeric", "numeric"))
        if(is.null(mirna_ids_ordered)) {
            my_order <- order(t$reads)
            mirna_ids_ordered <- t$mirna_id[my_order]
        }
        new_path <- tools::file_path_sans_ext(file)
        new_path <- paste(new_path, ".png", sep = "")
        png(filename = new_path, width = 1920, height = 1080)
        layout(matrix(1:3, 3, 1))
        ## print("TODO: check if the timestep has been split correctly")
        timestep <- strsplit(file, "mirna_expression_profile_")[[1]][2]
        timestep <- strsplit(timestep, ".tsv")[[1]][1]
        my_order <- match(mirna_ids_ordered, t$mirna_id)
        x_data <- 1:length(t$reads)
        y_data <- t$reads[my_order]
        mirna_dynamics[[length(mirna_dynamics) + 1]] <- y_data
        if(is.null(original_y_data)) {
            original_y_data <- y_data
        }
        plot(x_data, y_data, main = paste(timestep, "miRNA expression, reads"), ylim = c(0, max(original_y_data)))
        points(x_data, original_y_data, ylim = c(0, max(original_y_data)), col = "red", pch = "-")

        plot(x_data, log10(y_data), main = paste(timestep, "miRNA expression, reads (log10)"), ylim = c(-20, 0))
        points(x_data, log10(original_y_data), ylim = c(-20, 0), col = "red", pch = "-")

        if(i > 0) {
            previous_file <- files[i]
            previous_t <- read.table(previous_file, header = T, colClasses = c("numeric", "numeric"))
            previous_my_order <- match(mirna_ids_ordered, previous_t$mirna_id)
            previous_y_data <- previous_t$reads[previous_my_order]
            difference <- previous_y_data - y_data
            log_difference <- log10(difference)
            ## log_difference[log_difference == -Inf] = 0
            if(sum(is.nan(log_difference)) > 0) {
                print("warning: found a nan, probably this is due by having reached the machine precision during the simulation")
                log_difference[is.nan(log_difference)] = 0
                ## browser()
            }
            if(is.null(original_log_difference)) {
                original_log_difference <- log_difference
            }
            plot(x_data, log_difference, main = paste(timestep, "exchanged"), ylim = c(-20, 0))
            points(x_data, original_log_difference, ylim = c(-20, 0), col = "red", pch = "-")
        }
        dev.off()
        i <- i + 1
    }
    filename <- paste(expression_profile_copy_folder, "/animation.gif", sep = "")
    print(paste("generating", filename))
    system(paste("convert -delay 10 -loop 0 ", expression_profile_copy_folder, "/*.png ", filename, sep = ""))
    system(paste("open ", expression_profile_copy_folder, "/animation.gif", sep = ""))
    return(mirna_dynamics)
}

analyze_cluster_expression_profiles <- function(simulation_output_path)
{
    print("analyzing cluster expression profiles")
    expression_profile_folder <- paste(simulation_output_path, "cluster_expression_profiles", sep = "")
    expression_profile_copy_folder <- paste(simulation_output_path, "cluster_expression_profiles_copy", sep = "")
    system(paste("mkdir -p ", expression_profile_copy_folder, sep = ""))
    system(paste("rm ", expression_profile_copy_folder, "/cluster_expression_profile_*", sep = ""))
    system(paste("cp -r ", expression_profile_folder, "/* ", expression_profile_copy_folder, sep = ""))
    files <- list.files(path = expression_profile_copy_folder, pattern = "*.tsv", full.names = T, recursive = F)

    ## cluster_dynamics is the list of all the y_data
    cluster_dynamics <- list()
    cluster_addresses_ordered <- NULL
    original_y_data <- NULL
    ## original_log_difference <- NULL
    file_count <- length(files)
    max_timesteps <- 100
    for(i in seq_len(file_count)) {
        ## if(i <= 2) {
        ##     print(paste("skipping timestep", i))
        ##     next
        ## }
        print(paste(round(i/min(file_count, max_timesteps)*100), "%", sep = ""))
        file <- files[i]
        t <- read.table(file, header = T, colClasses = c("numeric", "numeric"))
        if(is.null(cluster_addresses_ordered)) {
            my_order <- order(t$relative_expression)
            cluster_addresses_ordered <- t$cluster_address[my_order]
        }
        ## print("TODO: check if the timestep has been split correctly")
        timestep <- strsplit(file, "cluster_expression_profile_")[[1]][2]
        timestep <- strsplit(timestep, ".tsv")[[1]][1]
        my_order <- match(cluster_addresses_ordered, t$cluster_address)
        x_data <- 1:length(t$relative_expression)
        y_data <- t$relative_expression[my_order]
        cluster_dynamics[[length(cluster_dynamics) + 1]] <- y_data
        if(is.null(original_y_data)) {
            original_y_data <- y_data
        }
        new_path <- tools::file_path_sans_ext(file)
        new_path <- paste(new_path, ".png", sep = "")
        rows <- 10

        png(filename = new_path, width = 1920, height = 1080)
        old_par <- par(mar = c(0,4,2,0))

        if(length(x_data) < 50) {
            rows <- 1
        }
        u <- seq_along(x_data)
        u_split <- split(u, ceiling(seq_along(u)/(length(u)/rows)))
        rows <- min(length(u_split), rows)

        ## we add an extra row for the title
        heights_vector <- rep(1, rows + 1)
        heights_vector[1] <- lcm(2)
        layout(matrix(seq(1, rows + 1), rows + 1, 1), heights = heights_vector)
        plot.new()
        title(main = paste(timestep, "cluster expression"))
        par(mar = c(2,4,2,0))

        for(i in seq_len(rows)){
            x_data_split <- x_data[u_split[[i]]]
            y_data_split <- y_data[u_split[[i]]]
            original_y_data_split <- original_y_data[u_split[[i]]]
            y_lim_split <- c(0, max(original_y_data_split))

            plot(x_data_split, y_data_split, ylim = y_lim_split, pch = ".")
            points(x_data_split, original_y_data_split, ylim = y_lim_split, col = "red", pch = "-")
            ## plot(x_data, log10(y_data), main = paste(timestep, "miRNA expression (log10)"), ylim = c(-20, 0))
            ## points(x_data, log10(original_y_data), ylim = c(-20, 0), col = "red", pch = "-")
        }

        ## ## probably useless
        ## par(mar = old_par)
        dev.off()
    }
    filename <- paste(expression_profile_copy_folder, "/animation.gif", sep = "")
    print(paste("generating", filename))
    system(paste("convert -delay 10 -loop 0 ", expression_profile_copy_folder, "/*.png ", filename, sep = ""))
    print("generating .gif")
    system(paste("open ", expression_profile_copy_folder, "/animation.gif", sep = ""))
}

generate_readable_dynamics_log <- function(simulation_output_path)
{
    ## note that we are reading data from the cluster_expression_profiles folder and not the cluster_expression_profiles_copy folder
    ## the animated .gifs are made out of the copied folder, to avoid accidental deletion by the C++ code, so be sure that you have not called
    ## the C++ code after having generated the animations, otherwise this function will agree with the data just produced by the C++ code and maybe not
    ## with the animations, unless you have generated them again

    extract_dynamics <- function(element)
    {
        expression_profile_folder <- paste(simulation_output_path, element, "_expression_profiles", sep = "")
        files <- list.files(path = expression_profile_folder, pattern = "*.tsv", full.names = T, recursive = F)
        ## element_dynamics is the list of all the y_data
        element_dynamics <- list()
        element_addresses_ordered <- NULL
        file_count <- length(files)
        for(i in seq_len(file_count)) {
            file <- files[i]
            t <- read.table(file, header = T, colClasses = c("numeric", "numeric"))
            if(is.null(element_addresses_ordered)) {
                my_order <- order(t$relative_expression)
                element_addresses_ordered <- t[[1]][my_order]
            }
            ## print("TODO: check if the timestep has been split correctly")
            timestep <- strsplit(file, paste(element, "_expression_profile_", sep = ""))[[1]][2]
            timestep <- strsplit(timestep, ".tsv")[[1]][1]
            my_order <- match(element_addresses_ordered, t[[1]])
            y_data <- t$relative_expression[my_order]
            element_dynamics[[length(element_dynamics) + 1]] <- y_data
        }
        return(element_dynamics)
    }

    mirna_dynamics <- extract_dynamics("mirna")
    cluster_dynamics <- extract_dynamics("cluster")

    if(length(mirna_dynamics) != length(cluster_dynamics)) {
        print(paste("error: length(mirna_dynamics) = ", length(mirna_dynamics), ", length(cluster_dynamics) = ", length(cluster_dynamics), sep = ""))
        stop("fatal error")
    }
    maximum_mirnas <- 5
    maximum_clusters <- 5
    if(length(mirna_dynamics[[1]]) > maximum_mirnas) {
        print("too many mirnas, not exporting the dynamics in a readable file")
        return()
    } else if(length(cluster_dynamics[[1]]) > maximum_clusters) {
        print("too many clusters, not exporting the dynamics in a readable file")
        return()
    }
    filename <- paste(simulation_output_path, "dynamics.txt", sep = "")
    s <- ""
    i <- 1

    get_dynamics_string <- function(element_dynamics)
    {
        state <- element_dynamics[[i]]
        state_string_formatted <- sapply(state, function(x)
        {
            if(x < 0.01) sprintf("%.e", x) else round(x, 3)
        })
        state_string_collapsed <- paste(state_string_formatted, collapse = ", ")
    }
    timesteps <- length(mirna_dynamics)
    for(i in 1:timesteps) {
        s <- paste(s, get_dynamics_string(mirna_dynamics), sep = "")
        s <- paste(s, "\n", sep = "")
        s <- paste(s, get_dynamics_string(cluster_dynamics), sep = "")
        s <- paste(s, "\n", sep = "\n")
    }
    out <- file(filename)
    writeLines(s, out)
    close(out)
    print(paste("saved dynamics into", filename))
}

analyze_dynamics_for_small_interaction_graphs <- function(simulation_output_path)
{
    print("loading the interaction matrix")
    expression_profile_folder <- paste(simulation_output_path, "mirna_expression_profiles", sep = "")
    files <- list.files(path = expression_profile_folder, pattern = "*.tsv", full.names = T, recursive = F)
    t <- read.table(files[1], header = T, colClasses = c("numeric", "numeric"))
    my_order <- order(t$relative_expression)
    mirna_ids_ordered <<- t$mirna_id[my_order]

    expression_profile_folder <- paste(simulation_output_path, "cluster_expression_profiles", sep = "")
    files <- list.files(path = expression_profile_folder, pattern = "*.tsv", full.names = T, recursive = F)
    t <- read.table(files[1], header = T, colClasses = c("numeric", "numeric"))
    my_order <- order(t$relative_expression)
    cluster_addresses_ordered <<- t$cluster_address[my_order]

    ## x_data <- 1:length(t$relative_expression)
    ## y_data <- t$relative_expression[my_order]

    maximum_mirnas <- 5
    maximum_clusters <- 5
    if(length(mirna_ids_ordered) > maximum_mirnas) {
        print("too many mirnas, not showing the interaction graph")
        return()
    } else if(length(cluster_addresses_ordered) > maximum_clusters) {
        print("too many clusters, not showing the interaction graph")
        return()
    }

    matrix_filename <- paste(simulation_output_path, "interaction_matrix.mat", sep = "")
    m <- read.table(matrix_filename, sep = "\t")
    m <- as.matrix(m)
    my_order <- match(mirna_ids_ordered, rownames(m))
    m <- m[my_order, ]
    colnames(m) <- lapply(colnames(m), function(x) as.numeric(strsplit(x, "X")[[1]][2]))
    my_order <- match(cluster_addresses_ordered, colnames(m))
    m <<- m[, my_order]
    new_maximized_device()
    ## 16/9 is my monitor size ratio
    plot(raster(m), asp = 9/16, col = sapply(seq(1, 0, length.out = max(m) + 1), function(x) rgb(x,x,x)))
    ## my_order_rows <- // need to export the matrix
}

analyze_probabilities <- function(patient_folder, simulation_output_path)
{
    new_maximized_device()
    layout(matrix(1:4, 2, 2, byrow = T))
    filename <- paste(simulation_output_path, "r_ic_values.tsv", sep = "")
    a <- read.table(filename, header = T, colClasses = c("numeric", "numeric", "numeric"))
    hist(log10(a$r_ic))

    filename <- paste(simulation_output_path, "r_ijk_values.tsv", sep = "")
    a <- read.table(filename, header = T, colClasses = c("numeric", "numeric", "numeric"))
    hist(log10(a$r_ijk))

    filename <- paste(simulation_output_path, "p_c_bound_values.tsv", sep = "")
    a <- read.table(filename, header = T, colClasses = c("numeric", "numeric"))
    aa <<- a
    hist(log10(a$p_c_bound))

    filename <- paste(simulation_output_path, "p_j_downregulated_given_c_bound_values.tsv", sep = "")
    a <- read.table(filename, header = T, colClasses = c("numeric", "numeric", "numeric"))
    hist(a$p_j_downregulated_given_c_bound_values)

    new_maximized_device()
    par(mfrow = c(1, 2))
    filename <- paste(simulation_output_path, "predicted_downregulation/p_j_downregulated_values_999999.tsv", sep = "")
    a <- read.table(filename, header = T, colClasses = c("numeric", "numeric"))
    hist(a$p_j_downregulated_values)

    filename <- paste(patient_folder, "gene_expression_profile_exported.tsv", sep = "")
    b <- read.table(filename, header = T, colClasses = c("numeric", "numeric"))

    d <- merge(a, b)
    d$rpm_downregulated <- d$rpm * d$p_j_downregulated_values
    hist(log(d$rpm_downregulated))

    new_maximized_device()
    print("analyzing downregulation")
    predicted_downregulation_folder <- paste(simulation_output_path, "predicted_downregulation", sep = "")
    predicted_downregulation_copy_folder <- paste(simulation_output_path, "predicted_downregulation_copy", sep = "")
    system(paste("mkdir -p ", predicted_downregulation_copy_folder, sep = ""))
    system(paste("rm ", predicted_downregulation_copy_folder, "/p_j_downregulated_values_*", sep = ""))
    system(paste("cp -r ", predicted_downregulation_folder, "/* ", predicted_downregulation_copy_folder, sep = ""))
    files <- list.files(path = predicted_downregulation_copy_folder, pattern = "*.tsv", full.names = T, recursive = F)

    genes_ordered <- NULL
    file_count <- length(files)
    max_timesteps <- 100
    for(i in seq_len(file_count)) {
        ## if(i > 1) {
        ##     print("not generating the animation")
        ##     break
        ## }
        print(paste(round(i/min(file_count, max_timesteps)*100), "%", sep = ""))
        file <- files[i]
        t <- read.table(file, header = T, colClasses = c("numeric", "numeric"))
        t <- merge(t, b)
        t$rpm_downregulated <- t$rpm * t$p_j_downregulated_values
        t$final_rpm <- t$rpm - t$rpm_downregulated
        if(is.null(genes_ordered)) {
            my_order <- order(t$rpm)
            genes_ordered <- t$gene_id[my_order]
        }
        ## print("TODO: check if the timestep has been split correctly")
        timestep <- strsplit(file, "predicted_downregulation_")[[1]][2]
        timestep <- strsplit(timestep, ".tsv")[[1]][1]
        my_order <- match(genes_ordered, t$gene_id)
        x_data <- 1:length(t$final_rpm)
        y_data <- t$final_rpm[my_order]
        original_y_data <- t$rpm[my_order]

        new_path <- tools::file_path_sans_ext(file)
        new_path <- paste(new_path, ".png", sep = "")
        rows <- 10

        png(filename = new_path, width = 1920, height = 1080)
        old_par <- par(mar = c(0,4,2,0))

        if(length(x_data) < 50) {
            rows <- 1
        }
        u <- seq_along(x_data)
        u_split <- split(u, ceiling(seq_along(u)/(length(u)/rows)))
        rows <- min(length(u_split), rows)

        ## we add an extra row for the title
        heights_vector <- rep(1, rows + 1)
        heights_vector[1] <- lcm(2)
        layout(matrix(seq(1, rows + 1), rows + 1, 1), heights = heights_vector)
        plot.new()
        ## print("TODO: check if the timestep has been split correctly")
        title(main = paste(timestep, "predicted gene expression after downregulation"))
        par(mar = c(2,4,2,0))

        for(i in seq_len(rows)){
            x_data_split <- x_data[u_split[[i]]]
            y_data_split <- y_data[u_split[[i]]]
            original_y_data_split <- original_y_data[u_split[[i]]]
            y_lim_split <- c(0, max(original_y_data_split))

            plot(x_data_split, y_data_split, ylim = y_lim_split, pch = ".")
            points(x_data_split, original_y_data_split, ylim = y_lim_split, col = "red", pch = ".")
        }
        dev.off()
    }
    filename <- paste(predicted_downregulation_copy_folder, "/animation.gif", sep = "")
    print(paste("generating", filename))
    system(paste("convert -delay 10 -loop 0 ", predicted_downregulation_copy_folder, "/*.png ", filename, sep = ""))
    print("generating .gif")
    system(paste("open ", predicted_downregulation_copy_folder, "/animation.gif", sep = ""))

    new_maximized_device()
    filename <- paste(predicted_downregulation_copy_folder, "/p_j_downregulated_values_999999.tsv", sep = "")
    t <- read.table(filename, header = T, colClasses = c("numeric", "numeric"))
    t <- merge(t, b)
    t$rpm_downregulated <- t$rpm * t$p_j_downregulated_values
    t$final_rpm <- t$rpm - t$rpm_downregulated
    my_order <- match(genes_ordered, t$gene_id)
    x_data <- 1:length(t$final_rpm)
    y_data <- t$final_rpm[my_order]
    ## cd stands for "cumulative downregulation"
    cd <- cumsum(y_data)
    plot(cd, main = "cumulative downregulation, genes order by ascending initial expression, lines at 1/4, 1/2, 3/4")
    abline(v = which.min(abs(cd - tail(cd, n = 1)/2)))
    abline(v = which.min(abs(cd - tail(cd, n = 1)/4)))
    abline(v = which.min(abs(cd - tail(cd, n = 1)*3/4)))
}

## export_predictions_using_ensembl_ids <- function(predicted_downregulation_folder, time_index, initial_rpm, final_rpm)
## {
##     path <- paste(predicted_downregulation_folder, "/predicted_downregulation_", index, ".tsv", sep = "")
##     predicted_downregulation <- data.frame()
##     stop("stop here!")
## }

analyze_predictions_of_perturbed_data <- function(patient_folder, simulation_output_paths, gene_ids = NULL, consider_only_specified_gene_ids = F, consider_only_gene_ids_of_current_experiment = F) ## , generate_images = T
{
    print("analyzing predictions of perturbed data")
    if(consider_only_specified_gene_ids) {
        if(is.null(gene_ids)) {
            print(paste("error: is.null(gene_ids) = ", is.null(gene_ids), ", consider_only_specified_gene_ids = ", consider_only_specified_gene_ids, sep = ""))
            stop("abort")
        }
    }
    filename <- paste(patient_folder, "gene_expression_profile_exported.tsv", sep = "")
    b <- read.table(filename, header = T, colClasses = c("numeric", "numeric"))
    gene_id_dictionary <- read.table("processed/gene_id_dictionary.tsv", header = T, colClasses = c("character", "numeric"))
    ## browser()
    all_gene_ids_cpp <- gene_id_dictionary$gene_id_cpp
    expressed_gene_ids_cpp <- b$gene_id
    not_expressed_gene_ids_cpp <- setdiff(all_gene_ids_cpp, expressed_gene_ids_cpp)
    b <- rbind(b, data.frame(gene_id = not_expressed_gene_ids_cpp, rpm = rep(0, length(not_expressed_gene_ids_cpp))))
    b <- b[order(b$gene_id), ]

    ## here I should use the equivalent of a static variable in C, but I do not know how to do it in R so I am using this temporary workaround
    if(!exists("genes_ordered_static_variable")) {
        genes_ordered_static_variable <<- NULL
    }
    ## number_of_genes_per_row <<- NULL

    message_about_gene_ids_already_printed <<- F

    plot_dataset <- function(dataset_index) {
        predicted_downregulation_folder <- paste(simulation_output_paths[[dataset_index]], "predicted_downregulation", sep = "")
        file <- paste(predicted_downregulation_folder, "/p_j_downregulated_values_999999.tsv", sep = "")
        t <- read.table(file, header = T, colClasses = c("numeric", "numeric"))
        t <- merge(t, b, all.y = T)
        na_count <- sum(is.na(t$p_j_downregulated_values))
        t$p_j_downregulated_values[is.na(t$p_j_downregulated_values)] <- rep(0, na_count)
        ## browser()
        t$rpm_downregulated <- t$rpm * t$p_j_downregulated_values
        ## browser()
        ## note that this RPMs do not sum to 1000000 because I want the to be directly comparable with the RPM values used for making the predictions
        t$final_rpm <- t$rpm - t$rpm_downregulated
        ## export_predictions_using_ensembl_ids(predicted_downregulation_folder = predicted_downregulation_folder, time_index = "999999", initial_rpm = t$rpm, final_rpm = t$final_rpm)
        ## if(!generate_images) {
        ##     return()
        ## }
        if(is.null(genes_ordered_static_variable)) {
            my_order <- order(t$rpm)
            ## browser()
            genes_ordered_static_variable <<- t$gene_id[my_order]
        }
        my_order <- match(genes_ordered_static_variable, t$gene_id)
        x_data <- 1:length(t$final_rpm)
        y_data <- t$final_rpm[my_order]
        is_the_id_in_gene_ids <- (t$gene_id[my_order] %in% gene_ids[[dataset_index]])
        original_y_data <- t$rpm[my_order]
        if(message_about_gene_ids_already_printed == F) {
            message_about_gene_ids_already_printed <<- T
            print(paste(length(is_the_id_in_gene_ids),
                        "genes in the interaction graph, of which",
                        sum(is_the_id_in_gene_ids),
                        "are considered in the gene_ids variable"))
        }
        genes_skipped_in_distance_based_prediction_filename <- paste(simulation_output_paths[[dataset_index]], "/genes_skipped_in_distance_based_prediction.tsv", sep = "")
        if(file.exists(genes_skipped_in_distance_based_prediction_filename)) {
            skipped <- read.table(genes_skipped_in_distance_based_prediction_filename, header = T, colClasses = c("numeric"))[[1]]
            ## print(paste("genes skipped in distance based prediction:", paste(skipped, collapse = ", ")))
            genes_skipped_in_distance_based_prediction <- rep(F, length(x_data))
            genes_skipped_in_distance_based_prediction[match(skipped, genes_ordered_static_variable)] <- T
            ## print(match(skipped, genes_ordered_static_variable))
        } else {
            genes_skipped_in_distance_based_prediction <- rep(F, length(x_data))
        }

        rows <- 10
        ## if(i > 1) {
        ##     par(new = T)
        ## }
        if(length(x_data) < 50) {
            rows <- 1
        }
        u <- seq_along(x_data)
        ## if(is.null(number_of_genes_per_row)) {
        number_of_genes_per_row <- floor(length(u) / rows)
            ## print("setting number_of_genes_per_row")
        ## }
        ## print(number_of_genes_per_row)
        u_split <- split(u, ceiling(seq_along(u) / number_of_genes_per_row))
        rows <- min(length(u_split), rows)

        ## we add an extra row for the title
        heights_vector <- rep(1, rows + 1)
        heights_vector[1] <- lcm(1)
        ## to put each image in the directory of each dataset
        ## new_path = paste(simulation_output_paths[[dataset_index]], "predicted_downregulation/p_j_downregulated_values_999999.png", sep = "")
        ## to gather all the images in out directory
        new_path = paste(patient_folder, "matchings_predictor_output/", basename(simulation_output_paths[[dataset_index]]), ".png", sep = "")
        png(filename = new_path, width = 1920, height = 1200)
        layout(matrix(seq(1, rows + 1), rows + 1, 1), heights = heights_vector)
        old_par <- par(mar = c(0,0,1,0))
        ## old_par <- par(mar = c(0,0,0,0))
        ## browser()
        plot.new()
        title(main = paste("predicted gene expression after downregulation", simulation_output_paths[[dataset_index]]))
        par(mar = c(0,0,1.3,0), mgp = c(3, 0.3, 0))

        for(i in seq_len(rows)){
            ## print(length(u_split[[i]]))
            x_data_split <- x_data[u_split[[i]]]
            y_data_split <- y_data[u_split[[i]]]
            is_the_id_in_gene_ids_split <- is_the_id_in_gene_ids[u_split[[i]]]
            original_y_data_split <- original_y_data[u_split[[i]]]
            genes_skipped_in_distance_based_prediction_split <- genes_skipped_in_distance_based_prediction[u_split[[i]]]
            l <- length(x_data[u_split[[1]]])
            x_lim_split <- c((i - 1) * l, (i - 1) * l + l - 1)
            y_lim_split <- c(0, max(original_y_data_split))

            ## plotting points corresponding to the selected gene_ids
            if(sum(is_the_id_in_gene_ids_split) > 0) {
                x_to_plot <- x_data_split[is_the_id_in_gene_ids_split == T]
                y_to_plot <- y_data_split[is_the_id_in_gene_ids_split == T]
                original_y_to_plot <- original_y_data_split[is_the_id_in_gene_ids_split == T]
                plot(x_to_plot, y_to_plot, xlim = x_lim_split, ylim = y_lim_split, col = "red", pch = 20)
                points(x_to_plot, original_y_to_plot, xlim = x_lim_split, ylim = y_lim_split, col = "red", pch = ".")
            } else {
                plot(c(0, 1), c(0, 1), type = "n", xlim = x_lim_split, ylim = y_lim_split)
            } ## , bty = "n", , xaxt = "n", yaxt = "n", ann = F

            ## plotting the remaining points (if consider_only_specified_gene_ids is FALSE)
            if(consider_only_specified_gene_ids == F && sum(!is_the_id_in_gene_ids_split) > 0) {
                x_to_plot <- x_data_split[is_the_id_in_gene_ids_split == F]
                y_to_plot <- y_data_split[is_the_id_in_gene_ids_split == F]
                original_y_to_plot <- original_y_data_split[is_the_id_in_gene_ids_split == F]
                points(x_to_plot, y_to_plot, xlim = x_lim_split, ylim = y_lim_split, col = "black", pch = 20)
                points(x_to_plot, original_y_to_plot, xlim = x_lim_split, ylim = y_lim_split, col = "black", pch = ".")
            }

            ## plotting vertical bars to mark eventual genes which have been skipped during the distance-based prediction procedure
            if(sum(genes_skipped_in_distance_based_prediction_split)) {
                x_to_plot <- x_data_split[genes_skipped_in_distance_based_prediction_split == T]
                abline(v = x_to_plot, col = "green")
            }
        }
        dev.off()
    }
    for(i in seq_along(simulation_output_paths)) {
        plot_dataset(i)
    }
}

compute_pairwise_distances <- function(simulation_output_paths, gene_ids = NULL, consider_only_specified_gene_ids = F, consider_relative_changes, allow_linear_correction = F, choice = NULL, plot_the_matrix = F)
{
    if(is.null(choice)) {
        scores <- list()
        ## compute_pairwise_distances(simulation_output_paths, gene_ids, consider_only_specified_gene_ids, consider_relative_changes, "ratio")
        ## scores[[length(scores) + 1]] <- compute_pairwise_distances(simulation_output_paths, gene_ids, consider_only_specified_gene_ids, consider_relative_changes, allow_linear_correction, "spearman")
        ## scores[[length(scores) + 1]] <- compute_pairwise_distances(simulation_output_paths, gene_ids, consider_only_specified_gene_ids, consider_relative_changes, allow_linear_correction, "average")
        ## scores[[length(scores) + 1]] <- compute_pairwise_distances(simulation_output_paths, gene_ids, consider_only_specified_gene_ids, consider_relative_changes, allow_linear_correction, "euclidean")
        scores[[length(scores) + 1]] <- compute_pairwise_distances(simulation_output_paths, gene_ids, consider_only_specified_gene_ids, consider_relative_changes, allow_linear_correction, "norm1")
        return(scores)
    }
    dataframes <- lapply(simulation_output_paths,
                         function(x) read.table(paste(x, "predicted_downregulation/p_j_downregulated_values_999999.tsv", sep = ""), header = T, colClasses = c("numeric", "numeric")))

    all_genes <- read.table("processed/gene_id_dictionary.tsv", header = T, colClasses = c("character", "numeric"))$gene_id_cpp
    for(i in 1:length(dataframes)) {
        expressed_genes <- dataframes[[i]]$gene_id
        not_expressed_genes <- setdiff(all_genes, expressed_genes)
        dataframes[[i]] <- rbind(dataframes[[i]], data.frame(gene_id = not_expressed_genes, p_j_downregulated_values = rep(0, length(not_expressed_genes))))
        dataframes[[i]] <- dataframes[[i]][order(dataframes[[i]]$gene_id), ]
    }

    ## if(!is.null(gene_ids)) {
    ##     dataframes <- lapply(dataframes, function(x) x[x$gene_id %in% gene_ids, ])
    ## }
    n <- length(dataframes)
    if(n < 2) {
        print(paste("error: n = ", n, sep = ""))
        stop("abort")
    }
    ## to speed up the function
    plot_only_essentials <- F
    if(n > 15) {
        plot_only_essentials <- T
    }
    if(!plot_only_essentials) {
        new_maximized_device()
        my_widths <- c(lcm(1), rep(1, n))
        my_heights <- c(lcm(1), rep(1, n))
        layout(matrix(1:(n+1)^2, n+1, n+1, byrow = T), widths = my_widths, heights = my_heights)
    }
    distance_matrix <<- matrix(rep(0, n^2), n, n)
    ## if n is too large only the level plots make sense

    my_distance <- function(x_i, x_j, is_the_id_in_gene_ids, consider_only_specified_gene_ids)
    {
        x_i_gene_ids <- x_i[is_the_id_in_gene_ids]
        x_j_gene_ids <- x_j[is_the_id_in_gene_ids]
        x_i_not_gene_ids <- x_i[!is_the_id_in_gene_ids]
        x_j_not_gene_ids <- x_j[!is_the_id_in_gene_ids]
        linear_model_correction <- F
        if(consider_only_specified_gene_ids && allow_linear_correction) {
            x_i_considered <- x_i_gene_ids
            x_j_considered <- x_j_gene_ids

            linear_model_correction <- T
            if(linear_model_correction && length(x_i_not_gene_ids) > 0) {
                df_regulated <- data.frame(x_i = x_i_gene_ids, x_j = x_j_gene_ids)
                df_not_regualted <- data.frame(x_i = x_i_not_gene_ids, x_j = x_j_not_gene_ids)
                model <- lm(x_j ~ x_i, data = df_not_regualted)
                predictions <- predict(model, newdata = df_regulated)
                ## x_j_gene_ids_corrected <- x_j_gene_ids - predictions
                x_j_gene_ids_corrected <- x_j_gene_ids - (predictions - x_i_gene_ids)

                x_j_considered <- x_j_gene_ids_corrected
            }
        } else {
            x_i_considered <- x_i
            x_j_considered <- x_j
        }
        if(choice == "ratio") {
            ## both_zero <- x_i_considered < 10^(-10) & x_j_considered < 10^(-10)
            ## x_i_considered <- x_i_considered[!both_zero]
            ## x_j_considered <- x_j_considered[!both_zero]
            ratios <- sapply(1:length(x_i_considered), function(k)
            {
                min(x_i_considered[k]/x_j_considered[k],
                    x_j_considered[k]/x_i_considered[k])
            })
            ## there is a nan when for an index both x_i_considered and x_j_considered are zero
            ratios <- ratios[!is.na(ratios)]
            if(!plot_only_essentials) {
                par(mar = c(2, 2, 1, 0))
                hist(ratios, main = "")
            }
            return(mean(ratios))
        } else if(choice == "spearman") {
            ## browser()
            correlation <- cor(x_i_considered, x_j_considered, method = "spearman")
            name_i <- basename(simulation_output_paths[[i]])
            name_j <- basename(simulation_output_paths[[j]])
            if(!plot_only_essentials) {
                par(mar = c(0, 0, 2, 0))
                ## red_color <- rgb(1, 0, 0, 0.01)
                red_color <- "red"
                ## old_par <- par(pty = "s")
                plot(x_i_gene_ids, x_j_gene_ids, col = red_color, main = paste("r_s = ", round(correlation, 2), " (", name_i, " vs ", name_j, ")"), cex.main = 0.8, xlab = "", ylab = "", xaxt = "n", yaxt = "n", pch = ".")
                points(x_i_not_gene_ids, x_j_not_gene_ids, col = "black", pch = ".")
                if(linear_model_correction) {
                    abline(model, col = "blue")
                    points(x_i_gene_ids, x_j_gene_ids_corrected, col = "blue", pch = ".")
                    abline(c(0, 0), c(1, 1), col = "green")
                }
                ## par(old_par)
            }
            return(correlation)
        } else if(choice == "average") {
            if(!plot_only_essentials) {
                par(mar = c(2, 2, 1, 0))
                ## hist((x_i_considered - x_j_considered)/length(x_i_considered), main = "")
                plot(density((x_i_considered - x_j_considered)/length(x_i_considered)), main = "")
                abline(v = 0, col = "red")
            }
            return(mean(x_i_considered - x_j_considered))
        } else if(choice == "euclidean") {
            ## browser()
            if(!plot_only_essentials) {
                par(mar = c(2, 2, 1, 0))
                hist((x_i_considered - x_j_considered)^2, main = "")
            }
            return(sqrt(sum((x_i_considered - x_j_considered)^2)))
        } else if(choice == "norm1") {
            if(!plot_only_essentials) {
                par(mar = c(2, 2, 1, 0))
                hist(abs(x_i_considered - x_j_considered), main = "")
            }
            return(sum(abs(x_i_considered - x_j_considered)))
        }
    }
    for(i in 0:n) {
        for(j in 0:n) {
            if(i == j) {
                if(!plot_only_essentials) {
                    par(mar = c(0, 0, 2, 0))
                    plot.new()
                }
                next
            }
            if(i == 0) {
                if(!plot_only_essentials) {
                    par(mar = c(0, 0, 0, 0))
                    plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
                    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = rgb(0.4, 0.4, 0.4, 0.5))
                    s <- basename(simulation_output_paths[[j]])
                    text(x = 0.5, y = 0.5, s, cex = 1.5)
                }
                next
            }
            ## par(mar = c(0, 0, 1, 0))
            if(j == 0) {
                if(!plot_only_essentials) {
                    par(mar = c(0, 0, 0, 0))
                    plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
                    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = rgb(0.4, 0.4, 0.4, 0.5))
                    s <- basename(simulation_output_paths[[i]])
                    text(x = 0.5, y = 0.5, s, cex = 1.5, srt = 90)
                }
                next
            }
            if(plot_only_essentials) {
                if(j > 1) {
                    next
                }
            }
            df_i <- dataframes[[i]]
            df_j <- dataframes[[j]]
            df_i <- df_i[order(df_i$gene_id), ]
            df_j <- df_j[order(df_j$gene_id), ]
            if(consider_relative_changes) {
                x_i <- df_i$p_j_downregulated_values
                x_j <- df_j$p_j_downregulated_values
            } else {
                if(!exists("initial_gene_expression_profile")) {
                    filename <- paste(patient_folder, "gene_expression_profile_exported.tsv", sep = "")
                    initial_gene_expression_profile <<- read.table(filename, header = T, colClasses = c("numeric", "numeric"))
                }
                df_i <- merge(df_i, initial_gene_expression_profile)
                df_i$rpm_downregulated <- df_i$rpm * df_i$p_j_downregulated_values
                x_i <- df_i$rpm_downregulated
                ## df_i$final_rpm <- df_i$rpm - df_i$rpm_downregulated
                ## x_i <- df_i$final_rpm

                df_j <- merge(df_j, initial_gene_expression_profile)
                df_j$rpm_downregulated <- df_j$rpm * df_j$p_j_downregulated_values
                x_j <- df_j$rpm_downregulated
                ## df_j$final_rpm <- df_j$rpm - df_j$rpm_downregulated
                ## x_j <- df_j$final_rpm
            }
            old_consider_only_specified_gene_ids <- consider_only_specified_gene_ids
            if(length(unlist(gene_ids)) == 0 ||
               length(gene_ids[[i]]) == 0 && length(gene_ids[[j]]) == 0) {
                is_the_id_in_gene_ids <- rep(T, length(x_i))
                consider_only_specified_gene_ids <- F
            } else {
                gene_ids_of_the_current_pair <- unique(unlist(c(gene_ids[[i]], gene_ids[[j]])))
                is_the_id_in_gene_ids <- (df_i$gene_id %in% gene_ids_of_the_current_pair)
                is_the_id_in_gene_ids_test <- (df_j$gene_id %in% gene_ids_of_the_current_pair)
                s <- sum(is_the_id_in_gene_ids != is_the_id_in_gene_ids_test)
                if(s > 0) {
                    print(paste("error: sum(is_the_id_in_gene_ids != is_the_id_in_gene_ids_test) = ", s, " > 0", sep = ""))
                    stop("abort")
                }
            }
            distance_matrix[i,j] <<- my_distance(x_i, x_j, is_the_id_in_gene_ids, consider_only_specified_gene_ids)
            consider_only_specified_gene_ids <- old_consider_only_specified_gene_ids
        }
    }
    ## print(round(distance_matrix, 2))
    if(plot_the_matrix) {
        new_maximized_device()
    }
    ## print(simulation_output_paths)
    my_labels <- unlist(lapply(simulation_output_paths, function(x) basename(x)))
    rownames(distance_matrix) <- my_labels
    colnames(distance_matrix) <- my_labels
    m <- distance_matrix
    ## browser()
    if(plot_only_essentials) {
        min_assigned_value <- min(m[2:(dim(m)[[1]]), 1])
        max_assigned_value <- max(m[2:(dim(m)[[1]]), 1])
    } else {
        m_without_diag <- m
        diag(m_without_diag) <- NA
        min_assigned_value <- min(m_without_diag, na.rm = T)
        max_assigned_value <- max(m_without_diag, na.rm = T)
    }
    diag(m) <- 1
    if(choice == "ratio") {
        color_order <- c("black", "white")
        if(plot_only_essentials) {
            m[, 2:(dim(m)[[2]])]  <- max_assigned_value
        }
        diag(m) <- max_assigned_value
    } else if(choice == "spearman") {
        color_order <- c("white", "black", "white")
        if(plot_only_essentials) {
            m[, 2:(dim(m)[[2]])]  <- max_assigned_value
        }
        diag(m) <- max_assigned_value
    } else if(choice == "average") {
        color_order <- c("blue", "white", "red")
    } else if(choice == "euclidean") {
        color_order <- c("white", "black")
    } else if(choice == "norm1") {
        color_order <- c("white", "black")
    }
    rgb.palette <- colorRampPalette(color_order, space = "rgb")
    my_colors <- rgb.palette(120)
    if(choice == "ratio") {
        my_positions1 <- seq(0, 1, length.out = 100)
    } else if(choice == "spearman") {
        my_positions1 <- seq(-1, 1, length.out = 100)
    } else if(choice == "average") {
        my_positions1 <- seq(min(m, -max(m)), max(m, -min(m)), length.out = 100)
    } else if(choice == "euclidean") {
        my_positions1 <- seq(0, max(m), length.out = 100)
    } else if(choice == "norm1") {
        my_positions1 <- seq(0, max(m), length.out = 100)
    }

    if(choice == "spearman") {
        color_order <- c("black", "white")
        rgb.palette <- colorRampPalette(color_order, space = "rgb")
        my_colors <- rgb.palette(120)
    }

    if(choice == "average") {
        ## my_positions2 <- my_positions1
        my_positions2 <- seq(min(m, -max(m)), max(m, -min(m)), length.out = 100)
    } else {
        m_without_diag <- m ## unused
        diag(m) <- NA
        my_positions2 <- seq(min(m, na.rm = T), max(m, na.rm = T), length.out = 100)
        if(length(unique(my_positions2)) == 1) {
            if(max(m, na.rm = T) != 0) {
                my_positions2 <- seq(0, max(m, na.rm = T), length.out = 100)
            } else {
                ## remember that the min is equal to max(m)
                my_positions2 <- seq(min(m + diag(max(m, na.rm = T), n)[, n:1], na.rm = T), 1, length.out = 100)
            }
        }
    }
    m <- m[, n:1]
    ## print(m)

    lattice.options(axis.padding=list(factor=0.5))
    my_levelplot <- function(my_positions) {
        obj <- levelplot(m, main = paste("matrix for distance \"", choice, "\"", sep = "" ),
                         xlab = "", ylab = "",
                         cuts = 100,
                         col.regions = my_colors,
                         at = my_positions,
                         colorkey = list(col = my_colors, at = my_positions),
                         scales = list(y = list(rot = 0), x = list(rot = 45), cex = 0.5)
                         ## panel = function(y, x, ...) {
                         ##     ltext(x = x, y = y, labels = round(m, 1), cex = 1, font = 2,
                         ##           fontfamily = "HersheySans")
                         ## }
                         )
        return(obj)
    }
    if(plot_the_matrix) {
        plot1 <- my_levelplot(my_positions1)
        plot2 <- my_levelplot(my_positions2)
        grid.arrange(plot1, plot2, ncol = 2)
    }

    ## we adjust the scores so that sorting them ascendignly we get the mirnas contributing the most
    unadjusted_scores <- m[2:n,n]
    if(choice == "ratio") {
        adjusted_scores <- unadjusted_scores
    } else if(choice == "spearman") {
        adjusted_scores <- abs(unadjusted_scores)
    } else if(choice == "average") {
        adjusted_scores <- 0 - abs(unadjusted_scores)
    } else if(choice == "euclidean") {
        adjusted_scores <- 0 - abs(unadjusted_scores)
    } else if(choice == "norm1") {
        adjusted_scores <- 0 - abs(unadjusted_scores)
    }
    return(adjusted_scores)
}

plot_rankings <- function(list_of_patients, list_of_rankings)
{
    if(length(list_of_patients) != length(list_of_rankings)) {
        print(paste("error: length(list_of_patients = ", length(list_of_patients, ", length(list_of_rankings) = ", length(list_of_rankings), sep = "")))
        stop("abort")
    }
    first_n_mirnas_to_consider <- 10
    rankings_count <- length(list_of_rankings)
    rankings_size <- length(list_of_rankings[[1]])
    new_maximized_device()
    layout(matrix(1:(rankings_count*rankings_size), rankings_count, rankings_size, byrow = T))
    for(i in 1:rankings_count) {
        rankings <- list_of_rankings[[i]]
        for(j in 1:rankings_size) {
            ## new code
            scores <- rankings[[j]]
            substituted <- sub("p__", "", names(scores), perl = T)
            substituted <- sub("__[ar][0-9]+__", "", substituted, perl = T)
            substituted <- as.numeric(substituted) + 1
            mirna_expression_profile_path <- paste("patients/", list_of_patients[[i]], "/mirna_expression_profile.tsv", sep = "")
            mirna_expression_profile <- read.table(mirna_expression_profile_path, header = T, colClasses = c("character", "numeric", "numeric"))
            mirna_expression_profile <- mirna_expression_profile[order(-mirna_expression_profile$reads), ]
            o <- order(scores)[1:first_n_mirnas_to_consider]
            ## browser()
            ranking_of_perturbed_indexes <- substituted[o]
            best_mirnas <- mirna_expression_profile$mirbase_id[ranking_of_perturbed_indexes]
            best_scores <- scores[o]
            print(best_mirnas)
            par(mar = c(0, 4, 1, 2), las = 1)
            plot(1:length(best_scores), best_scores, ylab = "")
            ## browser()

            ## old code
            ## scores <- rankings[[j]]
            ## ## note the "- 1", in C++ these values are 0-based
            ## o <- order(scores)[1:first_n_mirnas_to_consider] - 1
            ## s <- sort(scores)[1:first_n_mirnas_to_consider]
            ## ## par(mfrow = c(2,1))
            ## ranking_string <- paste(o, collapse = ", ")
            ## plot(o, s, xaxt = "n", xlab = "n-th strongest miRNA perturbed", ylab = "score", main = ranking_string)
            ## lines(o, s)
            ## abline(v = o, col = "gray", lty = "dashed")
            ## axis(1, at = o, labels = str(o))
            ## ## plot(s, xaxt = "n", main = "ordered scores", xlab = "miRNA perturbed", ylab = "score")
            ## ## lines(1:first_n_mirnas_to_consider, s)
            ## ## axis(1, at = 1:first_n_mirnas_to_consider, labels = paste(o))
            ## ## par(mfrow = c(1,1))
        }
    }
    end_browser_plot()
}

analyze_mirna_expression_profiles_of_perturbed_data <- function(simulation_output_paths)
{
    print("analyzing mirna expression profiles of perturbed data")
    mirna_ids_ordered <- NULL
    first_reference_y_data <- NULL
    last_reference_y_data <- NULL
    compare_mirna_expression_profile_to_reference <- function(simulation_output_path)
    {
        expression_profile_folder <- paste(simulation_output_path, "mirna_expression_profiles", sep = "")
        files <- list.files(path = expression_profile_folder, pattern = "*.tsv", full.names = T, recursive = F)

        first_file <- files[[1]]
        last_file <- files[[length(files)]]
        first_df <- read.table(first_file, header = T, colClasses = c("numeric", "numeric", "numeric"))
        last_df <- read.table(last_file, header = T, colClasses = c("numeric", "numeric", "numeric"))
        ## browser()
        if(is.null(first_reference_y_data)) {
            my_order <- order(first_df$reads)
            assign("mirna_ids_ordered", value = first_df$mirna_id[my_order], envir = parent.frame())
            assign("first_reference_y_data", value = first_df$reads[my_order], envir = parent.frame())

            my_order <- match(mirna_ids_ordered, last_df$mirna_id)
            assign("last_reference_y_data", value = last_df$reads[my_order], envir = parent.frame())
            ## browser()
            ## browser()
        } else {
            my_order <- match(mirna_ids_ordered, first_df$mirna_id)
            first_y_data <- first_df$reads[my_order]
            my_order <- match(mirna_ids_ordered, last_df$mirna_id)
            last_y_data <- last_df$reads[my_order]

            new_maximized_device()
            x_data <- 1:length(first_reference_y_data)
            y_lim <- c(0, max(first_reference_y_data))
            simulation_id <- basename(simulation_output_path)

            plot(x_data, first_reference_y_data, main = paste(simulation_id, "miRNA expression, reads"), ylim = y_lim, col = "black", pch = 20)
            points(x_data, last_reference_y_data, ylim = y_lim, col = "blue", pch = 20)
            points(x_data, first_y_data, ylim = y_lim, col = "red")
            points(x_data, last_y_data, ylim = y_lim, col = "green", pch = ".")
        }
    }

    if(length(simulation_output_paths) < 2) {
        stop(paste("length(simulation_output_paths) =", length(simulation_output_paths)))
    }
    for(simulation_output_path in simulation_output_paths) {
        ## print(simulation_output_path)
        compare_mirna_expression_profile_to_reference(simulation_output_path)
    }
}

patient_id <- "artificial_TCGA-CJ-4642"
## patient_id <- "artificial0"
## patient_id <- "artificial1"
patient_folder <- paste("patients/", patient_id, "/", sep = "")

simulation_id <- "original_data"
## simulation_id <- "______distance"
simulation_output_path <- paste(patient_folder, "matchings_predictor_output/", simulation_id, "/", sep = "")
close_all_devices()
## analyze_convergence(simulation_output_path)
## generate_readable_dynamics_log(simulation_output_path)
## analyze_mirna_expression_profiles(simulation_output_path)
## analyze_cluster_expression_profiles(simulation_output_path)
## analyze_dynamics_for_small_interaction_graphs(simulation_output_path)
## analyze_probabilities(patient_folder, simulation_output_path)

simulation_output_paths <- c()

## hela_patients <- c("artificial_TCGA-CJ-4642", "artificial_ENCFF360IHM-hela", "artificial_ENCFF495ZXC-hela", "artificial_ENCFF612ZIR-hela", "artificial_ENCFF729EQX-hela", "artificial_ENCFF806EYY-hela", "artificial_ENCFF902KUU-hela")
## hela_patients <- c("artificial_TCGA-CJ-4642", "artificial_ENCFF360IHM-hela", "artificial_ENCFF495ZXC-hela")

hela_patients <- c("artificial_ENCFF360IHM-hela", "artificial_ENCFF495ZXC-hela", "artificial_ENCFF612ZIR-hela", "artificial_ENCFF729EQX-hela", "artificial_ENCFF806EYY-hela", "artificial_ENCFF902KUU-hela")
## hela_patients <- c("artificial_ENCFF360IHM-hela", "artificial_ENCFF495ZXC-hela", "artificial_ENCFF612ZIR-hela")
## hela_patients <- c("artificial_ENCFF360IHM-hela", "artificial_ENCFF495ZXC-hela")
## hela_patients <- c("artificial_ENCFF360IHM-hela")

rankings <- list()

for(patient in hela_patients) {
    simulation_output_paths_for_current_patient <- c()
    patient_folder_for_current_patient <- paste("patients/", patient, "/", sep = "")
    ## list_of_perturbed_data <- c("original_data", "p__0__a500000__", "p__1__a500000__", "p__2__a500000__")
    list0 <- list.files(path = paste(patient_folder_for_current_patient, "matchings_predictor_output/", sep = ""), pattern = "original_data$", full.names = T, recursive = F)
    list1 <- list.files(path = paste(patient_folder_for_current_patient, "matchings_predictor_output/", sep = ""), pattern = "__a500000__$", full.names = T, recursive = F)
    list_of_perturbed_data <- c(list0, list1)
    ## browser()
    for(perturbed_data in list_of_perturbed_data) {
        ## new_path <- paste(patient_folder_for_current_patient, "matchings_predictor_output/", perturbed_data, "/", sep = "")
        new_path <- paste(perturbed_data, "/", sep = "")
        simulation_output_paths <- c(simulation_output_paths, new_path)
        simulation_output_paths_for_current_patient <- c(simulation_output_paths_for_current_patient, new_path)
    }
    ## browser()
    analyze_predictions_of_perturbed_data(patient_folder_for_current_patient,
                                          simulation_output_paths_for_current_patient,
                                          gene_ids = sapply(1:length(simulation_output_paths_for_current_patient), function(x) list()),
                                          consider_only_specified_gene_ids = F)

    ## rankings[[length(rankings) + 1]] <- compute_pairwise_distances(simulation_output_paths_for_current_patient,
    ##                                                                gene_ids = sapply(1:length(simulation_output_paths_for_current_patient), function(x) list()),
    ##                                                                consider_only_specified_gene_ids = T,
    ##                                                                consider_relative_changes = F,
    ##                                                                allow_linear_correction = F,
    ##                                                                plot_the_matrix = F)
}

plot_rankings(hela_patients, rankings)

## simulation_output_paths <- c(simulation_output_paths, simulation_output_path)
## simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", "g__87__r3__", "/", sep = ""))
## simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", "p__0__r100__", "/", sep = ""))
## simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", "p__0__r10__", "/", sep = ""))
## simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", "p__0__r1__", "/", sep = ""))

## simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", "______", "/", sep = ""))
## simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", "______distance", "/", sep = ""))
## ## simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", "p__0__r0.1__", "/", sep = ""))
## ## simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", "p__0__r-0.5__", "/", sep = ""))
## simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", "p__0__r-1__", "/", sep = ""))
## simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", "p__0__r-0.9__", "/", sep = ""))

## perturbed_datasets_count <- 10
## replicates <- 3
## for(i in seq_len(perturbed_datasets_count)) {
##     for(j in seq_len(replicates)) {
##         simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/g__86__r", i, "__", j - 1, "/", sep = ""))
##     }
## }

## ## first_n_mirnas_to_perturb <- 87
## first_n_mirnas_to_perturb <- 7
## ## delta <- 22
## delta <- 0
## for(i in seq(0 + delta, first_n_mirnas_to_perturb - 1 + delta)) {
##     ## for(j in c(0, 1)) {
##     for(j in c(0)) {
##         if(j == 0) {
##             simulation_id <- paste("p__", i, "__a500000__", sep = "")
##         } else {
##             simulation_id <- paste("p__", i, "__r-1__", sep = "")
##         }
##     }
##     simulation_output_paths <- c(simulation_output_paths, paste(patient_folder, "matchings_predictor_output/", simulation_id, "/", sep = ""))
## }

## ## for(path in simulation_output_paths) {
## ##     analyze_convergence(path)
## ## }

## mirna_expression_profile_filename <- paste(patient_folder, "mirna_expression_profile_exported.tsv", sep = "")
## mirnas <- read.table(mirna_expression_profile_filename, header = T, colClasses = c("numeric", "numeric"))
## mirnas_considered <- mirnas[order(-mirnas$rpm),][1:first_n_mirnas_to_perturb, "mirna_id"]
## gene_ids <- list()
## gene_ids[[length(gene_ids) + 1]] <- list()
## ## for(i in seq_len(perturbed_datasets_count * replicates)) {
## ##     gene_ids[[length(gene_ids) + 1]] <- list()
## ## }
## ## no mirna is perturbed for original data and gaussian data
## gene_ids[[1]] <- list()
## gene_ids[[2]] <- list()
## gene_ids[[3]] <- list()
## gene_ids[[4]] <- list()
## gene_ids[[5]] <- list()
## gene_ids[[6]] <- list()
## i <- length(gene_ids) + 1
## for(mirna_id in mirnas_considered) {
##     filename <- paste("interactions/mirna_id", mirna_id, "_interactions.rds", sep = "")
##     if(!file.exists(filename)) {
##         ## I should have used a SQL database here...
##         print(paste("generating target file for mirna_id =", mirna_id))
##         source("analyze_expression_profiles.r")
##         get_targets_for_mirna(mirna_id)
##         print("generated")
##     }
##     gene_ids_for_mirna <- readRDS(filename)
##     gene_ids[[i]] <- gene_ids_for_mirna
##     gene_ids[[i + 1]] <- gene_ids_for_mirna
##     ## gene_ids <- unique(c(gene_ids, gene_ids_for_mirna))
##     i <- i + 2
## }
## gene_ids[[2]] <- gene_ids[[5]]
## gene_ids[[3]] <- gene_ids[[5]]
## gene_ids[[4]] <- gene_ids[[5]]
## gene_ids[[4]] <- gene_ids[[6]]
## gene_ids[[5]] <- gene_ids[[6]]

## analyze_predictions_of_perturbed_data(patient_folder, simulation_output_paths, gene_ids = gene_ids, consider_only_specified_gene_ids = F)

## rankings[[length(rankings) + 1]] <- compute_pairwise_distances(simulation_output_paths, gene_ids = gene_ids, consider_only_specified_gene_ids = F, consider_relative_changes = T)
## rankings[[length(rankings) + 1]] <- compute_pairwise_distances(simulation_output_paths, gene_ids = gene_ids, consider_only_specified_gene_ids = T, consider_relative_changes = T, allow_linear_correction = T)
## rankings[[length(rankings) + 1]] <- compute_pairwise_distances(simulation_output_paths, gene_ids = gene_ids, consider_only_specified_gene_ids = T, consider_relative_changes = T, allow_linear_correction = F)

## analyze_mirna_expression_profiles_of_perturbed_data(simulation_output_paths)
## rankings[[length(rankings) + 1]] <- compute_pairwise_distances(simulation_output_paths, gene_ids = gene_ids, consider_only_specified_gene_ids = F, consider_relative_changes = F, allow_linear_correction = F)
## rankings[[length(rankings) + 1]] <- compute_pairwise_distances(simulation_output_paths, gene_ids = gene_ids, consider_only_specified_gene_ids = T, consider_relative_changes = F, allow_linear_correction = F)

## compute_pairwise_distances(simulation_output_paths)
end_browser_plot()
