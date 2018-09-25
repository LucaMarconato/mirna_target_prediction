source("my_device.r")
library(purrr)
library(raster)
library(RColorBrewer)

analyze_convergence <- function(patient_folder)
{    
    filename <- paste(patient_folder, "matchings_prediction.tsv", sep = "")
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

analyze_mirna_expression_profiles <- function(patient_folder)
{
    print("analyzing mirna expression profiles")
    expression_profile_folder <- paste(patient_folder, "mirna_expression_profiles", sep = "")
    expression_profile_copy_folder <- paste(patient_folder, "mirna_expression_profiles_copy", sep = "")
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
        t <- read.table(file, header = T, colClasses = c("numeric", "numeric"))
        if(is.null(mirna_ids_ordered)) {
            my_order <- order(t$relative_expression)
            mirna_ids_ordered <- t$mirna_id[my_order]
        }
        new_path <- tools::file_path_sans_ext(file)
        new_path <- paste(new_path, ".png", sep = "")
        png(filename = new_path, width = 1920, height = 1080)
        layout(matrix(1:3, 3, 1))
        timestep <- strsplit(file, "mirna_expression_profile_")[[1]][2]
        timestep <- strsplit(timestep, ".tsv")[[1]][1]
        my_order <- match(mirna_ids_ordered, t$mirna_id)
        x_data <- 1:length(t$relative_expression)
        y_data <- t$relative_expression[my_order]
        mirna_dynamics[[length(mirna_dynamics) + 1]] <- y_data
        if(is.null(original_y_data)) {
            original_y_data <- y_data
        }
        plot(x_data, y_data, main = paste(timestep, "miRNA expression"), ylim = c(0, max(original_y_data)))
        points(x_data, original_y_data, ylim = c(0, max(original_y_data)), col = "red", pch = "-")
        
        plot(x_data, log10(y_data), main = paste(timestep, "miRNA expression (log10)"), ylim = c(-20, 0))
        points(x_data, log10(original_y_data), ylim = c(-20, 0), col = "red", pch = "-")

        if(i > 0) {
            previous_file <- files[i]
            previous_t <- read.table(previous_file, header = T, colClasses = c("numeric", "numeric"))
            previous_my_order <- match(mirna_ids_ordered, previous_t$mirna_id)
            previous_y_data <- previous_t$relative_expression[previous_my_order]
            difference <- previous_y_data - y_data
            log_difference <- log10(difference)
            ## log_difference[log_difference == -Inf] = 0
            if(sum(is.nan(log_difference)) > 0) {
                print("warning: found a NaN, probably this is due by having reached the machine precision during the simulation")
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

analyze_cluster_expression_profiles <- function(patient_folder)
{
    print("analyzing cluster expression profiles")
    expression_profile_folder <- paste(patient_folder, "cluster_expression_profiles", sep = "")
    expression_profile_copy_folder <- paste(patient_folder, "cluster_expression_profiles_copy", sep = "")
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
        if(i <= 2) {
            print(paste("skipping timestep", i))
            next
        }
        print(paste(round(i/min(file_count, max_timesteps)*100), "%", sep = ""))
        file <- files[i]
        t <- read.table(file, header = T, colClasses = c("numeric", "numeric"))
        if(is.null(cluster_addresses_ordered)) {
            my_order <- order(t$relative_expression)
            cluster_addresses_ordered <- t$cluster_address[my_order]
        }
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

        ## probably useless
        par(mar = old_par)
        dev.off()
    }
    filename <- paste(expression_profile_copy_folder, "/animation.gif", sep = "")
    print(paste("generating", filename))
    system(paste("convert -delay 10 -loop 0 ", expression_profile_copy_folder, "/*.png ", filename, sep = ""))
    print("generating .gif")
    system(paste("open ", expression_profile_copy_folder, "/animation.gif", sep = ""))
}

generate_readable_dynamics_log <- function(patient_folder)
{
    ## note that we are reading data from the cluster_expression_profiles folder and not the cluster_expression_profiles_copy folder
    ## the animated .gifs are made out of the copied folder, to avoid accidental deletion by the C++ code, so be sure that you have not called
    ## the C++ code after having generated the animations, otherwise this function will agree with the data just produced by the C++ code and maybe not
    ## with the animations, unless you have generated them again

    extract_dynamics <- function(element)
    {
        expression_profile_folder <- paste(patient_folder, element, "_expression_profiles", sep = "")
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
    filename <- paste(patient_folder, "dynamics.txt", sep = "")
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

analyze_dynamics_for_small_interaction_graphs <- function(patient_folder)
{
    print("loading the interaction matrix")
    expression_profile_folder <- paste(patient_folder, "mirna_expression_profiles", sep = "")
    files <- list.files(path = expression_profile_folder, pattern = "*.tsv", full.names = T, recursive = F)
    t <- read.table(files[1], header = T, colClasses = c("numeric", "numeric"))
    my_order <- order(t$relative_expression)
    mirna_ids_ordered <<- t$mirna_id[my_order]

    expression_profile_folder <- paste(patient_folder, "cluster_expression_profiles", sep = "")
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

    matrix_filename <- paste(patient_folder, "interaction_matrix.mat", sep = "")
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

patient_id <- "TCGA-CJ-4642"
## patient_id <- "artificial0"
## patient_id <- "artificial1"
patient_folder <- paste("patients/", patient_id, "/", sep = "")
close_all_devices()
analyze_convergence(patient_folder)
generate_readable_dynamics_log(patient_folder)
analyze_mirna_expression_profiles(patient_folder)
analyze_cluster_expression_profiles(patient_folder)
analyze_dynamics_for_small_interaction_graphs(patient_folder)
