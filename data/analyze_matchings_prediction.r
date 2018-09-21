source("my_device.r")
library(purrr)

analyze_convergence <- function(patient_folder)
{    
    filename <- paste(patient_folder, "matchings_prediction.tsv", sep = "")
    a <<- read.table(filename, header = T, colClasses = c("numeric", "numeric", "numeric", "numeric", "numeric"))
    x_values <- seq(1, length(a[[1]]))
    
    new_maximized_device()
    layout(matrix(c(1,2,3,1,2,4,5,7,9,6,8,10), 4, 3, byrow = T))
    plot(x_values, a[[1]], main = "mirna cumulative scaling", pch = 20)
    grid()
    plot(x_values, a[[2]], main = "cluster cumulative scaling", pch = 20)
    grid()
    plot(x_values, a[[3]], main = "average mirna level", pch = 20)
    grid()
    plot(x_values, log10(a[[3]]), main = "average mirna level (log10)", pch = 20)
    grid()
    plot(x_values, a[[4]], main = "average cluster level", pch = 20)
    grid()
    plot(x_values, log10(a[[4]]), main = "average cluster level (log10)", pch = 20)
    grid()
    plot(x_values, a[[5]], main = "total exchange in the step", pch = 20)
    grid()
    plot(x_values, log10(a[[5]]), main = "total exchange in the step (log10)", pch = 20)
    grid()
    plot(x_values, cumsum(a[[5]]), main = "cumulative total exchange", pch = 20)
    grid()
    plot(x_values, log10(cumsum(a[[5]])), main = "cumulative total exchange (log10)", pch = 20)
    grid()
}

analyze_mirna_expression_profiles <- function(patient_folder) {
    expression_profile_folder <- paste(patient_folder, "mirna_expression_profiles", sep = "")
    files <- list.files(path = expression_profile_folder, pattern = "*.tsv", full.names = T, recursive = F)

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
    
    system(paste("convert -delay 10 -loop 0 ", expression_profile_folder, "/*.png ", expression_profile_folder, "/animation.gif", sep = ""))
    system(paste("open ", expression_profile_folder, "/animation.gif", sep = ""))
}

## patient_id <- "TCGA-CJ-4642"
patient_id <- "artificial0"
patient_folder <- paste("patients/", patient_id, "/", sep = "")
close_all_devices()
analyze_convergence(patient_folder)
analyze_mirna_expression_profiles(patient_folder)
