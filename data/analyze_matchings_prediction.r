source("my_device.r")

analyze_convergence <- function(patient_folder)
{    
    filename <- paste(patient_folder, "matchings_prediction.tsv", sep = "")
    a <<- read.table(filename, header = T, colClasses = c("numeric", "numeric", "numeric", "numeric", "numeric"))
    x_values <- seq(1, length(a[[1]]))
    
    new_maximized_device()
    layout(matrix(1:6, 2, 3, byrow = T))
    plot(x_values, a[[1]], main = "mirna cumulative scaling")
    plot(x_values, a[[2]], main = "cluster cumulative scaling")
    plot(x_values, a[[3]], main = "average mirna level")
    plot(x_values, a[[4]], main = "average cluster level")
    plot(x_values, a[[5]], main = "total exchange in the step")
    plot(x_values, cumsum(a[[5]]), main = "cumulative total exchange")
}

patient_id <- "TCGA-CJ-4642"
patient_folder <- paste("patients/", patient_id, "/", sep = "")
close_all_devices()
analyze_convergence(patient_folder)
