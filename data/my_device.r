library(tools)

## some of the code in this file works only on my computer.
## setting plot_to_browser_enabled <- F and modifying the resolution which appears in get_screen_physical_size() should make this code working for you
## anyway, it is easy to modify this code so that it will work also in your machine.

plot_to_browser_enabled <- T
## my_format is used only if plot_to_browser_enabled is TRUE
## my_format <- "png"
my_format <- NULL

## get_screen_resolution <- function() 
## {
##     return(c(2560, 1440))
##     ## or automatically
##     ## xdpyinfo_path <- "/opt/X11/bin/xdpyinfo"    
##     ## cmd <- sprintf("%s | grep dimensions | perl -pe 's/^.*?([0-9]+x[0-9]+).*/$1/g' | tr 'x' ' '", xdpyinfo_path)
##     ## output <- system(cmd, intern = T, ignore.stderr = T)
##     ## return(as.numeric(unlist(strsplit(output, split = " "))))
## }

## returns the screen size in inches
get_screen_physical_size <- function() 
{
    ## return(c(24 + 1/4 + 1/8, 13 + 1/32))
    ## return(c(30 + 1/4 + 1/8, 17 + 1/32))
    ## return(c(18.85, 10)) ## full hd monitor
    return(c(13, 7.85)) ## built-in retina 15" monitor
    ## or automatically, but I preferred to do it manually via a "binary search"
    ## cmd <- sprintf("%s | grep dimensions | perl -pe 's/^.*?(\\([0-9]+x[0-9]+)\ millimeters.*/$1/g' | tr -d '(' | tr 'x' ' '", xdpyinfo_path)
    ## output <- system(cmd, intern = T, ignore.stderr = T)
    ## v <- as.numeric(unlist(strsplit(output, split = " ")))
    ## return(v/25.4)
}

new_maximized_device <- function(output_format = "png")
{
    ## resolution <- get_screen_resolution()
    physical_size <- get_screen_physical_size()
    ## for some reasons the ratio of the result I get is different from the ratio of the resolution so I am using only one value
    ## physical_size[2] <- physical_size[1] * resolution[2]/resolution[1]
    ## the option unit or units = "px" does not work for me, so I am using inches
    ## print(physical_size)
    if(!plot_to_browser_enabled) {
        dev.new(width = physical_size[1], height = physical_size[2])
    } else {
        end_browser_plot()
        my_format <<- output_format
        new_browser_plot()
    }
}

my_system <- function(...)
{
  stopifnot(!any(names(list(...)) %in% "intern"))
  result <- base::system(..., intern = TRUE)
  return(result)
}

new_browser_plot <- function()
{
    ## old code, used with new_plot.sh (a script that I made for showing plots in the browser and update the view automatically, not used now)
    ## png(filename = "~/programming/web/r_image/plots/plot.png", width = res[1], height = res[2], units = "px", pointsize = 12)

    file_index <- my_system(paste("sh my_plot_name.sh", my_format))
    file_index <- as.numeric(file_index) + 1
    file_index <- sprintf("%05d", file_index)    
    full_path <- paste("~/programming/web/r_image/plots/plot", file_index, ".", my_format, sep = "")
    
    if(my_format == "pdf") {
        pdf_res <- c(16, 9.95)
        pdf(file = full_path, width = 10.5, height = 6.5)   
    } else if(my_format == "png") {
        res <- c(2560, 1600)
        ## res <- c(1920, 1080)
        png(filename = full_path, width = res[1], height = res[2], units = "px", pointsize = 12, res = 220)
    }
}

end_browser_plot <- function() 
{
    ## old code
    ## system("sh ~/programming/web/r_image/new_plot.sh ~/programming/web/r_image/plots/plot.png")
    if(plot_to_browser_enabled) {
        if(length(dev.list()) > 0) {
            for(l in dev.list()) {
                dev.off(l)
                file_index <- my_system(paste("sh my_plot_name.sh", my_format))
                file_index <- as.numeric(file_index)
                file_index <- sprintf("%05d", file_index)
                filename <- paste("plot", file_index, ".", my_format, sep = "")
                full_path <- paste("~/programming/web/r_image/plots/", filename, sep = "")
                if(my_format == "pdf") {
                    ## slow and no performance improvement
                    ## compactPDF(full_path, gs_quality = "screen")
                    cmd <- paste("open -a Safari.app", full_path)
                    system(cmd)   
                } else {
                    cmd = paste("~/programming/web/r_image/display_latest_plot.py", filename)
                    system(cmd)
                }
            }
        }   
    }
}

close_all_devices <- function() 
{
    if(length(dev.list()) > 0) {
        for(l in dev.list()) {
            dev.off(l)            
        }
    }
}
