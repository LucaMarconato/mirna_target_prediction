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
    ## return(c(12.8, 7.5))
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
