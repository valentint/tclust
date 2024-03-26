.onAttach <- function(lib, pkg){

    where <- match(paste("package:", pkg, sep = ""), search())
    ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
    ver <- as.character(ver)
    title <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Title")
    title <- as.character(title)
    packageStartupMessage(paste(title, " (version ", ver, ")\n", sep = ""))


##packageStartupMessage(" ________  ________  ________  ________  ___       ___  ___  ________      
##|\\   __  \\|\\   __  \\|\\   __  \\|\\   ____\\|\\  \\     |\\  \\|\\  \\|\\   ____\\     
##\\ \\  \\|\\  \\ \\  \\|\\  \\ \\  \\|\\ /\\ \\  \\___|\\ \\  \\    \\ \\  \\\\\\  \\ \\  \\___|_    
## \\ \\   _  _\\ \\  \\\\\\  \\ \\   __  \\ \\  \\    \\ \\  \\    \\ \\  \\\\\\  \\ \\_____  \\   
##  \\ \\  \\\\  \\\\ \\  \\\\\\  \\ \\  \\|\\  \\ \\  \\____\\ \\  \\____\\ \\  \\\\\\  \\|____|\\  \\  
##   \\ \\__\\\\ _\\\\ \\_______\\ \\_______\\ \\_______\\ \\_______\\ \\_______\\____\\_\\  \\ 
##    \\|__|\\|__|\\|_______|\\|_______|\\|_______|\\|_______|\\|_______|\\_________\\
##                                                               \\|_________|
##                                                                           
##Universidad de Valladolid (UVa) - https://www.uva.es
##")

}
