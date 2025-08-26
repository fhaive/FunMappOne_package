#' TODO input validation & output validation 
#' called by do_cluster
#' function calls 
#' @import 
#' @param x , matrix
#' @return transformed matrix #TODO check if this is correct and works
#' @export


range01 <- function(x){
    (x-min(x))/(max(x)-min(x))
    return(x)
    }