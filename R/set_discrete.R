#' TODO input validation 
#' helper function which is called in most of the clustering functions
#' function calls
#' @import 
#' @param mat_to_plot
#' @param continuous , what to plot, user input

#' @return boolean
#' @export


set_discrete <- function(continuous, mat_to_Plot){
    if (continuous=="discrete"){
      #Convert text values to numeric
      rNames <- rownames(mat_to_Plot)
      mat_to_Plot <- apply(mat_to_Plot, 2, as.numeric)
      rownames(mat_to_Plot) <- rNames

      mat_to_Plot[mat_to_Plot<0]=-1
      mat_to_Plot[mat_to_Plot>0]=1

      isDiscrete = T
    }else{
      isDiscrete = F
    }

    return(isDiscrete)
}