#' TODO input validation 
#' function calls collapse_paths, plot_grid, plotting
#' @import 
#' @param  pheno 
#' @param samplesID
#' @param KEGG_MAT
#' @param reduced_kegg_hierarchy
#' @param level selected kegg level to plot, user input
#' @param continuous what/ how to plot, user input
#' @param doGrouping group plot?, user input
#' @param aspectRatio plotting parameters, user input
#' @return list
#' @export


reset_cluster <- function (pheno, samplesID, KEGG_MAT, reduced_kegg_hierarchy, level, continuous, doGrouping, aspectRatio) {
    #indicate to front end to show loading
    #not included in here

    exp_ann = pheno
    clust_mat = NULL

    #based on the information contained in samplesID we decide what to plot
    #the whole kegg matrix or only the specified rows
    if("All" %in% samplesID){
      
      toPlotMap = KEGG_MAT
    }else{
      toPlotMap = KEGG_MAT[rownames(KEGG_MAT) %in% samplesID,]
    }

    #replace with lucas function
    toPlot = plotting(reduced_kegg_hierarchy,  toPlotMap, level, clust_mat, pheno, samplesID, continuous, doGrouping, aspectRatio)
    
    
    print("after PLOTS in reset cluster")

    hls <- NULL
    cls <- NULL
    clustered <- 1

    #return exp_ann?
    #not sure if hls, cls & clustered will be needed later one
    return (c(toPlot, hls, cls, clustered, exp_ann))
    #return(toPlot)
  }