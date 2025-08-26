#' @import grid
#' @param toPlotGenes output of plot grid genes/ do_gene_heat_map
#' @param kegg_nano_1 output of collapse_paths/ output of gene_heat_map
#' @param mat_to_Plot output of do_gene_heat_map
#' @param exp_ann output of do_cluster
#' @return plot
#' @export


grid_draw <- function(toPlotGenes, kegg_nano_1,  mat_to_Plot, exp_ann){

    

    nSamples = nrow(kegg_nano_1[[1]])
    length_gene  = max(unlist(lapply(X = colnames(mat_to_Plot),FUN = nchar)))
    nSamplesGene = nrow(mat_to_Plot)
    nGenes = ncol(mat_to_Plot)
    # #Get from list of grobs
    # print("class(gVars$toPlotGenes$gplotGene)")
    # print(class(gVars$toPlotGenes$gplotGene))
    # #print(gVars$toPlotGenes)
    # print(as.ggplot(gVars$toPlotGenes$gplotGene))
    width = function(nSamples, length_gene, nSamplesGene){
    if(!is.null(nSamples)){
      mwidth = (length_gene * 15) + (nSamplesGene * 20)
      print(paste("MY sample IS ---->", nSamplesGene))
      print(paste("MY path len IS ---->", length_gene))
      print(paste("MY WIDTH IS ---->", mwidth))

      mwidth = max(600, mwidth)
      return(mwidth)
    }else{
      return("auto")
    }
  }
  height = function(nGenes, exp_ann){
    if(!is.null(nGenes)){
      
      print("exp_ann")
      print(exp_ann)

      #if(is.null(exp_ann) || nrow(exp_ann)==0){
      #  exp_ann <- gVars$pheno
      #}
      mysize = (nGenes* 20 ) + 10 * max(sapply(exp_ann, nchar)) + 10 * max(sapply(exp_ann, nchar))
      # print(paste("MY HEIGHT IS ---->", mysize))
      # print(paste("MY path IS ---->", nGenes))

      mysize = min(max(600, mysize),30e3)
      return(mysize)
    }else{
      return("auto")
    }
}

# pdf("heatmap.plot", width=width(nSamples, length_gene, nSamplesGene), height=height(nGenes, exp_ann))
# grid::grid.draw(toPlotGenes)
# dev.off()

return(toPlotGenes)

}
