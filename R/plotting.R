#' TODO input validation 
#' helper plotting function which is called in most of the clustering functions
#' function calls  collapse_paths, set_discrete
#' @import 
#' @param reduced_kegg_hierarchy
#' @param toPlotMap
#' @param level
#' @param clust_mat
#' @param pheno
#' @param samplesID
#' @param continuous , what to plot, user input
#' @param doGrouping , what to plot, user input
#' @param aspectRatio , plotting parameter
#' @return 
#' @export


plotting <- function(reduced_kegg_hierarchy,  toPlotMap, level, clust_mat, pheno, samplesID, continuous, doGrouping, aspectRatio) {

    kegg_nano_1 <- collapse_paths(kegg_hierarchy = reduced_kegg_hierarchy,kegg_mat_cell = toPlotMap, collapse_level = as.numeric(level))
    #extract collapsed matrix and collapsed hierarachy
    mat <- kegg_nano_1[[1]]
    hier <- kegg_nano_1[[2]]

    #TODO add validation
    #shiny::validate(need(expr = ncol(mat)>0, message = "No result for the enrichment or the filters are too restrictive. Please enlarge your selection"))

    print("Matrix dimension -->")
    print(dim(mat))
    print("Hierarchy dimension -->")
    print(dim(hier))
    #plot the collapsed matrix

    if(is.null(clust_mat)){
      exp_ann = pheno#cbind(c(rep(1:5,5),1),rownames(mat))
      if("All" %in% samplesID == FALSE){
        exp_ann = exp_ann[exp_ann[,2] %in% samplesID,]
      }
    }

    #I think thats for testing only
    #xxx = gVars$exp_ann
    #save(mat, hier,xxx,file="demo/demo.RData")

    ############################   DISCRETIZE MAT if user chose "discrete and there are negative and positve values"
    ############################
    mat_to_Plot=mat

    isDiscrete = set_discrete(continuous, mat_to_Plot)

    ########################################################
    #plot_grid is build in R function
    if(doGrouping){
      print("grouping selected")
      toPlot = plot_grid(path_mat = mat_to_Plot,path_hier = hier,experiment_ann = exp_ann,discrete =  isDiscrete,level_col = max(1,as.numeric(level)-1),square_colors=c(),color_leg=c(),path_text_size = 12,treat_text_size = 12, asRatio=(aspectRatio))
    }else{
      print("grouping NOT selected")
      level_n = max(1,as.numeric(level)-1)
      fake_hier = hier
      fake_hier[,level_n] = rep("",nrow(fake_hier))

      toPlot = plot_grid(path_mat = mat_to_Plot,path_hier = fake_hier,experiment_ann = exp_ann ,discrete =  isDiscrete,level_col =level_n,square_colors=c(),color_leg=c(),path_text_size = 12,treat_text_size = 12, asRatio=(aspectRatio))

    }

    return(toPlot)





}