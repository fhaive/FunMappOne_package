 #' TODO input validation
 #' Need to add apropriate user feedback
#' function calls collapse paths, plot_grid_genes
#' @import dplyr, purrr, tibble
#' @param KEGG_MAT_GENES matrix returned from the enrichment function
#' @param selPath # the name of the pathway/term that the user can select. Eg. if levelGene=1 a possible value for selPath can be Metabolism if the enrichment has beed done for KEGG. Possible values are those lev1, lev2 and lev3 returned by the enriched pathways. we assume that the input is correct and handled by UI
#' @param levelGene, user input the level of the hierarchy. It can be 1,2,or 3.
#' @param kegg_hierarchy, hierarchy of kegg to be used, returned from enrichment function. It is the one already filtered by the pathways in the experiments
#' @param GList output of the function read excel. is a list containing dataframe in each posistion. the last one is the phedo table specifying how the experiment are grouped.
#' @param selScoreType, user input it is a string specifying which value for the genes we want to plot. Possible values list("logFC"="lfc", "P-Value"="pval", "Combined"="comb"). Default is lfc. These string must be unified with the one used in the frevious functions (e.g. enrichemnt function)
#' @param KEGG_MAT
#' @param pheno is the grouping (last dataframe of the excel file). It is the last object in list GList. we need to check if the other functions accept GList with this last element or not. in that case we need to split it and have in input two parameters. otherwise pheno is not required as parameter and can be used only GList
#' @param aspectRatio, user input
#' @param hls # i guess that needs to be returned by do_cluster & reset_cluster
#' @param cls # i guess that needs to be returned by do_cluster & reset_cluster
#' @return figure to be plotted /ggplot object
#' @export

do_gene_heat_map <- function(KEGG_MAT_GENES, kegg_hierarchy, GList, KEGG_MAT, pheno, hls, cls, aspectRatio, selPath = "Metabolism", levelGene=1,  selScoreType = "lfc", continuous = TRUE) {

    #if glist is dataframe with only one column you cannot run

    plotGeneMapCheck <- ifelse(any(purrr::map(GList, ncol)<2), FALSE, TRUE)
    # add user feedback
    if(is.null(plotGeneMapCheck) || !isTRUE(plotGeneMapCheck)){
      stop("Gene list contains too less columns")
      return(NULL)
    }

    #should com from gene_set_enrichment
    if(is.null(KEGG_MAT_GENES) || nrow(KEGG_MAT_GENES)==0){
      #shinyjs::info("Enrichment matrix is NULL or Empty!!!")
      stop("Enrichment matrix is NULL or Empty!!!")
      return(NULL)
    }

    #we assume term is correctly handled by GUI since there are too many options
    #has to be a pathway => str
    if(is.null(selPath) || selPath=="NONE"){
      stop("Please select a Pathway Term!")

      return(NULL)
    }

    if(!isTRUE((as.integer(levelGene) == 1) || (as.integer(levelGene) == 2) || (as.integer(levelGene) == 3))){
      stop("Level needs to be 1, 2 or 3, current value is not permitted")
      return(NULL)
    }

    if(!isTRUE((selScoreType == "lfc") || (selScoreType == "pval") || (selScoreType == "comb"))){
      stop("Please select a valid scoring value, options are logFC, p-value or combination")
      return(NULL)
    }

    selLevel <- as.integer(levelGene)
    #selTerm <- selPath
    #kegg_hierarchy <- hierarchy

    #Filter hierarchy by selected level and term
    kegg_hierarchy_sel <- dplyr::filter(kegg_hierarchy, kegg_hierarchy[,selLevel]==selPath)
    kegg_hierarchy_sel_pathIDs <- unique(as.character(kegg_hierarchy_sel$ID))

    pathIdMapIdx <- which(kegg_hierarchy_sel_pathIDs %in% colnames(KEGG_MAT_GENES))

    if(length(pathIdMapIdx)==0){

      stop("No enriched information for the selected TERM!")
      return(NULL)
    }else{
      kegg_hierarchy_sel_pathIDs <- kegg_hierarchy_sel_pathIDs[pathIdMapIdx]
    }

    KEGG_MAT_GENES_sel <- KEGG_MAT_GENES[, kegg_hierarchy_sel_pathIDs]

    #names of the gene enriched in the experiment
    KEGG_MAT_GENES_enriched_genes <- unique(unlist(KEGG_MAT_GENES_sel))

    #Gene tables uploaded by user
    #GList <- gVars$GList
    egGeneTable <- GList[[1]]

    #Score type selected by user
    #selScoreType <- input$selScoreType

    if(selScoreType %in% c("comb", "pval") && ncol(egGeneTable)<3){
      stop("Missing required scores! Options are comb or pval")
      return(NULL)
    }

    GListMod <- NULL
    if(selScoreType=="comb"){
      GListMod <- lapply(GList,
                         function(x){
                           x <- mutate(x, dscore=x[2] * -log10(x[3]))
                           x <- x[,c(1,4)]
                           return(x)
                         }
      )
    }else{
      valMapVec <- c("lfc"=2, "pval"=3)
      selColIDx <- unname(valMapVec[selScoreType])
      #print("selColIDx")
      #print(selColIDx)
      GListMod <- lapply(GList, function(x){
                                   x <- x[,c(1,selColIDx)]
                                   if(sum(class(x)=="character")==1){
                                     x1 = matrix(data = x, ncol = 2)
                                     colnames(x1) = names(x)
                                     return(x1)
                                   }
                                   return(x)
                                }
      )
    }

    if(is.null(GListMod)){
      stop("Unexpected response from score selection!")
      return(NULL)
    }

    ##Get the reduced data frame from combining all samples
    geneTableColNames <- colnames(egGeneTable)
    geneIDCol <- geneTableColNames[1]

    #Convert gene table matrices to data frames
    GListMod2 = lapply(GListMod, as.data.frame)

    reducedGeneTable <- purrr::reduce(GListMod2, dplyr::full_join, by=geneIDCol)

    reducedGeneTableSel <- dplyr::filter(reducedGeneTable, tolower(reducedGeneTable[,geneIDCol]) %in% tolower(KEGG_MAT_GENES_enriched_genes))

    reducedGeneTableSel <- tibble::column_to_rownames(reducedGeneTableSel, geneIDCol)

    colnames(reducedGeneTableSel) <- names(GList)


    if(ncol(reducedGeneTableSel)==0){
      stop("No result for the enrichment or the filters are too restrictive. Please enlarge your selection")
      return(NULL)
    }

    # matrix with genes on the column and experiments on the rows. Each cell contains the selected selScoreType value for each gene in that experiment
    mat_to_Plot <- t(as.matrix(reducedGeneTableSel))

    #KEGG_MAT <- gVars$KEGG_MAT
    #check if pathway is available
    pathIdMapIdx <- which(kegg_hierarchy_sel_pathIDs %in% colnames(KEGG_MAT))
    gridPlotTypeGoup <- NULL
    #TODO add user feedback
    if(length(pathIdMapIdx)==0){
      stop("Did not find TERM enrichment information in the pathway matrix!!!")
      return(NULL)
    }else{
      #print("############ Adding path information to gene plot ############")
      kegg_hierarchy_sel_pathIDs <- kegg_hierarchy_sel_pathIDs[pathIdMapIdx]
      KEGG_MAT_sel <- KEGG_MAT[, kegg_hierarchy_sel_pathIDs, drop=FALSE]

      res_collapsed <- collapse_paths(kegg_hierarchy=kegg_hierarchy_sel, kegg_mat_cell=KEGG_MAT_sel, collapse_level=as.numeric(selLevel))
      #extract collapsed matrix and collapsed hierarachy
      KEGG_MAT_sel_collaped <- res_collapsed[[1]]
      KEGG_MAT_sel_collaped_hier <- res_collapsed[[2]]

      #return(NULL)

      gridPlotTypeGoup <- factor(
        c(
          rep("Pathway", ncol(KEGG_MAT_sel_collaped)*nrow(mat_to_Plot)),
          rep("Genes", ncol(mat_to_Plot)*nrow(mat_to_Plot))
        ),
        levels=c("Pathway", "Genes")
      )
      mat_to_Plot <- cbind(KEGG_MAT_sel_collaped, mat_to_Plot)
    }



    ############################   DISCRETIZE MAT if user chose "discrete and there are negative and positve values"
    ############################


    #need to return ?
    nSamplesGene = nrow(mat_to_Plot)
    nGenes = ncol(mat_to_Plot)

    #i think nchar is a build in function??? yes it is
    #Get maximum character length for the columns
    MLG = max(unlist(lapply(X = colnames(mat_to_Plot),FUN = nchar)))
    #to return ?
    length_gene = MLG


    #Performing discretization
    #possible to put this into a subfunction
    isDiscrete = set_discrete(continuous, mat_to_Plot)

    if(!is.null(hls)){
      mat_to_Plot = mat_to_Plot[hls$order,]
    }

    if(!is.null(cls)){
      pheno = cbind(cls,names(cls))
    }

    toPlotGenes = plot_grid_genes(path_mat=mat_to_Plot, experiment_ann=pheno, gene_group=gridPlotTypeGoup, discrete=isDiscrete, square_colors=c(), color_leg=c(), treat_text_size=12, asRatio=(aspectRatio))
    # }

    #TODO if needed but sinced we are just stroing right now I guess we dont need it
    #check for display size, pop up message if too large
    #if(!is.null(gVars$nGenes) && ((gVars$nGenes* 20 ) + 10 * max(sapply(exp_ann, nchar))) > 30e3){
     # print("Exceeding dimensions")
     # shinyjs::info("Warning: too many functional categories, the map might not be readable. Download the PDF for better resolution image.")
    #}
    #print("Exiting doGeneHeatMap!")
    return(list(toPlotGenes = toPlotGenes, res_collapsed = res_collapsed, mat_to_Plot=mat_to_Plot))
  }

