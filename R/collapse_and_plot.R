# TODO Check if input variables can be improved
# TODO Test
#' Collapse and plot
#' Collapses paths and returns a plot object.
#' @param reduced_kegg_hierarchy The kegg hierarchy that will be used to collapse paths.
#' @param toPlotMap The map to be plotted.
#' @param experiment_ann Experiment annotation.
#' @param samplesID Row names for samples to be plotted. "All" to plot all of them.
#' @param collapse_level Must be between 1 and 3.
#' @param doGrouping Logical
#' @param aspectRatio Logical
#' @param continuous Logical. If FALSE values are discretized with negative, 0, and positive mapped to -1, 0, and 1.
#' @return Object to plot.
#' @keywords internal
collapse_and_plot <- function(
  reduced_kegg_hierarchy,
  toPlotMap,
  experiment_ann,
  samplesID = "All",
  collapse_level = 1,
  doGrouping = TRUE,
  aspectRatio = TRUE,
  continuous = FALSE) {

  if (is.null(reduced_kegg_hierarchy)) {
    stop("reduced_kegg_hierarchy is null.")
  }
  
  if (is.null(toPlotMap)) {
    stop("toPlotMap is null.")
  }
  
  if (is.null(experiment_ann)) {
    stop("experiment_ann is null.")
  }
  
  if (!content_selection_ok(samplesID)) {
    stop("Error in samplesID.")
  }
  
  if (!(collapse_level >= 1 && collapse_level <= 3)) {
    stop("collapse_level not between 1 and 3.")
  }
  
  if (!is.logical(aspectRatio)) {
    stop("aspectRatio not logical.")
  }
  
  if (!is.logical(continuous)) {
    stop("continuous not logical.")
  }

  if (! ("All" %in% samplesID)) {
    toPlotMap = toPlotMap[rownames(toPlotMap) %in% samplesID,]
  }

  kegg_nano_1 <- collapse_paths(kegg_hierarchy = reduced_kegg_hierarchy,
                                kegg_mat_cell = toPlotMap, collapse_level = collapse_level)
  #extract collapsed matrix and collapsed hierarachy
  mat <- kegg_nano_1[[1]]
  hier <- kegg_nano_1[[2]]

  if (ncol(mat) == 0) {
    stop("No result for the enrichment or the filters are too restrictive. Please enlarge your selection")
  }

  # Discretize mat if user choses so.
  if (!continuous){
    mat[mat<0]=-1
    mat[mat>0]=1
  }

  if(doGrouping){

    toPlot = plot_grid(
      path_mat = mat,
      path_hier = hier,
      experiment_ann =  experiment_ann,
      discrete = !continuous,
      level_col = max(1,as.numeric(collapse_level)-1),
      square_colors=c(),
      color_leg=c(),
      path_text_size = 12,
      treat_text_size = 12,
      asRatio = aspectRatio)

  }else{

    level_n = max(1,as.numeric(collapse_level)-1)
    fake_hier = hier
    fake_hier[,level_n] = rep("",nrow(fake_hier))

    toPlot = plot_grid(
      path_mat = mat,
      path_hier = fake_hier,
      experiment_ann = experiment_ann,
      discrete = !continuous,
      level_col = level_n,
      square_colors=c(),
      color_leg=c(),
      path_text_size = 12,
      treat_text_size = 12,
      asRatio = aspectRatio)
  }

  return (toPlot)

}
