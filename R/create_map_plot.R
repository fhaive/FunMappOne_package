# source("R/update_hierarchy.R")

# TODO Check if input variables can be improved
# TODO Test
#' @param lev1_content Level 1 content filter. Names of columns or "All".
#' @param lev2_content Level 2 content filter. Names of columns or "All".
#' @param lev3_content Level 3 content filter. Names of columns or "All".
#' @param kegg_hierarchy Kegg hierarchy that will be reduced by filters and plotted.
#' @param toPlotMap Map to be plotted.
#' @param experiment_ann Experiment annotations.
#' @param samplesID Samples to be plotted. Names of rows or "All".
#' @param collapse_level Level for collapsing content. Between 1 and 3.
#' @param doGrouping To group or not. if TRUE the original grouping from the input is used.
#' @param aspectRatio Logical.
#' @param continuous Logical. If FALSE values are discretized with negative, 0, and positive mapped to -1, 0, and 1.
#' @return Object to plot.
#' @export
create_map_plot <- function(
  lev1_content,
  lev2_content,
  lev3_content,
  kegg_hierarchy,
  toPlotMap,
  experiment_ann,
  samplesID,
  collapse_level,
  doGrouping,
  aspectRatio,
  continuous) {
  
  reduced_kegg_hierarchy =
    update_hierarchy(kegg_hierarchy = kegg_hierarchy ,lev1_content, lev2_content, lev3_content)

  toPlot <- collapse_and_plot(
    reduced_kegg_hierarchy,
    toPlotMap,
    experiment_ann,
    samplesID,
    collapse_level,
    doGrouping,
    aspectRatio,
    continuous)

  return (toPlot)

}
