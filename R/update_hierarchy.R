#' Checks validity of level content parameter.
#' @param lev_content Level content parameter to check.
#' @return True if ok.
#' @keywords internal
content_selection_ok <- function(lev_content) {
  if (is.null(lev_content)) {
    return (FALSE)
  }
  if (!is.character(lev_content)) {
    return (FALSE)
  }
  return (TRUE)
}


#' Updates the hierarchy.
#' Updates the hierarchy based on the specified content for the three levels.
#' @param kegg_hierarchy The hierarchy to update.
#' @param lev1_content Content selection for first level.
#' @param lev2_content Content selection for second level.
#' @param lev3_content Content selection for third level.
#' @return The hierarchy restricted to the columns that are specified by the content parameters.
#' @keywords internal
update_hierarchy=function(
  kegg_hierarchy,
  lev1_content = "All",
  lev2_content = "All",
  lev3_content = "All"){
  
  if (is.null(kegg_hierarchy)) {
    stop("kegg_hierarchy is NULL.")
  }
  
  if (!(content_selection_ok(lev1_content))) {
    stop("Error in lev1_content.")
  }
  if (!(content_selection_ok(lev2_content))) {
    stop("Error in lev2_content.")
  }
  if (!(content_selection_ok(lev3_content))) {
    stop("Error in lev3_content.")
  }
  
  # Updating the hierarchy
  if("All" %in% lev1_content){
    # All in lev1
    idx1 = 1:nrow(kegg_hierarchy)
  }else{
    idx1= which(kegg_hierarchy[,1] %in% lev1_content)
  }
  if("All" %in% lev2_content){
    # All in lev2
    idx2 = 1:nrow(kegg_hierarchy)
  }else{
    idx2= which(kegg_hierarchy[,2] %in% lev2_content)
  }
  if("All" %in% lev3_content){
    # All in lev3
    idx3 = 1:nrow(kegg_hierarchy)
  }else{
    idx3= which(kegg_hierarchy[,3] %in% lev3_content)
  }
  
  idx = intersect(idx1,intersect(idx2,idx3))
  
  return(kegg_hierarchy[idx,])
  
}