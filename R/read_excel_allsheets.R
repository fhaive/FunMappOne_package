#' read excel file as a list of dataframe
#'
#' @import readxl
#' @param filename is the path to the file
#' @param tibble boolean specifying if the content of each sheet should be read as tibble or dataframe
#' @return a list with a position for each sheet in the excel file
#' @export

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X){
    y = readxl::read_excel(filename, sheet = X)
    if(nrow(y)==0){
      stop("empty sheet")
    }
    y[,1] = as.character(as.vector(y[[1]]))
    y
  })

  nExperiments = length(x)-1

  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets

  groupTab = x[[length(x)]]
  #check that the excel sheets are not empty
  for(i in 1:length(x)){
    if(nrow(x[[i]])==0){
      stop("empty sheet")
    }
  }

  # check that last sheet contains a dataframe with two columns
  if(ncol(groupTab)!=2 || nrow(groupTab)!=nExperiments){
    stop("The last sheet of the excel file should contain a table with 2 columns and the same number of rows of the other sheets in the file")
  }

  #check that all the experiments ID are in the group table
  if(sum(names(x)[1:(length(x)-1)] %in% groupTab[,1]) != (length(x)-1)){
    stop("all the experiment sheet names have to be reported in the group table")
  }

  # check that the sheet names do not contains the keywork "All"
  if(any(grepl("all", names(x),ignore.case = TRUE))){
    stop("the string all cannot be used in the sheet names")
  }

  #check that the grouping column contains numbers and that has to be sequentials
  if(class(groupTab[,2])!="numeric"){
    stop("the second column of the group tab has to contain numbers")
  }

  if(sum(1:max(groupTab[,2]) %in% groupTab[,2]) != max(groupTab[,2])){
    stop("the numbers in the second column of the group tab should be consecutive")
  }

  x = x[1:length(x)-1]
  groupTab = cbind(groupTab[,2], groupTab[,1])
  groupTab = as.data.frame(groupTab)
  colnames(groupTab) = c(colnames(groupTab)[2],colnames(groupTab)[1])

  return(list(GList = x, pheno = groupTab))

}
