#' store the list of dataframe with the enriched pathways in an excel file
#'
#' @import readxl
#' @param filePath filename of the xlsx file
#' @param enriched_data_list list of dataframes given in output by the gene_set_enrichment_analysis funciton
#' @export

write_pathways_to_excel = function(filePath, enriched_data_list){

  if(file.exists(filePath)){
    stop("File already exist, please use a different file name.")
  }

  for (i in 1:length(enriched_data_list)){
   # print(i)
    if(nrow(enriched_data_list[[i]])>0){
      write.xlsx(x=enriched_data_list[[i]],file=filePath,sheetName=names(enriched_data_list)[i],append=T)
      XLConnect::xlcFreeMemory()
    }

  }
}
