#TODO: refactor this shitty naming
kegg_mat_p <- function(EnrichDatList,hierarchy) {
  # ##cat("Inside kegg_mat_p")
  conmp_names = names(EnrichDatList)
  # #cat("Inside kegg_mat_p GO")
  names_un = unique(hierarchy$ID)
  kegg_mat_cell <- matrix(data = NA, nrow = length(conmp_names),ncol = length(names_un), dimnames = list(conmp_names,names_un))

  for (i in 1:length(EnrichDatList)){
    kegg = EnrichDatList[[i]]
    if(!gtools::invalid(kegg)){
      if (nrow(kegg)>0)
        for(kg in 1:nrow(kegg)){
          if(kegg$annID[kg] %in% colnames(kegg_mat_cell) ){
            kegg_mat_cell[conmp_names[[i]],kegg$annID[kg]] <- kegg$pValueAdj[kg]
          }else{
            1==1
            ##print(paste(kegg$annID[kg],"not found in the provided kegg hierarchy"))
          }
        }
    }
  }


  not_empty <- colSums(!is.na(kegg_mat_cell)) >0
  #print(not_empty)
  if (sum(not_empty) == 1){
    # #print("before")
    #print(kegg_mat_cell)

    col_name_un <- colnames(kegg_mat_cell)[not_empty]
    row_nams <- rownames(kegg_mat_cell)
    kegg_mat_cell <- kegg_mat_cell[,colSums(!is.na(kegg_mat_cell)) >0]
    kegg_mat_cell <- matrix(kegg_mat_cell, ncol=1)
    colnames(kegg_mat_cell) <- col_name_un
    rownames(kegg_mat_cell) <- row_nams
    # #print("after")
    #print(class(kegg_mat_cell))
  } else {
      kegg_mat_cell <- kegg_mat_cell[,colSums(!is.na(kegg_mat_cell)) >0]
  }


  #kegg_mat_cell <- kegg_mat_cell[,colSums(!is.na(kegg_mat_cell)) < nrow(kegg_mat_cell)]
  return(kegg_mat_cell)
}

kegg_mat_fc <- function(EnrichDatList,hierarchy,GList, summ_fun=median) {
  ##cat("Inside kegg_mat_fc")
  conmp_names = names(EnrichDatList)
  names_un = unique(hierarchy$ID)
  kegg_mat_cell <- matrix(data = NA, nrow = length(conmp_names),ncol = length(names_un), dimnames = list(conmp_names,names_un))

  ###print(head(kegg_mat_cell))

  for (i in 1:length(EnrichDatList)){
    kegg = EnrichDatList[[i]]
    if(!gtools::invalid(kegg)){
      if (nrow(kegg)>0)
        for(kg in 1:nrow(kegg)){
          if(kegg$annID[kg] %in% colnames(kegg_mat_cell) ){
            genes_in_path <- unlist(strsplit(as.character(kegg$gID[kg]),","))
            MM = GList[[i]]
            summFC <- summ_fun(as.numeric(MM[tolower(MM[,1]) %in% tolower(genes_in_path),2]))
            kegg_mat_cell[conmp_names[[i]],kegg$annID[kg]] <- summFC
          }else{
            ##print(paste(kegg$annID[kg],"not found in the provided kegg hierarchy"))
          }
        }
    }
  }

  not_empty <- colSums(!is.na(kegg_mat_cell)) >0
  #print(not_empty)
  if (sum(not_empty) == 1){
    #print("before")
    #print(kegg_mat_cell)

    col_name_un <- colnames(kegg_mat_cell)[not_empty]
    row_nams <- rownames(kegg_mat_cell)
    kegg_mat_cell <- kegg_mat_cell[,colSums(!is.na(kegg_mat_cell)) >0]
    kegg_mat_cell <- matrix(kegg_mat_cell, ncol=1)
    colnames(kegg_mat_cell) <- col_name_un
    rownames(kegg_mat_cell) <- row_nams
    #print("after")
    #print(class(kegg_mat_cell))
  } else {
      kegg_mat_cell <- kegg_mat_cell[,colSums(!is.na(kegg_mat_cell)) >0]
  }
  #kegg_mat_cell <- kegg_mat_cell[,colSums(kegg_mat_cell,na.rm = T)>0]
  return(kegg_mat_cell)
}

#kegg_mat_genes <- function(EnrichDatList, hierarchy, GList, summ_fun=median) {
kegg_mat_genes <- function(EnrichDatList, hierarchy) {
  #print("Inside kegg_mat_genes")
  conmp_names = names(EnrichDatList)
  #print("length(conmp_names)")
  #print(length(conmp_names))
  #print("conmp_names")
  #print(conmp_names)
  names_un = unique(hierarchy$ID)
  #print("length(names_un)")
  #print(length(names_un))
  #print("names_un")
  #print(names_un)
  kegg_mat_cell <- matrix(data=list(), nrow=length(conmp_names), ncol=length(names_un), dimnames=list(conmp_names,names_un))
  #kegg_mat_cell <- as.list(rep(NA, length(conmp_names)*length(names_un)))
  #dim(kegg_mat_cell) <- c(length(conmp_names), length(names_un))
  #colnames(kegg_mat_cell) <- names_un
  #rownames(kegg_mat_cell) <- conmp_names
  #print("dim(kegg_mat_cell)")
  #print(dim(kegg_mat_cell))

  ###print(head(kegg_mat_cell))
  for (i in 1:length(EnrichDatList)){
    kegg = EnrichDatList[[i]]
    if(!gtools::invalid(kegg)){
      if (nrow(kegg)>0)
        for(kg in 1:nrow(kegg)){
          if(kegg$annID[kg] %in% colnames(kegg_mat_cell)){
            genes_in_path <- unlist(strsplit(as.character(kegg$gID[kg]),","))
            ##print("str(genes_in_path)")
            ##print(str(genes_in_path))
            kegg_mat_cell[[conmp_names[[i]],kegg$annID[kg]]] <- genes_in_path
          }else{
            #print(paste(kegg$annID[kg],"not found in the provided kegg hierarchy"))
          }
        }
    }
  }

  #not_empty <- colSums(!is.na(kegg_mat_cell)) >0
  not_empty <- colSums(apply(kegg_mat_cell, c(1,2), function(x){!is.null(unlist(x))})) >0
  #print(not_empty)
  # if (sum(not_empty) == 1){
  #   #print("before")
  #   #print(kegg_mat_cell)
  #
  #   col_name_un <- colnames(kegg_mat_cell)[not_empty]
  #   row_nams <- rownames(kegg_mat_cell)
  #   #kegg_mat_cell <- kegg_mat_cell[,colSums(!is.na(kegg_mat_cell)) >0]
  #   kegg_mat_cell <- kegg_mat_cell[,not_empty]
  #   kegg_mat_cell <- matrix(kegg_mat_cell, ncol=1)
  #   colnames(kegg_mat_cell) <- col_name_un
  #   rownames(kegg_mat_cell) <- row_nams
  #   #print("after")
  #   #print(class(kegg_mat_cell))
  # } else {
  #     #kegg_mat_cell <- kegg_mat_cell[,colSums(!is.na(kegg_mat_cell)) >0]
  #     kegg_mat_cell <- kegg_mat_cell[,not_empty]
  # }

  ##Commented previous if else case added drop=FALSE to preserve matrix structure if only one column remains
  kegg_mat_cell <- kegg_mat_cell[,not_empty, drop=FALSE]

  #print("dim(kegg_mat_cell)")
  #print(dim(kegg_mat_cell))
  #kegg_mat_cell <- kegg_mat_cell[,colSums(kegg_mat_cell,na.rm = T)>0]
  return(kegg_mat_cell)
}
