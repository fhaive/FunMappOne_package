#' convert gene identifiers
#'
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @import AnnotationDbi
#' @param GList is a list containing dataframes whose first column is the gene identifiers
#' @param organism is a string specifying the organism type. Possible values: Human, Mouse, Rat. Default value = Human
#' @param annType is a string specifying the type of gene identifier. Possible values: SYMBOL, ENSEMBL, ENTREZID. Default value = SYMBOL
#' @return the original GList with gene identifier converted#'
#' @export
#'
convert_genes = function(GList, organism = "Human", annType = "SYMBOL"){

  if((organism %in% c("Human", "Mouse","Rat")==FALSE)){
    stop("Organism value must be one between Human, Mouse or Rat!")
  }

  if((annType %in% c("SYMBOL", "ENSEMBL","ENTREZID")==FALSE)){
    stop("annType value must be one between SYMBOL, ENSEMBL or ENTREZID!")
  }

  orgLibs <- list("Human"=org.Hs.eg.db, "Mouse"=org.Mm.eg.db, "Rat" = org.Rn.eg.db)
  orgDB <- orgLibs[[organism]]

  if(annType == "SYMBOL"){
    tmp = GList[[1]]

    tryCatch(
      {
        tmp <- AnnotationDbi::select(orgDB, keys=as.character(tmp[,1]), columns="SYMBOL", keytype=annType)
        return(GList)
      },error=function(cond){
        stop("the genes do not maps the organism or you are specifying the wrong annotation type")
      }
    )
  }

  for(i in 1:length(GList)){
    M=GList[[i]]
    genes = M[,1]

    tryCatch(
      {
        selectAnnDF <- AnnotationDbi::select(orgDB, keys=genes, columns="SYMBOL", keytype=annType)

      },error=function(cond){
        stop("the genes do not maps the organism or you are specifying the wrong annotation type")
      }
    )


    toRem = which(is.na(selectAnnDF[,2]))
    if(length(toRem)>0){
      selectAnnDF = selectAnnDF[-toRem,]
    }

    #remove duplicates, if any
    selectAnnDF = selectAnnDF[!duplicated(selectAnnDF$SYMBOL),]

    M = M[M[,1] %in% selectAnnDF[,1],]
    M = as.matrix(M)

    # macke sure they match by the key type
    # matches = match(selectAnnDF[,1],M[,1])
    # 
    # M[matches,1] = selectAnnDF[,2]
    # rownames(M) = M[,1]
    #update gene symbols
    GList[[i]] = M

  }

  return(GList)
}
