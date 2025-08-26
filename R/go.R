#TODO: refactor this mess
filterGO <- function(EnrichDatList,go_type="BP"){
  goTerm =c()
  for(i in 1:length(EnrichDatList)){
    df = EnrichDatList[[i]]
    # idx = which(df$ONTOLOGY %in% go_type)
    # if(length(idx)>0){
    #   EnrichDatList[[i]] = df[idx,]
    #   goTerm = union(goTerm,df[idx,"GOID"])
    # }else{
    #   EnrichDatList[[i]] = NULL
    # }
    goTerm = union(goTerm,df[,"annID"])

  }

  XX = AnnotationDbi::select(x = GO.db,columns = c("GOID","ONTOLOGY"),keys = goTerm)
  toRem = which(is.na(XX[,2]))
  if(length(toRem)>0){
    toRemGo = XX[toRem,1]
    for(i in 1:length(EnrichDatList)){
      df = EnrichDatList[[i]]
      toRem2 = which(df$annID %in% toRemGo)
      if(length(toRem2)>0){
        EnrichDatList[[i]] = EnrichDatList[[i]][-toRem2,]
      }
    }
  }

  if(length(toRem)>0){
    goTerm = XX[-toRem,1]

  }

  return(list(EnrichDatList=EnrichDatList,goTerm=goTerm))
}
