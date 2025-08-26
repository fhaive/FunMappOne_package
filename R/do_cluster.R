#' TODO input validation
#' TODO update description
#' function calls range01
#' @import
#' @param MKEGG_MAT, KEGG_MAT
#' @param Distance, distance measure to be used, user input
#' @param ClusterMethod, clustering method to be used, user input, from R function values are "ward.D, Ward.D2, "single, complete, average, mcquitty,median, centroid
#' @param nc, user input #TODO check what this actually is/ does
#' @param samplesID
#' @return list
#' @export

do_cluster <- function (KEGG_MAT, samplesID = c("All"),nc=1, Distance="euclidean", ClusterMethod="ward.D") {
    #indicate to front end to show loading

    #TODO add 3rd option
    #create matrix
    if(!isTRUE((Distance %in% "euclidean") || (Distance %in% "jaccard") || (Distance %in% "comb"))){

      stop("Unknown Distance value")
      return(NULL)
    }

    M=KEGG_MAT

    D = matrix(data = 0,nrow = nrow(M),ncol = nrow(M))
    rownames(D) = colnames(D) = rownames(M)
    for(i in 1:(nrow(M)-1)){
      pi= colnames(M)[!is.na(M[i,])]
      for(j in (i+1):nrow(M)){
        pj= colnames(M)[!is.na(M[j,])]
        D[i,j] = D[j,i] = length(intersect(pi,pj))/length(union(pi,pj))
      }
    }
    idx= which(rowSums(!is.na(M)) == 0)
    #print(idx)

    if(length(idx)>0){
      D[idx,] = 0
      D[,idx] = 0
      D[idx,idx] = 1
    }
    diag(D) = 1

    D1 = matrix(data = 0,nrow = nrow(M),ncol = nrow(M))
    rownames(D1) = colnames(D1) = rownames(M)

    for(i in 1:(nrow(M)-1)){
      idx1 = which(!is.na(M[i,]))
      #print(idx1)
      for(j in (i+1):nrow(M)){
        idx2 = which(!is.na(M[j,]))
        #print(idx2)
        idx = intersect(idx1,idx2)
        #print(idx)

        if(length(idx)>0){
          D1[i,j] = D1[j,i] = as.numeric(dist(t(cbind(M[i,idx],M[j,idx]))))
        }

      }
    }

    D = 1-D

    #View(D)

    # print(class(D1))
    # print(dim(D1))


    D1 = range01(D1)
    # print(class(D1))
    # print(dim(D1))

    #View(D1)

    if(Distance %in% "euclidean"){
      DD = D1
    }else{
      if(Distance %in% "jaccard"){
        DD = D
      }else{
        DD = (D + D1)/2
      }
    }

    print(class(DD))
    print(dim(DD))
    #View(DD)

    hls = hclust(as.dist(DD),method = ClusterMethod)
    hls <- hls
    #plot(hls)

    #save plot ggplot

    plot(hls,xlab="", sub="",hang = -1)
    if(as.numeric(nc)>1){
      rect.hclust(tree = hls,k = as.numeric(nc))
    }



    #print(rownames(gVars$M))
    #print(hls$order)
    clust_mat = KEGG_MAT[hls$order,]


    cat("CHECK ROWNAMES -->")
    #print(rownames(clust_mat))
    nClust = as.numeric(nc)

    cls = cutree(hls,k=nClust)
    cls = cls[hls$order]
    #i guess that means we need to return cls
    #gVars$cls <- cls

    cat("CLUSTERING RESULTS -->")
    #print(cls)

    #return exp_ann
    exp_ann = cbind(cls,names(cls))#nrow(gVars$KEGG_MAT)

    if("All" %in% samplesID == FALSE){
      exp_ann = exp_ann[exp_ann[,2] %in% samplesID,]
    }

    #print(gVars$exp_ann)
    #toPlotMap = clust_mat


    #toPlot = plotting(reduced_kegg_hierarchy,  toPlotMap, level, clust_mat, pheno, samplesID, continuous, doGrouping, aspectRatio)

    #print("after PLOTS in clustering")
    #clustered <- clustered + 1


    return (list(clust_mat=clust_mat, cls=cls, hls=hls, exp_ann=exp_ann))


  }
