#collapse a matrix of pathways X treatments based on the provided hierarchy and the deired level
#return a list of (matrix, subHierarchy) one element for each distinct collapsing level
# the value is computed summarising row elements with the provided function col_fun (default is median)
collapse_paths <- function (kegg_hierarchy,kegg_mat_cell, collapse_level=1,col_fun=function(x){median(x,na.rm = T)}) {
	#save(kegg_hierarchy,kegg_mat_cell,collapse_level,col_fun, file="demo/collapse_path.RData")

	#select pathways that are present in the matrix
	hierarchy_sub <- kegg_hierarchy[kegg_hierarchy$ID %in% colnames(kegg_mat_cell),]
	#order the matrix columns following the provided hyerarchy
	mat_cell_ord <- kegg_mat_cell[,as.character(hierarchy_sub$ID)]

	#get the list of categories on which we will collapse our pathways
	#keg_lev <- unique(hierarchy_sub[,collapse_level])
	if(collapse_level==1){
		keg_lev <- unique(hierarchy_sub[,collapse_level])
	}else{
		if(collapse_level==3){
			keg_lev <- hierarchy_sub[,collapse_level]
		}
		else{
			keg_lev <- unique(hierarchy_sub[,1:collapse_level])
			keg_lev = keg_lev[,2]
		}
	}


	#create a new matrix with columns corresponding to collapse categories
	if(is.null(nrow(mat_cell_ord))){
		mat_cell_ord = matrix(mat_cell_ord,ncol=1,nrow=length(mat_cell_ord),dimnames = list(as.character(names(mat_cell_ord)),hierarchy_sub$ID))
	}
	lev_mat <- matrix(NA,ncol = length(keg_lev),nrow = nrow(mat_cell_ord),dimnames = list(rownames(mat_cell_ord),keg_lev))

	#for each category take all pathways in that category and summaryze their values for each row (sample)
	for (lev in keg_lev){ #controllare che il padre e' arricchito. se e' arricchito ci metti il padre (ci metti asterisco), altrimenti fai il summary
		#take sub-hierarchy rooted in lev
		path_lev <- unique(hierarchy_sub[hierarchy_sub[,collapse_level] == lev, ]$ID)
		#extract columns with pathways in the defined sub-hierarchy
		mat_cell_sub <- as.matrix(mat_cell_ord[,colnames(mat_cell_ord) %in% path_lev])
		#summarize vaules over rows using the provided col_fun function
		summary <- unlist(apply(mat_cell_sub,1,FUN = col_fun))
		#store summarized column for culumn lev
		## if duplicate column names in mat it was assigning only the first
		# lev_mat[,lev] <- summary

		lev_mat[,colnames(lev_mat) %in% lev] <- summary
	}
	#create a nemed list element containing collapsed matrix and collapsed hierarchy as well
	lev_hierarchy <- unique(as.data.frame(hierarchy_sub[,1:collapse_level]))
	colnames(lev_hierarchy) <- colnames(hierarchy_sub)[1:collapse_level]

	return(list(lev_mat,lev_hierarchy))
}
