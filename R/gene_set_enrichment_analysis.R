#' Compute uncorrected p-values from Bonferroni-corrected p-values by
#' multiplying each p-value by the total number of p-values
#'
#' @param pvalues The array of p-values to un-correct
#' @return The uncorrected p-values
#' @keywords internal
undo_bonferroni_correction <- function(pvalues) {
	unadjusted_pvalues <- pvalues / length(pvalues)
	return(unadjusted_pvalues)
}

#' Perform enrichment analysis for a collection of gene sets
#'
#' @param gene_sets A list of many gene symbols lists
#' @param organism The name of the organism
#' @param pathway_database The pathway database to use for the enrichment
#' @param only_annotated_genes Whether to consider only genes known to belong to a pathway
#' @param pvalue_correction Correction of the enrichment p-values
#' @param pvalue_threshold Significance threshold
#' @param return_only_significant Whether to return only significant enriched genes
#' @param min_intersection_size Smallest amount of genes matching a pathway
#' @param output_map_type Type of output values in the heatmap returned
#' @param aggregation_function
#' @return hierarchy The whole pathway hierarchy
#' @return lvl1_h The first level of the collapsed hierarchy
#' @return lvl2_h The second level of the collapsed hierarchy
#' @return lvl3_h The third level of the collapsed hierarchy
#' @return KEGG_MAT
#' @return KEGG_MAT_GENES
#' @return enriched_data_list Enriched genes lists
#' @export
gene_set_enrichment_analysis <- function(
	gene_sets, organism="mmusculus", pathway_database="KEGG", only_annotated_genes=F,
	pvalue_correction="fdr", pvalue_threshold=0.05, return_only_significant=T,
	min_intersection_size=0, output_map_type="FCPV",aggregation_function = "mean") {
	#----- CORE ENRICHMENT ANALYSIS
	if(pvalue_correction == "none") {
		enriched_gene_sets = lapply(
			gene_sets,
			enrich,
			pathway_database,
			organism,
			1.0,
			"bonferroni",
			FALSE,
			min_intersection_size,
			only_annotated_genes
		)
		# TODO: probably it's better to move this into enrich function
		# TODO: filter significant genes only when return_only_significant is TRUE
		for(i in 1:length(enriched_gene_sets)) {
			egs = enriched_gene_sets[[i]]
			egs$pValueAdj = undo_bonferroni_correction(egs$pValueAdj)
			egs$pValue = undo_bonferroni_correction(egs$pValue)
			significant_genes = egs$pValueAdj <= pvalue_threshold
			enriched_gene_sets[[i]] = egs[significant_genes, ]
		}
	} else {
		enriched_gene_sets = lapply(
			gene_sets,
			enrich,
			pathway_database,
			organism,
			pvalue_threshold,
			pvalue_correction,
			return_only_significant,
			min_intersection_size,
			only_annotated_genes
		)

	}

	#-----
	# TODO: move this into it's own function
	#----- HIERARCHY CONSTRUCTION
	if(pathway_database == "REAC") {
		if(organism == "hsapiens") {
			hierarchy = load_human_reactome_hierarchy()
		} else if (organism == "mmusculus") {
			hierarchy = load_mouse_reactome_hierarchy()
		} else if (organism == "rnorvegicus") {
			hierarchy = load_rat_reactome_hierarchy()
		} else {
			stop(paste0(
				"Invalid value for parameter organism, ",
				"got {organism}, expected \"hsapiens\", ",
				"\"mmusculus\" or \"rnorvegicus\""
			))
		}
	} else if (pathway_database == "KEGG") {
		hierarchy = load_kegg_hierarchy()
	} else if (pathway_database == "GO:BP") {
		makeGOGraph(ont = "bp") -> geograph
		igraph.from.graphNEL(geograph) -> igraphgeo
		igraphgeo = as.undirected(igraphgeo)
		root = "GO:0008150"
	} else if (pathway_database == "GO:CC") {
		makeGOGraph(ont = "cc") -> geograph
		igraph.from.graphNEL(geograph) -> igraphgeo
		igraphgeo = as.undirected(igraphgeo)
		root="GO:0005575"
	} else if (pathway_database == "GO:MF") {
		#EnrichDatList = all_GO_MF
		makeGOGraph(ont = "mf") -> geograph
		igraph.from.graphNEL(geograph) -> igraphgeo
		igraphgeo = as.undirected(igraphgeo)
		root="GO:0003674"
	} else {
		error(paste0(
			"Invalid value for parameter pathway_database, ",
			"got {pathway_database}, expected \"REAC\", ",
			"\"KEGG\", \"GO:BP\", \"GO:CC\" or \"GO:MF\""
		 ))

	}
	if (pathway_database %in% c("GO:BP", "GO:CC", "GO:MF")) {
		go_type = substr(pathway_database, 4, 5)
		#find the list of term into  the enriched list
		res2 = filterGO(enriched_gene_sets, go_type = go_type)
		enriched_gene_sets = res2$EnrichDatList
		go_terms = res2$goTerm
		#print("compute GO hierarchy")

		#Compute all the shortest path between the root to the go_terms
		asp = all_shortest_paths(graph = igraphgeo,from = root, to = go_terms)

		#reduce the shortest path to length 3 to build a 3 level hierarchy
		go_hierarchy = matrix("",nrow = length(asp$res),ncol = 3)

		for(i in 1:length(asp$res)){
			nn = names(asp$res[[i]])
			nn = nn[2:length(nn)]
			if(length(nn)<3){
				nn2 = rep(nn[length(nn)],3)
				nn2[1:length(nn)] = nn
			}else{
				nn2 = c(nn[1:2],nn[length(nn)])
			}
			go_hierarchy[i,] =nn2
		}

		go_hierarchy = unique(go_hierarchy)
		colnames(go_hierarchy) = c("level1","level2","level3")
		go_hierarchy = as.data.frame(go_hierarchy)
		go_hierarchy$ID = go_hierarchy[,3]

		idx = which(is.na(go_hierarchy[,1]))
		if(length(idx)>0){
			go_hierarchy = go_hierarchy[-idx,]
		}
		go_hierarchy[,1] = Term(GOTERM[as.character(go_hierarchy[,1])])
		go_hierarchy[,2] = Term(GOTERM[as.character(go_hierarchy[,2])])
		go_hierarchy[,3] = Term(GOTERM[as.character(go_hierarchy[,3])])

		# colnames(hier_names) = c("Level1","Level2","Pathway")
		hierarchy = go_hierarchy
		#print(gVars$hierarchy)
	}
	#-----
	#----- COMPUTE SCORE MATRIX FOR HEATMAP
	NCol = ncol(gene_sets[[1]])
	if(output_map_type == "PVAL"){
		#print("only genes")
		M = kegg_mat_p(enriched_gene_sets,hierarchy = hierarchy)
		if(length(enriched_gene_sets)==1){
		  MM = matrix(M, nrow = 1, ncol = length(M))
		  colnames(MM) = names(M)
		  rownames(MM) = names(enriched_gene_sets)
		  M = as.data.frame(MM)
		}
	}else{

		if(NCol==1){
			M = kegg_mat_p(enriched_gene_sets,hierarchy = hierarchy)
      if(length(enriched_gene_sets)==1){
        MM = matrix(M, nrow = 1, ncol = length(M))
        colnames(MM) = names(M)
        rownames(MM) = names(enriched_gene_sets)
        M = as.data.frame(MM)
      }
		}else{
			M1 = kegg_mat_p(enriched_gene_sets,hierarchy = hierarchy)
			M2 = kegg_mat_fc(
				EnrichDatList=enriched_gene_sets,
				hierarchy=hierarchy,
				GList=gene_sets,
				summ_fun=get(aggregation_function)
			)

			# if(input$MapValueType == "PVAL"){
			#   M = M1
			# }
			if(output_map_type == "FC"){
				M = M2
			}
			if(output_map_type == "FCPV"){
				M = M2 * -log(M1)
			}

			if(length(enriched_gene_sets)>1){
			  rownames(M) = rownames(M1)
			  colnames(M) = colnames(M1)
			}else{
			  MM = matrix(M, nrow = 1, ncol = length(M))
			  colnames(MM) = names(M)
			  rownames(MM) = names(enriched_gene_sets)
			  M = as.data.frame(MM)
			}

		}
	}

	KEGG_MAT = M
	#------
	#------ HIERARCHY PLOT
	#TODO: heavy rename!!!
	local_hierarchy <- collapse_paths(
		kegg_hierarchy=hierarchy,
		kegg_mat_cell=KEGG_MAT,
		collapse_level=3
	)
	mat <- local_hierarchy[[1]]
	hier <- local_hierarchy[[2]]

	lev1_h = as.list(c("All",unique(as.character(hier[,1]))))

	lev2_h = list("All" = "All")
	for(i in unique(hier[,1])){
		lev2_h[[i]] = as.character(unique(hier[which(hier[,1] %in% i),2]))
	}

	lev3_h = list("All" = "All")
	for(i in unique(hier[,2])){
		lev3_h[[i]] = as.character(unique(hier[which(hier[,2] %in% i),3]))
	}

	#Compute genes matrix
	Mgenes <- kegg_mat_genes(EnrichDatList=enriched_gene_sets, hierarchy=hierarchy)
	#------

	return(list(
		hierarchy=hierarchy,
		lvl1_h=lev1_h,
		lvl2_h=lev2_h,
		lvl3_h=lev3_h,
		KEGG_MAT=KEGG_MAT,
		KEGG_MAT_GENES=Mgenes,
		enriched_data_list = enriched_gene_sets
	))
}
