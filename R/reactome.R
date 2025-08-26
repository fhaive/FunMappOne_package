load_mouse_reactome_hierarchy <- function() {
	data("mm_reactome_hierarchy")
	org = "mmu"
	reactome_hierarchy = mm_reactome_hierarchy
	reactome_hierarchy$ID = unlist(reactome_hierarchy$Pathway)
	reactome_hierarchy[,1] = mouse_map[reactome_hierarchy[,1],2]
	reactome_hierarchy[,2] = mouse_map[unlist(reactome_hierarchy[,2]),2]
	reactome_hierarchy[,3] = mouse_map[unlist(reactome_hierarchy[,3]),2]
	return(unique(reactome_hierarchy))
}

load_human_reactome_hierarchy <- function() {
	data("hm_reactome_hierarchy")
	org = "hsa"
	reactome_hierarchy = hm_reactome_hierarchy
	reactome_hierarchy$ID = unlist(reactome_hierarchy$Pathway)
	reactome_hierarchy[,1] = human_map[reactome_hierarchy[,1],2]
	reactome_hierarchy[,2] = human_map[unlist(reactome_hierarchy[,2]),2]
	reactome_hierarchy[,3] = human_map[unlist(reactome_hierarchy[,3]),2]
	return(unique(reactome_hierarchy))
}

load_rat_reactome_hierarchy <- function() {
	data("rat_reactome_hierarchy")
	org = "Rat"
	reactome_hierarchy = rat_reactome_hierarchy
	reactome_hierarchy$ID = unlist(reactome_hierarchy$Pathway)
	reactome_hierarchy[,1] = rat_map[reactome_hierarchy[,1],2]
	reactome_hierarchy[,2] = rat_map[unlist(reactome_hierarchy[,2]),2]
	reactome_hierarchy[,3] = rat_map[unlist(reactome_hierarchy[,3]),2]
	return(unique(reactome_hierarchy))
}
