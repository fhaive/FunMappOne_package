save_map_plot <- function(filename, plot_obj) {
	if(!("ggplot" %in% class(plot_obj))) {
		stop(
			paste0(
				"Unexpected class for the plot object: ", 
				paste0(class(p), collapse=" ")
			)
		) 
	}
	if(!(grepl("pdf$", filename))) {
		filename <- paste0(filename, ".pdf")
	}
	pdf(filename)
	print(plot_obj)
	dev.off()
}
