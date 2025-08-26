library("FunMappOnePackage")

inputs = read_excel_allsheets("Mouse_Test_Data.xlsx", tibble = FALSE)
GList = inputs[[1]]
pheno = inputs[[2]]
converted_GList = convert_genes(GList, organism = "Mouse", annType = "SYMBOL")

enrichedList = gene_set_enrichment_analysis(gene_sets = GList,
                             organism="mmusculus",
                             pathway_database="KEGG",
                             only_annotated_genes=F,
                             pvalue_correction="fdr",
                             pvalue_threshold=0.05,
                             return_only_significant=F,
                             min_intersection_size = 0)

# make function that store the enriched list in a excel file. One table in each sheet. Name sheet with list name
library(xlsx)

file_for_excel_write_test <- "temp.xlsx"
if (file.exists(file_for_excel_write_test)) {
  file.remove(file_for_excel_write_test)
}
write_pathways_to_excel(filePath = file_for_excel_write_test,
                        enriched_data_list = enrichedList$enriched_data_list)

#missing values, samplesID, nc, Distance, ClusterMethod (user Inputs), they are assumed to be validated by frontend/ enforced by frontend

plotting_params = do_cluster(enrichedList$KEGG_MAT,
                        samplesID = c("All"),
                        nc=3,
                        Distance="comb",
                        ClusterMethod="ward.D") #returns list object of 4

#FOR do_gene_heat_map hls, cls, pheno & kegg hierary are missing (also from enrichment function??)
heatmap_params = do_gene_heat_map(KEGG_MAT_GENES = enrichedList$KEGG_MAT_GENES,
                                  kegg_hierarchy = enrichedList$hierarchy,
                                  GList =  GList,
                                  KEGG_MAT = enrichedList$KEGG_MAT,
                                  pheno=pheno, #missing
                                  hls = plotting_params$hls,
                                  cls = plotting_params$cls,
                                  aspectRatio=F, #missing
                                  selPath = "Metabolism",
                                  levelGene=1,
                                  selScoreType = "lfc") #returns list

# make a new function
grid_draw_plot = grid_draw(toPlotGenes = heatmap_params$toPlotGenes,
                           kegg_nano_1 = heatmap_params$res_collapsed,
                           mat_to_Plot = heatmap_params$mat_to_Plot,
                           exp_ann = plotting_params$exp_ann)


plot_count <- 1

toPlot <- base::list()

toPlot[[plot_count]] <- create_map_plot(
                          lev1_content = "All",
                          lev2_content = "All",
                          lev3_content = "All",
                          experiment_ann = pheno,
                          kegg_hierarchy = enrichedList$hierarchy,
                          toPlotMap = enrichedList$KEGG_MAT,
                          samplesID = "All",
                          collapse_level = 1,
                          doGrouping = TRUE,
                          aspectRatio = F,
                          continuous = FALSE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
                          lev1_content = "All",
                          lev2_content = "All",
                          lev3_content = "All",
                          experiment_ann = plotting_params$exp_ann,
                          kegg_hierarchy = enrichedList$hierarchy,
                          toPlotMap = plotting_params$clust_mat,
                          samplesID = "All",
                          collapse_level = 1,
                          doGrouping = FALSE,
                          aspectRatio = F,
                          continuous = FALSE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = "All",
  lev2_content = "All",
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = "All",
  collapse_level = 1,
  doGrouping = FALSE,
  aspectRatio = FALSE,
  continuous = FALSE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = "All",
  lev2_content = "All",
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = "All",
  collapse_level = 1,
  doGrouping = FALSE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = "All",
  lev2_content = "All",
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = "All",
  collapse_level = 2,
  doGrouping = FALSE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = "All",
  lev2_content = "All",
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = "All",
  collapse_level = 2,
  doGrouping = TRUE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = "All",
  lev2_content = "All",
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = "All",
  collapse_level = 3,
  doGrouping = FALSE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = "All",
  lev2_content = "All",
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = c("CSA","WY","CSPT","DM","OA","DIDP","CA"),
  collapse_level = 2,
  doGrouping = FALSE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = "All",
  lev2_content = c("Immune diseases","Lipid metabolism","Immune system","Transcription","Aging"),
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = c("CSA","WY","CSPT","DM","OA","DIDP","CA"),
  collapse_level = 2,
  doGrouping = FALSE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = "All",
  lev2_content = c("Immune diseases","Lipid metabolism","Immune system","Transcription","Aging"),
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = c("CSA","WY","CSPT","DM","OA","DIDP","CA"),
  collapse_level = 2,
  doGrouping = TRUE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = "All",
  lev2_content = c("Immune diseases","Lipid metabolism","Immune system","Transcription","Aging"),
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = c("CSA","WY","CSPT","DM","OA","DIDP","CA"),
  collapse_level = 3,
  doGrouping = TRUE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = c("Organismal Systems","Human Diseases"),
  lev2_content = "All",
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = c("CSA","WY","CSPT","DM","OA","DIDP","CA"),
  collapse_level = 1,
  doGrouping = FALSE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = c("Organismal Systems","Human Diseases"),
  lev2_content = "All",
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = c("CSA","WY","CSPT","DM","OA","DIDP","CA"),
  collapse_level = 2,
  doGrouping = TRUE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = c("Organismal Systems","Human Diseases"),
  lev2_content = c("Immune diseases","Lipid metabolism","Immune system","Transcription","Aging","Circulatory system","Substance dependence"),
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = c("CSA","WY","CSPT","DM","OA","DIDP","CA"),
  collapse_level = 2,
  doGrouping = TRUE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = c("Organismal Systems","Human Diseases"),
  lev2_content = c("Immune diseases","Lipid metabolism","Immune system","Transcription","Aging","Circulatory system","Substance dependence"),
  lev3_content = "All",
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = c("CSA","WY","CSPT","DM","OA","DIDP","CA"),
  collapse_level = 3,
  doGrouping = TRUE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))

plot_count <- plot_count + 1
toPlot[[plot_count]] <- create_map_plot(
  lev1_content = c("Organismal Systems","Human Diseases"),
  lev2_content = c("Immune diseases","Lipid metabolism","Immune system","Transcription","Aging","Circulatory system","Substance dependence"),
  lev3_content = c("Rheumatoid arthritis", "Nicotine addiction", "Alcoholism"),
  experiment_ann = plotting_params$exp_ann,
  kegg_hierarchy = enrichedList$hierarchy,
  toPlotMap = plotting_params$clust_mat,
  samplesID = c("CSA","WY","CSPT","DM","OA","DIDP","CA"),
  collapse_level = 3,
  doGrouping = TRUE,
  aspectRatio = F,
  continuous = TRUE)
print(paste("After create_map_plot call",plot_count))
graphics::plot(toPlot[[plot_count]])
print(paste("After cgraphics::plot call",plot_count))
