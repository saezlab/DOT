source("process_data.R")

vis_id <- "1142243F"
sc_data <- load_sc_wu('TNBC')
vis_data <- load_vis_tnbc(vis_id)

dot_vis <- DOT::setup.srt(vis_data$counts, vis_data$meta[, c("Col", "Row")], radius = sqrt(2)+0.1)
dot_sc <- DOT::setup.ref(sc_data$counts, sc_data$meta[, "celltype_major"], max_genes = 10000)

dot <- DOT::create.DOT(dot_vis, dot_sc)
dot <- DOT::run.DOT.lowresolution(dot, ratios_weight = 0.25, max_size = 20)