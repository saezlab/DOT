source('process_data.R')

ref_type <- "brain"
vis_index <- "1"

libd_data <- load_libd_experiment(ref_type, vis_index)

dot_ref <- DOTr::setup.ref(ref_data = libd_data$ref$counts, ref_annotations = libd_data$ref$spatial$layer, max_genes = 10000)

dot_srt <- DOTr::setup.srt(srt_data = libd_data$srt$counts, srt_coords = libd_data$srt$spatial[, c('Row', 'Col')])

dot <- DOTr::create.DOT(dot_srt, dot_ref)

dot <- DOTr::run.DOT.lowresolution(dot, ratios_weight = ifelse(ref_type == 'aggregated', 0.25, 1), max_size = 20)
