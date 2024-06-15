source("process_data.R")

st_id <- "G2"
sc_data <- load_sc_wu('HER2+')
st_data <- load_st_her2p(st_id)

dot_st <- DOTr::setup.srt(st_data$counts, st_data$meta[, c("x", "y")], radius = sqrt(2)+0.1)
dot_sc <- DOTr::setup.ref(sc_data$counts, sc_data$meta[, "celltype_major"], max_genes = 10000)

dot <- DOTr::create.DOT(dot_st, dot_sc)
dot <- DOTr::run.DOT.lowresolution(dot, ratios_weight = 0.25, max_size = 200)
