#!/usr/bin/env Rscript
## ─────────────────────────────────────────────────────────────────────────────
## Spearman correlation heatmap: lncRNAs × DEGs axonales
## Reproduce y mejora la figura de correlación del paper
## Uso: Rscript spearman_heatmap_publication.R
## ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

## ── Parámetros ────────────────────────────────────────────────────────────────
COUNTS_FILE   <- "lncRNAcountsTPM.txt"   ## TPM de lncRNAs (12 muestras)
MRNA_FILE     <- "mRNATPMcounts.txt"     ## TPM de mRNAs (mismas 12 muestras)
RHO_THRESHOLD <- 0.4                     ## |rho| mínimo para filtrar pares
PADJ_CUTOFF   <- 0.1                     ## Para filtrar DEGs significativos
OUT_DIR       <- "outputs/figures"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ── Orden de muestras (12): E13.5 + E15.5 × CH + retina × replicados ─────────
## Ajusta este vector al orden real de tus columnas
SAMPLE_ORDER <- c(
  "e13_5_CH1_both.txt", "e13_5_CH2_both.txt", "e13_5_CH4_both.txt",
  "htp19_BTA_both.txt", "htp20_BTA_both.txt", "htp21_BTA_both.txt",
  "e15_5_CH1_both.txt", "e15_5_CH2_both.txt", "e15_5_CH3_both.txt",
  "htp55_STA_both.txt", "htp56_STA_both.txt",
  "e15_5_ret1_both.txt"
)

## ── Cargar datos ──────────────────────────────────────────────────────────────
message("Cargando datos...")
lnc_raw  <- read.table(COUNTS_FILE,  header = TRUE, row.names = 1, sep = "\t")
mrna_raw <- read.table(MRNA_FILE, header = TRUE, row.names = 1, sep = "\t")

## Seleccionar y ordenar muestras comunes
shared_samples <- intersect(SAMPLE_ORDER, intersect(colnames(lnc_raw), colnames(mrna_raw)))
lnc_mat  <- as.matrix(lnc_raw[,  shared_samples])
mrna_mat <- as.matrix(mrna_raw[, shared_samples])

message(sprintf("lncRNAs: %d | mRNAs: %d | muestras: %d",
                nrow(lnc_mat), nrow(mrna_mat), ncol(lnc_mat)))

## ── Función: correlación Spearman por pares ───────────────────────────────────
compute_spearman_matrix <- function(lnc_mat, mrna_mat, rho_thr = 0.4) {
  ## Solo lncRNAs con varianza > 0
  lnc_var <- apply(lnc_mat, 1, var)
  lnc_use <- lnc_mat[lnc_var > 0, ]
  
  ## Solo mRNAs con varianza > 0
  mrn_var <- apply(mrna_mat, 1, var)
  mrn_use <- mrna_mat[mrn_var > 0, ]
  
  message(sprintf("Calculando %d × %d correlaciones...",
                  nrow(lnc_use), nrow(mrn_use)))
  
  rho_mat <- cor(t(lnc_use), t(mrn_use), method = "spearman")
  
  ## Filtrar: solo filas/columnas con al menos un |rho| >= threshold
  row_keep <- apply(abs(rho_mat), 1, max, na.rm = TRUE) >= rho_thr
  col_keep <- apply(abs(rho_mat), 2, max, na.rm = TRUE) >= rho_thr
  
  rho_mat[row_keep, col_keep]
}

rho_mat <- compute_spearman_matrix(lnc_mat, mrna_mat, RHO_THRESHOLD)
message(sprintf("Matriz filtrada: %d lncRNAs × %d mRNAs", nrow(rho_mat), ncol(rho_mat)))

## ── Anotación de columnas (mRNAs): IPSI vs CONTRA ────────────────────────────
## Carga la clasificación desde el Excel de candidatos
## Si no tienes el CSV, define manualmente los vectores
CONTRA_GENES <- c(
  "Sema5a", "Islr2", "Plxna4", "Rgs7bp", "Nrp2", "Ephb1", "Robo2",
  "Bcl11a", "Bcl11b", "Neurod2", "Fezf1", "Ago2", "Ago3", "Smg1",
  "Nova1", "Kmt2a", "Epha4", "Rian", "Tenm2"
)
IPSI_GENES <- c(
  "H19", "Lhx1", "Grin2b", "Sema3f", "Nr2f2", "Tenm3", "Robo3",
  "Plxna1", "Gda", "Tenm1"
)

col_group <- case_when(
  colnames(rho_mat) %in% CONTRA_GENES ~ "CONTRA",
  colnames(rho_mat) %in% IPSI_GENES   ~ "IPSI",
  TRUE                                 ~ "other"
)

## ── Colores consistentes con el paper ────────────────────────────────────────
col_anno_colors <- list(
  group = c(CONTRA = "#C9375A", IPSI = "#3B8BD4", other = "#888780")
)

## Anotación superior (mRNAs)
top_anno <- HeatmapAnnotation(
  group = col_group,
  col   = col_anno_colors,
  annotation_name_side = "left",
  show_legend = TRUE
)

## Anotación lateral (lncRNAs): destacar los 5 lncRNAs principales
KEY_LNCRNAS <- c("H19", "Mirg", "Meg3", "Rian", "Lhx1os")
row_highlight <- ifelse(rownames(rho_mat) %in% KEY_LNCRNAS, "key", "other")

left_anno <- rowAnnotation(
  lncRNA = row_highlight,
  col = list(lncRNA = c(key = "#E85D24", other = "#D3D1C7")),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  width = unit(3, "mm")
)

## ── Escala de colores ─────────────────────────────────────────────────────────
col_fun <- colorRamp2(
  breaks = c(-1, -0.4, 0, 0.4, 1),
  colors = c("#053061", "#4393C3", "white", "#D6604D", "#67001F")
)

## ── Heatmap principal ─────────────────────────────────────────────────────────
## Solo mostrar nombres de los lncRNAs clave; el resto queda en blanco
row_labels <- ifelse(rownames(rho_mat) %in% KEY_LNCRNAS, rownames(rho_mat), "")

ht <- Heatmap(
  rho_mat,
  name = "Spearman\nrho",
  
  ## Colores
  col = col_fun,
  
  ## Filas (lncRNAs)
  row_labels         = row_labels,
  row_names_gp       = gpar(fontsize = 8, fontface = "italic"),
  show_row_dend      = TRUE,
  cluster_rows       = TRUE,
  clustering_method_rows = "ward.D2",
  left_annotation    = left_anno,
  
  ## Columnas (mRNAs)
  column_names_gp    = gpar(fontsize = 7),
  show_column_names  = ncol(rho_mat) <= 60,  ## solo si hay pocos genes
  top_annotation     = top_anno,
  cluster_columns    = TRUE,
  clustering_method_columns = "ward.D2",
  
  ## Apariencia
  cell_fun = function(j, i, x, y, width, height, fill) {
    ## Marcar con asterisco si |rho| > 0.6
    if (abs(rho_mat[i, j]) > 0.6)
      grid.text("*", x, y, gp = gpar(fontsize = 6, col = "black"))
  },
  
  heatmap_legend_param = list(
    title       = "Spearman rho",
    title_gp    = gpar(fontsize = 9),
    labels_gp   = gpar(fontsize = 8),
    legend_height = unit(3, "cm"),
    at = c(-1, -0.5, 0, 0.5, 1)
  ),
  
  use_raster = TRUE,
  raster_quality = 3
)

## ── Guardar ──────────────────────────────────────────────────────────────────
out_pdf <- file.path(OUT_DIR, "spearman_lncrna_degs_publication.pdf")
pdf(out_pdf, width = 12, height = 8)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
message("Guardado: ", out_pdf)

## ── Tabla de pares significativos ────────────────────────────────────────────
rho_long <- rho_mat %>%
  as.data.frame() %>%
  rownames_to_column("lncRNA") %>%
  pivot_longer(-lncRNA, names_to = "mRNA", values_to = "rho") %>%
  filter(abs(rho) >= RHO_THRESHOLD) %>%
  mutate(
    direction = ifelse(rho > 0, "positive", "negative"),
    group_mrna = case_when(
      mRNA %in% CONTRA_GENES ~ "CONTRA",
      mRNA %in% IPSI_GENES   ~ "IPSI",
      TRUE                   ~ "other"
    )
  ) %>%
  arrange(desc(abs(rho)))

out_tsv <- file.path(OUT_DIR, "spearman_significant_pairs.tsv")
write_tsv(rho_long, out_tsv)
message(sprintf("Pares significativos (|rho| >= %.1f): %d → %s",
                RHO_THRESHOLD, nrow(rho_long), out_tsv))

## ── Resumen por lncRNA ────────────────────────────────────────────────────────
summary_tbl <- rho_long %>%
  group_by(lncRNA, group_mrna) %>%
  summarise(
    n_pairs   = n(),
    mean_rho  = round(mean(rho), 3),
    n_pos     = sum(rho > 0),
    n_neg     = sum(rho < 0),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = group_mrna, values_from = c(n_pairs, mean_rho),
              values_fill = 0)

print(summary_tbl)
write_tsv(summary_tbl, file.path(OUT_DIR, "spearman_summary_by_lncrna.tsv"))
