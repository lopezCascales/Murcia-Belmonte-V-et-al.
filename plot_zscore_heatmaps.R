#!/usr/bin/env Rscript
## ─────────────────────────────────────────────────────────────────────────────
## Z-score heatmaps para interacciones RNA-RNA lncRNA × mRNA
## Genera las figuras principales del análisis RNAplex
## Uso: Rscript plot_zscore_heatmaps.R
## ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

RESULTS_CSV <- "outputs/rna_interactions/zscore_full_results.csv"
DINUC_CSV   <- "outputs/rna_interactions/zscore_dinucleotide_results.csv"
OUT_DIR     <- "outputs/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ── Cargar resultados ─────────────────────────────────────────────────────────
message("Cargando resultados Z-score...")
df <- read_csv(RESULTS_CSV, show_col_types = FALSE)

## Cargar dinucleotídico si existe
has_dinuc <- file.exists(DINUC_CSV)
if (has_dinuc) {
  dinuc <- read_csv(DINUC_CSV, show_col_types = FALSE) %>%
    select(query, target, dinuc_z, survives) %>%
    rename(lncRNA = query)
}

## ── Paleta de colores del paper ───────────────────────────────────────────────
CONTRA_COL <- "#C9375A"
IPSI_COL   <- "#3B8BD4"
KEY_LNCS   <- c("H19", "Mirg", "Meg3", "Rian", "Lhx1os")

## ── Preparar matriz para heatmap ─────────────────────────────────────────────
## Pivotar a matriz lncRNA × target
zscore_wide <- df %>%
  select(query, target, zscore, group) %>%
  pivot_wider(names_from = target, values_from = zscore, id_cols = query) %>%
  column_to_rownames("query")

zscore_mat <- as.matrix(zscore_wide)
zscore_mat[is.na(zscore_mat)] <- 0

## Clampar a [-5, 0] para visualización (solo interacciones negativas son relevantes)
zscore_mat_clamped <- pmax(zscore_mat, -5)

## ── Anotación de columnas: grupos de genes ───────────────────────────────────
## Obtener el grupo de cada target del dataframe original
target_groups <- df %>%
  distinct(target, group) %>%
  deframe()

col_groups <- target_groups[colnames(zscore_mat_clamped)]
col_groups[is.na(col_groups)] <- "unknown"

## Clasificar en categorías biológicas
make_category <- function(g) {
  case_when(
    str_detect(toupper(g), "CONTRA") & !str_detect(toupper(g), "TF|EPIG") ~ "CONTRA axon",
    str_detect(toupper(g), "IPSI")   & !str_detect(toupper(g), "TF|EPIG") ~ "IPSI axon",
    str_detect(toupper(g), "TF|TRANSCRI") ~ "TF nuclear",
    str_detect(toupper(g), "EPIG|CHROM")  ~ "Epigenetic",
    str_detect(toupper(g), "CPE|RBP|RISC|NMD") ~ "CPE/RBP",
    TRUE ~ "other"
  )
}

col_category <- make_category(col_groups)

cat_colors <- c(
  "CONTRA axon" = CONTRA_COL,
  "IPSI axon"   = IPSI_COL,
  "TF nuclear"  = "#888780",
  "Epigenetic"  = "#B4B2A9",
  "CPE/RBP"     = "#3B6D11",
  "other"       = "#D3D1C7"
)

top_anno <- HeatmapAnnotation(
  category = col_category,
  col = list(category = cat_colors),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = 8)
)

## ── Significancia: asteriscos ─────────────────────────────────────────────────
## Construir matriz booleana de significancia
sig_wide <- df %>%
  mutate(sig = zscore < -1.96) %>%
  select(query, target, sig) %>%
  pivot_wider(names_from = target, values_from = sig, id_cols = query) %>%
  column_to_rownames("query")

sig_mat <- as.matrix(sig_wide)[rownames(zscore_mat_clamped), colnames(zscore_mat_clamped)]
sig_mat[is.na(sig_mat)] <- FALSE

## ── Colores del heatmap ───────────────────────────────────────────────────────
col_fun_zscore <- colorRamp2(
  breaks = c(-5, -2, -1, 0),
  colors = c("#67001F", "#D6604D", "#FDDBC7", "white")
)

## ── FIGURA 1a: Moléculas de guía axonal ──────────────────────────────────────
axon_targets <- df %>%
  filter(str_detect(toupper(group), "AXON|CONTRA|IPSI")) %>%
  filter(!str_detect(toupper(group), "TF|EPIG|NUCLEAR")) %>%
  pull(target) %>% unique()

mat_axon <- zscore_mat_clamped[KEY_LNCS[KEY_LNCS %in% rownames(zscore_mat_clamped)],
                                axon_targets[axon_targets %in% colnames(zscore_mat_clamped)]]

ht_axon <- Heatmap(
  mat_axon,
  name = "Z-score",
  col  = col_fun_zscore,
  top_annotation = top_anno[, colnames(mat_axon)],
  
  row_names_gp    = gpar(fontsize = 9, fontface = "italic"),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  
  cluster_rows    = FALSE,  ## orden fijo de lncRNAs
  cluster_columns = TRUE,
  clustering_method_columns = "ward.D2",
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (isTRUE(sig_mat[rownames(mat_axon)[i], colnames(mat_axon)[j]]))
      grid.text("*", x, y, gp = gpar(fontsize = 8, col = "white", fontface = "bold"))
  },
  
  column_title = "Axonal guidance & adhesion molecules",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-5, -4, -3, -2, -1, 0),
    labels = c("≤−5", "−4", "−3", "−2", "−1", "0"),
    direction = "horizontal"
  )
)

## ── FIGURA: scatter mono vs dinuc Z-score ────────────────────────────────────
if (has_dinuc) {
  scatter_data <- dinuc %>%
    left_join(df %>% select(query, target, zscore) %>% rename(lncRNA = query, mono_z = zscore),
              by = c("lncRNA", "target")) %>%
    mutate(
      status = case_when(
        survives & lncRNA == "H19" ~ "H19 survives",
        survives                   ~ "other survives",
        TRUE                       ~ "artifact"
      ),
      label = ifelse(abs(mono_z - dinuc_z) > 2 | abs(dinuc_z) > 4, 
                     paste0(lncRNA, "×", target), NA_character_)
    )
  
  p_scatter <- ggplot(scatter_data, aes(x = mono_z, y = dinuc_z)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
    geom_hline(yintercept = -1.96, linetype = "dotted", color = CONTRA_COL, linewidth = 0.4) +
    geom_vline(xintercept = -1.96, linetype = "dotted", color = IPSI_COL,   linewidth = 0.4) +
    geom_point(aes(color = status, shape = lncRNA), size = 2.5, alpha = 0.85) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15,
                    segment.size = 0.3, segment.color = "gray70") +
    scale_color_manual(
      values = c("H19 survives" = "#E85D24",
                 "other survives" = "#3B6D11",
                 "artifact" = "#B4B2A9"),
      name = NULL
    ) +
    scale_shape_manual(values = c(H19 = 16, Mirg = 17, Meg3 = 15, Rian = 18),
                       name = "lncRNA") +
    labs(
      x = "Mononucleotide Z-score",
      y = "Dinucleotide Z-score",
      title = "Mononucleotide vs dinucleotide shuffle",
      subtitle = sprintf("N=%d pairs tested | %d/%d survive dinucleotide (%.0f%%)",
                         nrow(scatter_data),
                         sum(scatter_data$survives, na.rm = TRUE),
                         nrow(scatter_data),
                         100 * mean(scatter_data$survives, na.rm = TRUE))
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "right",
          plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8, color = "gray40"))
  
  ggsave(file.path(OUT_DIR, "scatter_mono_vs_dinuc.pdf"),
         p_scatter, width = 7, height = 5, device = "pdf")
  message("Guardado: scatter_mono_vs_dinuc.pdf")
}

## ── FIGURA: IPSI vs CONTRA preference (Fig 3 del paper) ─────────────────────
preference_data <- df %>%
  filter(str_detect(toupper(group), "AXON|CONTRA|IPSI")) %>%
  filter(!str_detect(toupper(group), "TF|EPIG|NUCLEAR")) %>%
  mutate(laterality = ifelse(str_detect(toupper(group), "CONTRA"), "CONTRA", "IPSI")) %>%
  group_by(query, laterality) %>%
  summarise(
    mean_z  = mean(zscore, na.rm = TRUE),
    se_z    = sd(zscore, na.rm = TRUE) / sqrt(n()),
    n_sig   = sum(zscore < -1.96, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(query = factor(query, levels = KEY_LNCS))

p_pref <- ggplot(preference_data, aes(x = laterality, y = mean_z,
                                       fill = laterality, size = n_sig)) +
  geom_col(aes(fill = laterality), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
                width = 0.15, linewidth = 0.4) +
  geom_point(aes(size = n_sig), shape = 21, fill = "white",
             color = "gray30", stroke = 0.4) +
  scale_fill_manual(values = c(CONTRA = CONTRA_COL, IPSI = IPSI_COL)) +
  scale_size_continuous(range = c(2, 8), name = "# significant\ninteractions") +
  facet_wrap(~query, nrow = 1) +
  labs(
    x = NULL, y = "Mean Z-score ± SE",
    title = "CONTRA preference: lncRNAs interact preferentially with CONTRA targets",
    subtitle = "Axonal guidance/adhesion molecules only"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position   = "right",
    strip.text        = element_text(face = "italic", size = 9),
    plot.title        = element_text(size = 10, face = "bold"),
    plot.subtitle     = element_text(size = 8, color = "gray40"),
    axis.text.x       = element_text(size = 8),
    panel.spacing     = unit(0.8, "lines")
  )

ggsave(file.path(OUT_DIR, "contra_preference_barplot.pdf"),
       p_pref, width = 10, height = 4.5, device = "pdf")
message("Guardado: contra_preference_barplot.pdf")

## ── Guardar heatmap axonal ────────────────────────────────────────────────────
pdf(file.path(OUT_DIR, "zscore_heatmap_axonal_guidance.pdf"), width = 14, height = 5)
draw(ht_axon, heatmap_legend_side = "bottom")
dev.off()
message("Guardado: zscore_heatmap_axonal_guidance.pdf")

message("\nTodos los outputs guardados en: ", OUT_DIR)
