######################################################################
##  Figures for all 9 lncRNAs: original 5 + 4 new
##  Handles different target coverage (64 vs 10 targets)
##  Run AFTER rna_zscore_new_rnaplex.py completes
######################################################################

library(tidyverse)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(conflicted)
conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename)

dir.create("outputs/rna_interactions/figures", recursive=TRUE, showWarnings=FALSE)

col_ipsi   <- "#4BBFCF"
col_contra <- "#C0395A"

lnc_col <- c(
  H19        = "#1D9E75",
  Lhx1os     = "#9FE1CB",
  Mirg       = "#C0395A",
  Meg3       = "#D4537E",
  Rian       = "#E89BB0",
  Kcnq1ot1   = "#8B0000",
  `Mir9-3hg` = "#FF6B6B",
  Rpph1      = "#006994",
  Arhgap20os = "#4FC3F7"
)

lnc_order <- c("H19","Lhx1os","Rpph1","Arhgap20os",
                "Mirg","Meg3","Rian","Kcnq1ot1","Mir9-3hg")

lnc_class_map <- c(
  H19="lncAI", Lhx1os="lncAI", Rpph1="lncAI", Arhgap20os="lncAI",
  Mirg="lncAC", Meg3="lncAC", Rian="lncAC", Kcnq1ot1="lncAC", `Mir9-3hg`="lncAC"
)

reclassify <- tribble(
  ~target,   ~new_group,
  "Lhx2",    "IPSI_TF",   "Unc5c",  "CONTRA",
  "Plxna1",  "CONTRA",    "Pou4f1", "CONTRA_TF",
  "Sox4",    "CONTRA_TF", "Sox11",  "CONTRA_TF",
  "Isl2",    "CONTRA_TF", "Cxcr4",  "SHARED",
  "Igf2",    "CIS_LOCUS"
)

NUCLEAR <- c("Neurod2","Fezf1","Lhx1","Zic4","Nr2f2","Atoh7",
             "Bcl11a","Bcl11b","Kmt2a","Kmt2b","Hdac4","Ep300",
             "Pou4f2","Sox4","Sox11","Isl2","Pou4f1","Foxd1",
             "Zic2","Lhx2","Nr2f1")
RBP_CPE <- c("Ago2","Ago3","Nova1","Qk","Smg1","Tnrc6a","Tnrc6b","Cnot3")

## ── Load original 5 lncRNAs ──────────────────────────────────────────
orig_csv <- "outputs/rna_interactions/zscore_full_results.csv"
if (!file.exists(orig_csv)) stop("Run rna_zscore_full.py first")
res_orig <- read.csv(orig_csv, stringsAsFactors=FALSE)
cat("Original results:", nrow(res_orig), "pairs\n")

## ── Load new 4 lncRNAs ───────────────────────────────────────────────
new_csv <- "outputs/rna_interactions/zscore_new_lncrnas_results.csv"
if (!file.exists(new_csv)) stop("Run rna_zscore_new_rnaplex.py first")
res_new <- read.csv(new_csv, stringsAsFactors=FALSE)
cat("New results:", nrow(res_new), "pairs\n")
cat("New lncRNAs:", paste(unique(res_new$query), collapse=", "), "\n")

## ── Find shared targets (only targets present in BOTH datasets) ──────
shared_targets <- intersect(
  unique(res_orig$target),
  unique(res_new$target)
)
cat("Shared targets:", length(shared_targets), "->", paste(shared_targets, collapse=", "), "\n")

## ── Combine — restrict orig to shared targets for fair comparison ────
res_orig_shared <- dplyr::filter(res_orig, target %in% shared_targets)
res_all <- dplyr::bind_rows(res_orig_shared, res_new) %>%
  dplyr::left_join(reclassify, by="target") %>%
  dplyr::mutate(
    group = dplyr::coalesce(new_group, group),
    mol_type = dplyr::case_when(
      target %in% NUCLEAR  ~ "NUCLEAR",
      target %in% RBP_CPE  ~ "RBP_CPE",
      group %in% c("CIS_LOCUS","SHARED") ~ group,
      group %in% c("IPSI_heatmap","IPSI_new","IPSI","IPSI_TF") ~ "AXONAL_IPSI",
      group %in% c("CONTRA_heatmap","CONTRA_new","CONTRA","CONTRA_TF") ~ "AXONAL_CONTRA",
      TRUE ~ "OTHER"
    ),
    specific   = zscore < -2.0,
    lnc_class  = dplyr::recode(query, !!!lnc_class_map),
    query      = factor(query, levels=lnc_order)
  ) %>%
  dplyr::select(-new_group) %>%
  dplyr::filter(!is.na(query))

res_axon <- dplyr::filter(res_all, mol_type %in% c("AXONAL_IPSI","AXONAL_CONTRA"))
cat("\nAxonal pairs for shared targets:", nrow(res_axon), "\n")

## Print significant hits from new lncRNAs
cat("\n=== Significant hits (Z < -2) from NEW lncRNAs ===\n")
new_sig <- dplyr::filter(res_all,
                          query %in% c("Kcnq1ot1","Mir9-3hg","Rpph1","Arhgap20os"),
                          specific)
if (nrow(new_sig) > 0) {
  print(new_sig %>% dplyr::select(query, target, group, zscore, pval) %>%
          dplyr::arrange(query, zscore))
} else {
  cat("No significant interactions found for new lncRNAs\n")
  cat("(200nt query screening — expected for highly truncated sequences)\n")
}

## ════════════════════════════════════════════════════════════════════
## FIG A: Heatmap all 9 lncRNAs — shared targets only
## ════════════════════════════════════════════════════════════════════
mat <- res_axon %>%
  dplyr::select(query, target, zscore) %>%
  tidyr::pivot_wider(names_from=target, values_from=zscore,
                     values_fn=mean) %>%
  tibble::column_to_rownames("query") %>%
  as.matrix()

mat_cl <- pmax(mat, -5)

## Annotations
col_ann <- res_axon %>%
  dplyr::distinct(target, mol_type) %>%
  dplyr::mutate(group=ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA")) %>%
  dplyr::select(target, group) %>%
  tibble::column_to_rownames("target")

row_ann <- data.frame(
  lnc_class = lnc_class_map[rownames(mat)],
  row.names  = rownames(mat)
) %>% dplyr::filter(!is.na(lnc_class))

row_ann <- row_ann[rownames(row_ann) %in% rownames(mat), , drop=FALSE]

## Row colors using lnc_col
row_cols <- lnc_col[rownames(mat)]
row_cols <- row_cols[!is.na(row_cols)]

sig_mat <- ifelse(mat < -2, "*", "")

hm_col <- colorRampPalette(
  c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#B2182B"))(100)

ann_colors <- list(
  group     = c(IPSI=col_ipsi, CONTRA=col_contra),
  lnc_class = c(lncAI=col_ipsi, lncAC=col_contra)
)

n_shared <- ncol(mat)
cat(sprintf("\nHeatmap: %d lncRNAs x %d shared targets\n",
            nrow(mat), n_shared))

pheatmap::pheatmap(
  mat_cl,
  color         = hm_col,
  breaks        = seq(-5, 3, length.out=101),
  annotation_col = col_ann,
  annotation_row = row_ann,
  annotation_colors = ann_colors,
  cluster_rows  = FALSE,
  cluster_cols  = TRUE,
  border_color  = NA,
  fontsize_row  = 9,
  fontsize_col  = 8,
  display_numbers = sig_mat,
  number_color  = "black",
  fontsize_number = 9,
  main = sprintf(
    "RNA-RNA interaction Z-score — all 9 lncRNAs x %d shared axonal targets (* p<0.05)",
    n_shared
  ),
  filename = "outputs/rna_interactions/figures/FigNew_heatmap_all9.pdf",
  width    = max(8, n_shared * 0.5 + 2),
  height   = 6
)
cat("Saved: FigNew_heatmap_all9.pdf\n")

## ════════════════════════════════════════════════════════════════════
## FIG B: CONTRA preference — all 9 lncRNAs
## ════════════════════════════════════════════════════════════════════
sum9 <- res_axon %>%
  dplyr::mutate(group_plot=ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA")) %>%
  dplyr::group_by(query, lnc_class, group_plot) %>%
  dplyr::summarise(
    mean_z = mean(zscore, na.rm=TRUE),
    se_z   = sd(zscore, na.rm=TRUE) / sqrt(dplyr::n()),
    n_sig  = sum(specific, na.rm=TRUE),
    n      = dplyr::n(),
    .groups="drop"
  )

## Wilcoxon per lncRNA
wt9 <- res_axon %>%
  dplyr::mutate(group_plot=ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA")) %>%
  dplyr::group_by(query, lnc_class) %>%
  dplyr::summarise(
    p_val = tryCatch(
      wilcox.test(zscore[group_plot=="IPSI"],
                  zscore[group_plot=="CONTRA"],
                  alternative="greater")$p.value,
      error=function(e) NA_real_
    ),
    .groups="drop"
  ) %>%
  dplyr::mutate(
    sig_label = dplyr::case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      !is.na(p_val) ~ "n.s.",
      TRUE ~ ""
    )
  )

cat("\nWilcoxon results (all 9 lncRNAs):\n")
print(wt9 %>% dplyr::select(query, lnc_class, p_val, sig_label))

pB <- ggplot(sum9, aes(x=group_plot, y=mean_z, colour=query, group=query)) +
  geom_hline(yintercept=0, colour="grey80", linewidth=0.3) +
  geom_hline(yintercept=-2, colour="grey40", linewidth=0.5, linetype="dashed") +
  geom_errorbar(aes(ymin=mean_z-se_z, ymax=mean_z+se_z),
                width=0.12, alpha=0.5, linewidth=0.7) +
  geom_line(alpha=0.7, linewidth=0.9) +
  geom_point(aes(size=n_sig), alpha=0.95) +
  ggrepel::geom_text_repel(
    data=dplyr::filter(sum9, mean_z < -0.5 & group_plot=="CONTRA"),
    aes(label=paste0(query,"\n(",n_sig,"/",n,")")),
    size=2.5, nudge_x=0.2, box.padding=0.3, max.overlaps=12
  ) +
  scale_colour_manual(values=lnc_col, name="lncRNA") +
  scale_size_continuous(range=c(2,7), breaks=c(0,1,2,3),
                        name="n significant") +
  scale_x_discrete(labels=c(IPSI="Ipsi targets", CONTRA="Contra targets")) +
  facet_wrap(~lnc_class, nrow=1,
    labeller=labeller(lnc_class=c(
      lncAI="lncAI (ipsi: H19, Lhx1os, Rpph1, Arhgap20os)",
      lncAC="lncAC (contra: Mirg, Meg3, Rian, Kcnq1ot1, Mir9-3hg)"
    ))
  ) +
  labs(
    x=NULL, y="Mean Z-score (+/- SE)",
    title="CONTRA preference: all 9 axonal lncRNAs",
    subtitle=paste0(
      "Compared on ", n_shared, " shared targets\n",
      "New lncRNAs (Kcnq1ot1, Mir9-3hg, Rpph1, Arhgap20os) screened with 200nt query"
    )
  ) +
  theme_bw(base_size=10) +
  theme(
    panel.grid.major.x = element_blank(),
    strip.background   = element_rect(fill="grey93"),
    strip.text         = element_text(face="bold", size=8),
    legend.position    = "right"
  )

ggsave("outputs/rna_interactions/figures/FigNew_contra_preference_9lncrnas.pdf",
       pB, width=12, height=5)
ggsave("outputs/rna_interactions/figures/FigNew_contra_preference_9lncrnas.svg",
       pB, width=12, height=5)
cat("Saved: FigNew_contra_preference_9lncrnas.pdf\n")

## ════════════════════════════════════════════════════════════════════
## FIG C: Boxplot all 9 lncRNAs
## ════════════════════════════════════════════════════════════════════
res_box <- res_axon %>%
  dplyr::mutate(group_plot=ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA"))

pC <- ggplot(res_box, aes(x=group_plot, y=zscore, fill=group_plot)) +
  geom_hline(yintercept=0, colour="grey80", linewidth=0.3) +
  geom_hline(yintercept=-2, colour="grey40", linewidth=0.5, linetype="dashed") +
  geom_boxplot(alpha=0.3, outlier.shape=NA, width=0.5,
               linewidth=0.5, colour="grey30") +
  geom_jitter(aes(colour=group_plot, size=specific),
              width=0.16, alpha=0.7) +
  geom_text(data=wt9,
            aes(x=1.5, y=Inf, label=sig_label),
            inherit.aes=FALSE, vjust=1.5, size=4.5, colour="grey25") +
  ggrepel::geom_text_repel(
    data=dplyr::filter(res_box, specific),
    aes(label=target, colour=group_plot),
    size=2.0, box.padding=0.2, max.overlaps=5, show.legend=FALSE
  ) +
  scale_fill_manual(values=c(IPSI=col_ipsi, CONTRA=col_contra), guide="none") +
  scale_colour_manual(values=c(IPSI=col_ipsi, CONTRA=col_contra), guide="none") +
  scale_size_manual(values=c("FALSE"=1.5,"TRUE"=3.5),
                    labels=c("FALSE"="n.s.","TRUE"="p<0.05"),
                    name="Interaction") +
  scale_x_discrete(labels=c(IPSI="Ipsi\ntargets", CONTRA="Contra\ntargets")) +
  scale_y_continuous(limits=c(NA, 5.5)) +
  facet_wrap(~query, nrow=2) +
  labs(
    x=NULL, y="Z-score",
    title="RNA-RNA interaction: all 9 lncRNAs vs ipsi/contra targets",
    subtitle=paste0(
      "Shared targets only (n=", n_shared, ")  |  ",
      "* p<0.05  ** p<0.01  *** p<0.001 (Wilcoxon, ipsi > contra)\n",
      "New lncRNAs screened with 200nt 5' query — full-length analysis pending IntaRNA"
    )
  ) +
  theme_bw(base_size=9) +
  theme(
    panel.grid.major.x = element_blank(),
    strip.background   = element_rect(fill="grey93"),
    strip.text         = element_text(face="bold", size=8),
    legend.position    = "bottom"
  )

ggsave("outputs/rna_interactions/figures/FigNew_boxplot_9lncrnas.pdf",
       pC, width=16, height=8)
ggsave("outputs/rna_interactions/figures/FigNew_boxplot_9lncrnas.svg",
       pC, width=16, height=8)
cat("Saved: FigNew_boxplot_9lncrnas.pdf\n")

cat("\n=== DONE ===\n")
cat("New figures:\n")
cat("  FigNew_heatmap_all9.pdf\n")
cat("  FigNew_contra_preference_9lncrnas.pdf\n")
cat("  FigNew_boxplot_9lncrnas.pdf\n")
cat(sprintf("\nNote: new lncRNAs compared on %d shared targets only\n", n_shared))
cat("Full analysis pending IntaRNA installation\n")
