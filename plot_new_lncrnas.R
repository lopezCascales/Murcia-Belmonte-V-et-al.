######################################################################
##  Figures for all 9 lncRNAs — original 5 + 4 new
##  Uses zscore_new_lncrnas_axonal_results.csv
##  Run: source("outputs/plot_new_lncrnas.R")
######################################################################

lopez_cascales@BM-581122AEE16D:~/LABMEMBERS/MAYTE$ python3 rna_zscore_new_lncrnas_axonal
python3: can't open file '/home/lopez_cascales/LABMEMBERS/MAYTE/rna_zscore_new_lncrnas_axonal': [Errno 2] No such file or directory
lopez_cascales@BM-581122AEE16D:~/LABMEMBERS/MAYTE$ python3 rna_zscore_new_lncrnas_axonal.py
=================================================================
RNA-RNA Z-score — AXONAL TARGETS ONLY (no nuclear TFs)
39 targets | N=20 shuffles
=================================================================
  Rpph1: 319 nt — FULL sequence
  Kcnq1ot1: 83436 nt -> 200 nt (pending IntaRNA)
  Mir9-3hg: 37711 nt -> 200 nt (pending IntaRNA)
  Arhgap20os: 43073 nt -> 200 nt (pending IntaRNA)

Timing test: Rpph1 x Plxna1...
  TIMEOUT — Rpph1 too slow even for short targets

─────────────────────────────────────────────────────────────────
Rpph1 (lncAI) [FULL]
─────────────────────────────────────────────────────────────────
    Gda            Z=+0.00  n.s.     [IPSI_heatmap]  (3610s)
  **Tenm2          Z=-2.37  p<0.05   [IPSI_heatmap]  (3613s)
  **Boc            Z=-2.52  p<0.05   [IPSI_new]  (1237s)
  **Sema3f         Z=-2.01  p<0.05   [IPSI_new]  (2727s)
    Plxna1         Z=+0.00  n.s.     [IPSI_new]  (3780s)
    Rgs7bp         Z=+0.00  n.s.     [CONTRA_new]  (3780s)
  **Dll1           Z=-2.04  p<0.05   [CONTRA_heatmap]  (1958s)
  -> 4/39 significant

─────────────────────────────────────────────────────────────────
Kcnq1ot1 (lncAC) [200nt pending IntaRNA]
─────────────────────────────────────────────────────────────────
    Gda            Z=+0.00  n.s.     [IPSI_heatmap]  (3610s)
  **Boc            Z=-5.26  p<0.05   [IPSI_new]  (636s)
    Plxna1         Z=+0.00  n.s.     [IPSI_new]  (3780s)
    Rgs7bp         Z=+0.00  n.s.     [CONTRA_new]  (3780s)
  -> 1/39 significant

─────────────────────────────────────────────────────────────────
Mir9-3hg (lncAC) [200nt pending IntaRNA]
─────────────────────────────────────────────────────────────────
    Gda            Z=+0.00  n.s.     [IPSI_heatmap]  (3608s)
  **Boc            Z=-2.65  p<0.05   [IPSI_new]  (552s)
  **Sema3f         Z=-3.16  p<0.05   [IPSI_new]  (1088s)
    Plxna1         Z=+0.00  n.s.     [IPSI_new]  (3780s)
  **Islr2          Z=-3.17  p<0.05   [CONTRA_new]  (1807s)
  **Epha4          Z=-6.22  p<0.05   [CONTRA_new]  (3149s)
    Rgs7bp         Z=+0.00  n.s.     [CONTRA_new]  (3780s)
  -> 4/39 significant

─────────────────────────────────────────────────────────────────
Arhgap20os (lncAI) [200nt pending IntaRNA]
─────────────────────────────────────────────────────────────────
    Gda            Z=+0.00  n.s.     [IPSI_heatmap]  (3608s)
  **Tenm2          Z=-2.04  p<0.05   [IPSI_heatmap]  (1541s)
  **Boc            Z=-2.38  p<0.05   [IPSI_new]  (563s)
    Plxna1         Z=+0.00  n.s.     [IPSI_new]  (3780s)
  **Islr2          Z=-3.08  p<0.05   [CONTRA_new]  (1866s)
    Rgs7bp         Z=+0.00  n.s.     [CONTRA_new]  (3780s)
  **Dll1           Z=-4.42  p<0.05   [CONTRA_heatmap]  (828s)
  -> 4/39 significant

Saved: outputs/rna_interactions/zscore_new_lncrnas_axonal_results.csv
Total: 7273.0 min

=================================================================
SUMMARY
=================================================================

Rpph1 [FULL]: 4/39 sig | 1 CONTRA | 3 IPSI
  Boc            Z=-2.520 [IPSI_new]
  Tenm2          Z=-2.371 [IPSI_heatmap]
  Dll1           Z=-2.036 [CONTRA_heatmap]
  Sema3f         Z=-2.014 [IPSI_new]

Kcnq1ot1 [200nt*]: 1/39 sig | 0 CONTRA | 1 IPSI
  Boc            Z=-5.256 [IPSI_new]

Mir9-3hg [200nt*]: 4/39 sig | 2 CONTRA | 2 IPSI
  Epha4          Z=-6.219 [CONTRA_new]
  Islr2          Z=-3.173 [CONTRA_new]
  Sema3f         Z=-3.156 [IPSI_new]
  Boc            Z=-2.651 [IPSI_new]

Arhgap20os [200nt*]: 4/39 sig | 2 CONTRA | 2 IPSI
  Dll1           Z=-4.425 [CONTRA_heatmap]
  Islr2          Z=-3.077 [CONTRA_new]
  Boc            Z=-2.377 [IPSI_new]
  Tenm2          Z=-2.044 [IPSI_heatmap]

* 200nt query = screening only; full analysis requires IntaRNA

Next: source('outputs/plot_new_lncrnas.R')


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

lnc_order_ai <- c("H19","Lhx1os","Rpph1","Arhgap20os")
lnc_order_ac <- c("Mirg","Meg3","Rian","Kcnq1ot1","Mir9-3hg")
lnc_order    <- c(lnc_order_ai, lnc_order_ac)

lnc_class_map <- c(
  H19="lncAI", Lhx1os="lncAI", Rpph1="lncAI", Arhgap20os="lncAI",
  Mirg="lncAC", Meg3="lncAC", Rian="lncAC",
  Kcnq1ot1="lncAC", `Mir9-3hg`="lncAC"
)

full_analysis <- c("H19","Lhx1os","Mirg","Meg3","Rian","Rpph1")

reclassify <- tribble(
  ~target,   ~new_group,
  "Plxna1",  "CONTRA",
  "Unc5c",   "CONTRA",
  "Dll1",    "CONTRA",
  "Cxcr4",   "SHARED",
  "Igf2",    "CIS_LOCUS"
)

## ── Load data ────────────────────────────────────────────────────────
orig <- read.csv("outputs/rna_interactions/zscore_full_results.csv",
                 stringsAsFactors=FALSE)

new_csv <- "outputs/rna_interactions/zscore_new_lncrnas_axonal_results.csv"
if (!file.exists(new_csv)) stop("Run rna_zscore_new_lncrnas_axonal.py first")
res_new <- read.csv(new_csv, stringsAsFactors=FALSE)

## ── Shared targets (only those tested in new analysis) ───────────────
new_targets <- unique(res_new$target)
tested_new  <- new_targets[!new_targets %in% c("")]

## Filter orig to shared targets
orig_shared <- dplyr::filter(orig, target %in% tested_new)

## Combine
res_all <- dplyr::bind_rows(orig_shared, res_new) %>%
  dplyr::left_join(reclassify, by="target") %>%
  dplyr::mutate(
    group = dplyr::coalesce(new_group, group),
    mol_type = dplyr::case_when(
      group %in% c("IPSI_heatmap","IPSI_new","IPSI") ~ "AXONAL_IPSI",
      group %in% c("CONTRA_heatmap","CONTRA_new","CONTRA") ~ "AXONAL_CONTRA",
      group == "CIS_LOCUS" ~ "CIS_LOCUS",
      group == "SHARED"    ~ "SHARED",
      TRUE ~ "OTHER"
    ),
    specific  = zscore < -2.0,
    lnc_class = dplyr::recode(query, !!!lnc_class_map),
    is_full   = query %in% full_analysis,
    query     = factor(query, levels=lnc_order)
  ) %>%
  dplyr::select(-new_group) %>%
  dplyr::filter(!is.na(query))

res_axon <- dplyr::filter(res_all,
                           mol_type %in% c("AXONAL_IPSI","AXONAL_CONTRA"))

n_targets <- length(unique(res_axon$target))
cat("Shared axonal targets:", n_targets, "\n")
cat("lncRNAs:", paste(levels(res_all$query), collapse=", "), "\n\n")

## Print sig hits
cat("=== Significant hits (Z < -2) ===\n")
res_axon %>%
  dplyr::filter(specific) %>%
  dplyr::arrange(query, zscore) %>%
  dplyr::select(query, target, group, zscore, pval) %>%
  print()
=== Significant hits (Z < -2) ===
  query target          group    zscore   pval
1         H19   Ago3 CONTRA_heatmap -7.452621 p<0.05
2         H19   Smg1 CONTRA_heatmap -6.987588 p<0.05
3         H19 Rgs7bp     CONTRA_new -4.550226 p<0.05
4         H19 Grin2b   IPSI_heatmap -4.185124 p<0.05
5         H19 Sema5a     CONTRA_new -3.518047 p<0.05
6         H19 Plxna4     CONTRA_new -3.238867 p<0.05
7         H19   Ago2 CONTRA_heatmap -3.147020 p<0.05
8         H19  Nova1 CONTRA_heatmap -2.801071 p<0.05
9         H19  Tenm3       IPSI_new -2.179523 p<0.05
10        H19  Islr2     CONTRA_new -2.033321 p<0.05
11     Lhx1os Plxna4     CONTRA_new -3.041892 p<0.05
12     Lhx1os   Dll1         CONTRA -2.075479 p<0.05
13      Rpph1    Boc       IPSI_new -2.520000 p<0.05
14      Rpph1  Tenm2   IPSI_heatmap -2.371000 p<0.05
15      Rpph1   Dll1         CONTRA -2.036000 p<0.05
16      Rpph1 Sema3f       IPSI_new -2.014000 p<0.05
17 Arhgap20os   Dll1         CONTRA -4.425000 p<0.05
18 Arhgap20os  Islr2     CONTRA_new -3.077000 p<0.05
19 Arhgap20os    Boc       IPSI_new -2.377000 p<0.05
20 Arhgap20os  Tenm2   IPSI_heatmap -2.044000 p<0.05
21       Mirg Plxna1         CONTRA -5.386019 p<0.05
22       Mirg  Robo2 CONTRA_heatmap -3.442232 p<0.05
23       Mirg  Epha4     CONTRA_new -2.860030 p<0.05
24       Mirg Rgs7bp     CONTRA_new -2.786733 p<0.05
25       Mirg Grin2b   IPSI_heatmap -2.532797 p<0.05
26       Mirg     Qk CONTRA_heatmap -2.448932 p<0.05
27       Mirg  Robo3   IPSI_heatmap -2.160656 p<0.05
28       Meg3 Rgs7bp     CONTRA_new -6.739728 p<0.05
29       Meg3   Ago2 CONTRA_heatmap -3.354120 p<0.05
30       Meg3   Ago3 CONTRA_heatmap -2.723905 p<0.05
31       Meg3  Robo3   IPSI_heatmap -2.577250 p<0.05
32       Meg3  Epha4     CONTRA_new -2.319895 p<0.05
33       Meg3  Nova1 CONTRA_heatmap -2.232671 p<0.05
34       Meg3     Qk CONTRA_heatmap -2.024260 p<0.05
35       Rian  Tenm2   IPSI_heatmap -6.043745 p<0.05
36       Rian   Smg1 CONTRA_heatmap -3.317335 p<0.05
37       Rian Sema3f       IPSI_new -3.032335 p<0.05
38       Rian Grin2b   IPSI_heatmap -2.982716 p<0.05
39       Rian   Ago3 CONTRA_heatmap -2.660647 p<0.05
40       Rian  Tenm3       IPSI_new -2.216683 p<0.05
41   Kcnq1ot1    Boc       IPSI_new -5.256000 p<0.05
42   Mir9-3hg  Epha4     CONTRA_new -6.219000 p<0.05
43   Mir9-3hg  Islr2     CONTRA_new -3.173000 p<0.05
44   Mir9-3hg Sema3f       IPSI_new -3.156000 p<0.05
45   Mir9-3hg    Boc       IPSI_new -2.651000 p<0.05

## ════════════════════════════════════════════════════════════════════
## FIG A: Heatmap all 9 lncRNAs
## ════════════════════════════════════════════════════════════════════
mat <- res_axon %>%
  dplyr::select(query, target, zscore) %>%
  dplyr::group_by(query, target) %>%
  dplyr::summarise(zscore=mean(zscore), .groups="drop") %>%
  tidyr::pivot_wider(names_from=target, values_from=zscore) %>%
  tibble::column_to_rownames("query") %>%
  as.matrix()

mat_cl <- pmax(mat, -5)

col_ann <- res_axon %>%
  dplyr::distinct(target, mol_type) %>%
  dplyr::mutate(group=ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA")) %>%
  dplyr::select(target, group) %>%
  tibble::column_to_rownames("target")

row_ann <- data.frame(
  lnc_class   = lnc_class_map[rownames(mat)],
  full_analysis = ifelse(rownames(mat) %in% full_analysis, "Full","200nt*"),
  row.names = rownames(mat)
)

sig_mat <- ifelse(mat < -2, "*", "")

hm_col <- colorRampPalette(
  c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#B2182B"))(100)

ann_col <- list(
  group         = c(IPSI=col_ipsi, CONTRA=col_contra),
  lnc_class     = c(lncAI=col_ipsi, lncAC=col_contra),
  full_analysis = c(Full="grey20", `200nt*`="grey75")
)

pheatmap::pheatmap(
  mat_cl,
  color             = hm_col,
  breaks            = seq(-5, 3, length.out=101),
  annotation_col    = col_ann,
  annotation_row    = row_ann,
  annotation_colors = ann_col,
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  border_color      = NA,
  fontsize_row      = 9,
  fontsize_col      = 8,
  display_numbers   = sig_mat,
  number_color      = "black",
  fontsize_number   = 9,
  gaps_row          = 4,   ## gap between lncAI and lncAC
  main = sprintf(
    "RNA-RNA interaction Z-score — 9 lncRNAs x %d axonal targets (* p<0.05, clamped -5)\n* = full 3'UTR analysis; grey = 200nt 5' screening",
    n_targets
  ),
  filename = "outputs/rna_interactions/figures/FigNew_heatmap_all9.pdf",
  width    = max(10, n_targets * 0.55 + 3),
  height   = 7
)
cat("Saved: FigNew_heatmap_all9.pdf\n")

## ════════════════════════════════════════════════════════════════════
## FIG B: CONTRA preference — all 9 lncRNAs
## ════════════════════════════════════════════════════════════════════
sum9 <- res_axon %>%
  dplyr::mutate(group_plot=ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA")) %>%
  dplyr::group_by(query, lnc_class, is_full, group_plot) %>%
  dplyr::summarise(
    mean_z = mean(zscore, na.rm=TRUE),
    se_z   = sd(zscore, na.rm=TRUE)/sqrt(dplyr::n()),
    n_sig  = sum(specific, na.rm=TRUE),
    n      = dplyr::n(),
    .groups="drop"
  )

## Wilcoxon per lncRNA
wt9 <- res_axon %>%
  dplyr::mutate(group_plot=ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA")) %>%
  dplyr::group_by(query, lnc_class, is_full) %>%
  dplyr::summarise(
    p_val = tryCatch(
      wilcox.test(zscore[group_plot=="IPSI"],
                  zscore[group_plot=="CONTRA"],
                  alternative="greater")$p.value,
      error=function(e) NA_real_),
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

cat("\nWilcoxon (ipsi Z > contra Z):\n")
print(wt9 %>% dplyr::select(query, lnc_class, is_full, p_val, sig_label))

pB <- ggplot(sum9, aes(x=group_plot, y=mean_z, colour=query, group=query)) +
  geom_hline(yintercept=0, colour="grey80", linewidth=0.3) +
  geom_hline(yintercept=-2, colour="grey40", linewidth=0.5, linetype="dashed") +
  geom_errorbar(aes(ymin=mean_z-se_z, ymax=mean_z+se_z),
                width=0.12, alpha=0.5, linewidth=0.8) +
  geom_line(aes(linetype=is_full), alpha=0.8, linewidth=1.0) +
  geom_point(aes(size=n_sig, shape=is_full), alpha=0.95) +
  ggrepel::geom_text_repel(
    data=dplyr::filter(sum9, mean_z < -0.8 & group_plot=="CONTRA"),
    aes(label=paste0(query," (",n_sig,")")),
    size=2.5, nudge_x=0.2, box.padding=0.3, max.overlaps=12
  ) +
  scale_colour_manual(values=lnc_col, name="lncRNA") +
  scale_linetype_manual(values=c("TRUE"="solid","FALSE"="dashed"),
                        labels=c("TRUE"="Full analysis","FALSE"="200nt screening"),
                        name="Analysis") +
  scale_shape_manual(values=c("TRUE"=16,"FALSE"=1),
                     labels=c("TRUE"="Full","FALSE"="200nt*"),
                     name="Analysis") +
  scale_size_continuous(range=c(2,7), name="n significant") +
  scale_x_discrete(labels=c(IPSI="Ipsi targets", CONTRA="Contra targets")) +
  facet_wrap(~lnc_class, nrow=1,
    labeller=labeller(lnc_class=c(
      lncAI="lncAI (ipsi axons)\nH19, Lhx1os, Rpph1, Arhgap20os",
      lncAC="lncAC (contra axons)\nMirg, Meg3, Rian, Kcnq1ot1, Mir9-3hg"
    ))
  ) +
  labs(x=NULL, y="Mean Z-score (+/- SE)",
       title="CONTRA preference: all 9 axonal lncRNAs",
       subtitle=paste0(
         "Solid lines = full 3'UTR analysis | Dashed = 200nt 5' screening\n",
         n_targets, " shared axonal targets (no nuclear TFs)"
       )) +
  theme_bw(base_size=10) +
  theme(panel.grid.major.x=element_blank(),
        strip.background=element_rect(fill="grey93"),
        strip.text=element_text(face="bold",size=8),
        legend.position="right")

ggsave("outputs/rna_interactions/figures/FigNew_contra_preference_9lncrnas.pdf",
       pB, width=13, height=5)
ggsave("outputs/rna_interactions/figures/FigNew_contra_preference_9lncrnas.svg",
       pB, width=13, height=5)
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
  geom_jitter(aes(colour=group_plot, size=specific,
                  alpha=is_full),
              width=0.16) +
  geom_text(data=wt9,
            aes(x=1.5, y=Inf, label=sig_label),
            inherit.aes=FALSE, vjust=1.5, size=4.5, colour="grey25") +
  geom_text(data=wt9,
            aes(x=1.5, y=Inf,
                label=ifelse(!is.na(p_val),sprintf("p=%.3f",p_val),"")),
            inherit.aes=FALSE, vjust=3.2, size=2.5, colour="grey45") +
  ggrepel::geom_text_repel(
    data=dplyr::filter(res_box, specific, is_full),
    aes(label=target, colour=group_plot),
    size=2.2, box.padding=0.2, max.overlaps=6, show.legend=FALSE
  ) +
  ## Mark screened (non-full) points differently
  geom_point(data=dplyr::filter(res_box, specific, !is_full),
             aes(colour=group_plot),
             shape=1, size=4, stroke=1.2, alpha=0.7,
             show.legend=FALSE) +
  scale_fill_manual(values=c(IPSI=col_ipsi, CONTRA=col_contra), guide="none") +
  scale_colour_manual(values=c(IPSI=col_ipsi, CONTRA=col_contra), guide="none") +
  scale_size_manual(values=c("FALSE"=1.5,"TRUE"=3.5),
                    labels=c("FALSE"="n.s.","TRUE"="p<0.05"),
                    name="Interaction") +
  scale_alpha_manual(values=c("TRUE"=0.85,"FALSE"=0.45),
                     labels=c("TRUE"="Full analysis","FALSE"="200nt screening"),
                     name="Analysis") +
  scale_x_discrete(labels=c(IPSI="Ipsi\ntargets", CONTRA="Contra\ntargets")) +
  scale_y_continuous(limits=c(NA, 5.5)) +
  facet_wrap(~query, nrow=2,
             labeller=labeller(query=function(x) {
               ifelse(x %in% full_analysis, x, paste0(x,"\n(200nt*)"))
             })) +
  labs(x=NULL, y="Z-score",
       title="RNA-RNA interaction: all 9 lncRNAs vs axonal targets",
       subtitle=paste0(
         n_targets, " shared axonal targets (nuclear TFs excluded)  |  ",
         "* p<0.05  ** p<0.01  *** p<0.001\n",
         "Solid = full analysis | Transparent/open circles = 200nt 5' screening"
       )) +
  theme_bw(base_size=9) +
  theme(panel.grid.major.x=element_blank(),
        strip.background=element_rect(fill="grey93"),
        strip.text=element_text(face="bold",size=8),
        legend.position="bottom")

ggsave("outputs/rna_interactions/figures/FigNew_boxplot_9lncrnas.pdf",
       pC, width=16, height=8)
ggsave("outputs/rna_interactions/figures/FigNew_boxplot_9lncrnas.svg",
       pC, width=16, height=8)
cat("Saved: FigNew_boxplot_9lncrnas.pdf\n")

## ════════════════════════════════════════════════════════════════════
## FIG D: Summary barplot — n significant per lncRNA
## ════════════════════════════════════════════════════════════════════
sig_summary <- res_axon %>%
  dplyr::mutate(group_plot=ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA")) %>%
  dplyr::filter(specific) %>%
  dplyr::count(query, lnc_class, is_full, group_plot) %>%
  dplyr::mutate(
    query = factor(query, levels=lnc_order),
    group_plot = factor(group_plot, levels=c("CONTRA","IPSI"))
  )

pD <- ggplot(sig_summary,
             aes(x=query, y=n, fill=group_plot, alpha=is_full)) +
  geom_col(position="stack", colour="grey30", linewidth=0.4) +
  geom_text(aes(label=n), position=position_stack(vjust=0.5),
            size=3, colour="white", fontface="bold") +
  scale_fill_manual(values=c(IPSI=col_ipsi, CONTRA=col_contra),
                    name="Target type") +
  scale_alpha_manual(values=c("TRUE"=0.95,"FALSE"=0.55),
                     labels=c("TRUE"="Full analysis","FALSE"="200nt screening"),
                     name="Analysis") +
  scale_x_discrete(labels=function(x)
    ifelse(x %in% full_analysis, x, paste0(x,"\n(200nt*)"))) +
  geom_vline(xintercept=4.5, linetype="dashed", colour="grey50") +
  annotate("text", x=2.5, y=Inf, vjust=1.5, label="lncAI (ipsi)",
           colour=col_ipsi, fontface="bold", size=3.5) +
  annotate("text", x=7, y=Inf, vjust=1.5, label="lncAC (contra)",
           colour=col_contra, fontface="bold", size=3.5) +
  labs(x=NULL, y="Number of significant interactions (Z < -2)",
       title="Significant RNA-RNA interactions per lncRNA",
       subtitle=paste0(
         "Filled bars = full analysis | Pale bars = 200nt 5' screening\n",
         "IPSI (blue) vs CONTRA (red) axonal targets"
       )) +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(size=8),
        legend.position="right")

ggsave("outputs/rna_interactions/figures/FigNew_sig_barplot.pdf",
       pD, width=11, height=5)
ggsave("outputs/rna_interactions/figures/FigNew_sig_barplot.svg",
       pD, width=11, height=5)
cat("Saved: FigNew_sig_barplot.pdf\n")

cat("\n=== DONE ===\n")
cat("Figures:\n")
cat("  FigNew_heatmap_all9.pdf\n")
cat("  FigNew_contra_preference_9lncrnas.pdf\n")
cat("  FigNew_boxplot_9lncrnas.pdf\n")
cat("  FigNew_sig_barplot.pdf\n")
