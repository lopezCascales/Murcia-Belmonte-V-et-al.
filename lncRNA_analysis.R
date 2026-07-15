R

R version 4.3.3 (2024-02-29) -- "Angel Food Cake"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> 
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
── Attaching core tidyverse packages ────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.2.1     ✔ readr     2.2.0
✔ forcats   1.0.1     ✔ stringr   1.6.0
✔ ggplot2   4.0.2     ✔ tibble    3.3.1
✔ lubridate 1.9.5     ✔ tidyr     1.3.2
✔ purrr     1.2.1     
── Conflicts ──────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
ℹ Use the conflicted package to force all conflicts to become errors
[conflicted] Will prefer dplyr::filter over any other package.
[conflicted] Will prefer dplyr::select over any other package.
[conflicted] Will prefer dplyr::rename over any other package.
> lnc_class_map <- c(
  lnc_class_map <- c(="lncAI", Rpph1="lncAI", Arhgap20os="lncAI",
  H19="lncAI", Lhx1os="lncAI", Rpph1="lncAI", Arhgap20os="lncAI",r9-3hg`="lncAC"
  Mirg="lncAC", Meg3="lncAC", Rian="lncAC", Kcnq1ot1="lncAC", `Mir9-3hg`="lncAC"
)
reclassify <- tribble(
reclassify <- tribble(p,
  ~target,   ~new_group,  "Unc5c",  "CONTRA",
  "Lhx2",    "IPSI_TF",   "Unc5c",  "CONTRA",F",
  "Plxna1",  "CONTRA",    "Pou4f1", "CONTRA_TF",
  "Sox4",    "CONTRA_TF", "Sox11",  "CONTRA_TF",
  "Isl2",    "CONTRA_TF", "Cxcr4",  "SHARED",
  "Igf2",    "CIS_LOCUS"
)
NUCLEAR <- c("Neurod2","Fezf1","Lhx1","Zic4","Nr2f2","Atoh7",
NUCLEAR <- c("Neurod2","Fezf1","Lhx1","Zic4","Nr2f2","Atoh7",",
             "Bcl11a","Bcl11b","Kmt2a","Kmt2b","Hdac4","Ep300",
             "Pou4f2","Sox4","Sox11","Isl2","Pou4f1","Foxd1",
             "Zic2","Lhx2","Nr2f1")"Qk","Smg1","Tnrc6a","Tnrc6b","Cnot3")
RBP_CPE <- c("Ago2","Ago3","Nova1","Qk","Smg1","Tnrc6a","Tnrc6b","Cnot3")
## ── Load original 5 lncRNAs ──────────────────────────────────────────
## ── Load original 5 lncRNAs ──────────────────────────────────────────
orig_csv <- "outputs/rna_interactions/zscore_full_results.csv"")
if (!file.exists(orig_csv)) stop("Run rna_zscore_full.py first")
res_orig <- read.csv(orig_csv, stringsAsFactors=FALSE)
cat("Original results:", nrow(res_orig), "pairs\n")
## ── Load new 4 lncRNAs ───────────────────────────────────────────────
## ── Load new 4 lncRNAs ───────────────────────────────────────────────
new_csv <- "outputs/rna_interactions/zscore_new_lncrnas_results.csv"")
if (!file.exists(new_csv)) stop("Run rna_zscore_new_rnaplex.py first")
res_new <- read.csv(new_csv, stringsAsFactors=FALSE)
cat("New results:", nrow(res_new), "pairs\n")y), collapse=", "), "\n")
cat("New lncRNAs:", paste(unique(res_new$query), collapse=", "), "\n")
Original results: 320 pairs
New results: 40 pairs
New lncRNAs: Kcnq1ot1, Mir9-3hg, Rpph1, Arhgap20os 
> shared_targets <- intersect(
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
Shared targets: 10 -> Gda, Tenm2, Grin2b, Ago2, Ago3, Nova1, Smg1, Plxna1, Plxna4, Rgs7bp 
Error in `dplyr::mutate()`:
ℹ In argument: `query = factor(query, levels = lnc_order)`.
Caused by error:
! object 'lnc_order' not found
Run `rlang::last_trace()` to see where the error occurred.
> orig_csv <- "outputs/rna_interactions/zscore_full_results.csv"
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
} cat("(200nt query screening — expected for highly truncated sequences)\n")20os"),
Original results: 320 pairs
New results: 40 pairs
New lncRNAs: Kcnq1ot1, Mir9-3hg, Rpph1, Arhgap20os 
Shared targets: 10 -> Gda, Tenm2, Grin2b, Ago2, Ago3, Nova1, Smg1, Plxna1, Plxna4, Rgs7bp 
Error in `dplyr::mutate()`:
ℹ In argument: `query = factor(query, levels = lnc_order)`.
Caused by error:
! object 'lnc_order' not found
Run `rlang::last_trace()` to see where the error occurred.
> lnc_order <- c("H19","Lhx1os","Rpph1","Arhgap20os",
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
> NUCLEAR <- c("Neurod2","Fezf1","Lhx1","Zic4","Nr2f2","Atoh7",
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
Original results: 320 pairs
New results: 40 pairs
New lncRNAs: Kcnq1ot1, Mir9-3hg, Rpph1, Arhgap20os 
Shared targets: 10 -> Gda, Tenm2, Grin2b, Ago2, Ago3, Nova1, Smg1, Plxna1, Plxna4, Rgs7bp 
> res_orig_shared <- dplyr::filter(res_orig, target %in% shared_targets)
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
} cat("(200nt query screening — expected for highly truncated sequences)\n")

Axonal pairs for shared targets: 54 

=== Significant hits (Z < -2) from NEW lncRNAs ===
No significant interactions found for new lncRNAs
(200nt query screening — expected for highly truncated sequences)
> mat <- res_axon %>%
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
cat("Saved: FigNew_heatmap_all9.pdf\n"),igures/FigNew_heatmap_all9.pdf",ets (* p<0.05)",

Heatmap: 9 lncRNAs x 6 shared targets
There were 21 warnings (use warnings() to see them)
Saved: FigNew_heatmap_all9.pdf
> sum9 <- res_axon %>%
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
print(wt9 %>% dplyr::select(query, lnc_class, p_val, sig_label))
Warning message:
There were 4 warnings in `dplyr::summarise()`.
The first warning was:
ℹ In argument: `p_val = tryCatch(...)`.
ℹ In group 3: `query = Rpph1`, `lnc_class = "lncAI"`.
Caused by warning in `wilcox.test.default()`:
! cannot compute exact p-value with ties
ℹ Run `dplyr::last_dplyr_warnings()` to see the 3 remaining warnings. 

Wilcoxon results (all 9 lncRNAs):
# A tibble: 9 × 4
  query      lnc_class p_val sig_label
  <fct>      <chr>     <dbl> <chr>    
1 H19        lncAI     0.2   n.s.     
2 Lhx1os     lncAI     0.1   n.s.     
3 Rpph1      lncAI     0.252 n.s.     
4 Arhgap20os lncAI     1     n.s.     
5 Mirg       lncAC     0.2   n.s.     
6 Meg3       lncAC     0.5   n.s.     
7 Rian       lncAC     0.9   n.s.     
8 Kcnq1ot1   lncAC     1     n.s.     
9 Mir9-3hg   lncAC     0.909 n.s.     
> pB <- ggplot(sum9, aes(x=group_plot, y=mean_z, colour=query, group=query)) +
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
cat("Saved: FigNew_contra_preference_9lncrnas.pdf\n")a_preference_9lncrnas.svg",
Saved: FigNew_contra_preference_9lncrnas.pdf
> res_box <- res_axon %>%
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
cat("Full analysis pending IntaRNA installation\n")red targets only\n", n_shared))aRNA"
There were 11 warnings (use warnings() to see them)
Warning messages:
1: Removed 1 row containing non-finite outside the scale range (`stat_boxplot()`). 
2: Removed 1 row containing missing values or values outside the scale range (`geom_point()`). 
Saved: FigNew_boxplot_9lncrnas.pdf

=== DONE ===
New figures:
  FigNew_heatmap_all9.pdf
  FigNew_contra_preference_9lncrnas.pdf
  FigNew_boxplot_9lncrnas.pdf

Note: new lncRNAs compared on 6 shared targets only
Full analysis pending IntaRNA installation
> 
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
[conflicted] Removing existing preference.
[conflicted] Will prefer dplyr::filter over any other package.
[conflicted] Removing existing preference.
[conflicted] Will prefer dplyr::select over any other package.
[conflicted] Removing existing preference.
[conflicted] Will prefer dplyr::rename over any other package.
> orig <- read.csv("outputs/rna_interactions/zscore_full_results.csv",
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
Shared axonal targets: 39 
lncRNAs: H19, Lhx1os, Rpph1, Arhgap20os, Mirg, Meg3, Rian, Kcnq1ot1, Mir9-3hg 

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
> mat <- res_axon %>%
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
cat("Saved: FigNew_heatmap_all9.pdf\n") 3),res/FigNew_heatmap_all9.pdf",
There were 21 warnings (use warnings() to see them)
Saved: FigNew_heatmap_all9.pdf
> sum9 <- res_axon %>%
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
Warning message:
There were 4 warnings in `dplyr::summarise()`.
The first warning was:
ℹ In argument: `p_val = tryCatch(...)`.
ℹ In group 3: `query = Rpph1`, `lnc_class = "lncAI"`, `is_full = FALSE`.
Caused by warning in `wilcox.test.default()`:
! cannot compute exact p-value with ties
ℹ Run `dplyr::last_dplyr_warnings()` to see the 3 remaining warnings. 

Wilcoxon (ipsi Z > contra Z):
# A tibble: 9 × 5
  query      lnc_class is_full  p_val sig_label
  <fct>      <chr>     <lgl>    <dbl> <chr>    
1 H19        lncAI     FALSE   0.0187 *        
2 Lhx1os     lncAI     FALSE   0.157  n.s.     
3 Rpph1      lncAI     FALSE   0.947  n.s.     
4 Arhgap20os lncAI     FALSE   0.916  n.s.     
5 Mirg       lncAC     FALSE   0.114  n.s.     
6 Meg3       lncAC     FALSE   0.0569 n.s.     
7 Rian       lncAC     FALSE   0.584  n.s.     
8 Kcnq1ot1   lncAC     FALSE   0.750  n.s.     
9 Mir9-3hg   lncAC     FALSE   0.788  n.s.     
> pB <- ggplot(sum9, aes(x=group_plot, y=mean_z, colour=query, group=query)) +
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
Saved: FigNew_contra_preference_9lncrnas.pdf
> res_box <- res_axon %>%
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
cat("Saved: FigNew_boxplot_9lncrnas.pdf\n")gNew_boxplot_9lncrnas.svg",
Warning messages:
1: Removed 37 rows containing non-finite outside the scale range (`stat_boxplot()`). 
2: Removed 37 rows containing missing values or values outside the scale range (`geom_point()`). 
Warning messages:
1: Removed 37 rows containing non-finite outside the scale range (`stat_boxplot()`). 
2: Removed 37 rows containing missing values or values outside the scale range (`geom_point()`). 
Saved: FigNew_boxplot_9lncrnas.pdf
> sig_summary <- res_axon %>%
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
Saved: FigNew_sig_barplot.pdf
> cat("\n=== DONE ===\n")
cat("Figures:\n")
cat("  FigNew_heatmap_all9.pdf\n")
cat("  FigNew_contra_preference_9lncrnas.pdf\n")
cat("  FigNew_boxplot_9lncrnas.pdf\n")
cat("  FigNew_sig_barplot.pdf\n")

=== DONE ===
Figures:
  FigNew_heatmap_all9.pdf
  FigNew_contra_preference_9lncrnas.pdf
  FigNew_boxplot_9lncrnas.pdf
  FigNew_sig_barplot.pdf
> 
