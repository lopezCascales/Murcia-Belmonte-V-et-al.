######################################################################
##  FIGURES v4 — corrected
##  - Fixed pheatmap annotation error (Category levels)
##  - Plxna1 reclassified to CONTRA (Kuwajima 2012)
##  - Dinucleotide results integrated
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
col_rbp    <- "#7B3F9E"
col_cis    <- "#EF9F27"
col_shared <- "#9B59B6"
lnc_col    <- c(H19="#1D9E75", Lhx1os="#9FE1CB",
                Mirg="#C0395A", Meg3="#D4537E", Rian="#E89BB0")

NUCLEAR_TF  <- c("Neurod2","Fezf1","Lhx1","Zic4","Nr2f2","Atoh7")
NUCLEAR_EPI <- c("Bcl11a","Bcl11b","Kmt2a","Kmt2b","Hdac4","Ep300")
RBP_CPE     <- c("Ago2","Ago3","Nova1","Qk","Smg1","Tnrc6a","Tnrc6b","Cnot3")
INTRACELL   <- c("Ctnnb1","Apc2")

tf_lens <- c(Zic2=611,Pou4f2=1729,Lhx2=1014,Nr2f1=1197,Foxd1=3317,
             Isl2=757,Sox4=2782,Sox11=6959,Pou4f1=2868,
             Neurod2=1773,Fezf1=696,Lhx1=1360,Zic4=2813,
             Nr2f2=2716,Atoh7=677,Bcl11a=3266,Bcl11b=5000,
             Kmt2a=4550,Kmt2b=306,Hdac4=4571,Ep300=1232)

tf_class_map <- c(
  Zic2="IPSI_TF", Lhx2="IPSI_TF", Nr2f1="IPSI_TF",
  Lhx1="IPSI_TF", Zic4="IPSI_TF", Nr2f2="IPSI_TF",
  Pou4f2="CONTRA_TF", Isl2="CONTRA_TF", Sox4="CONTRA_TF",
  Sox11="CONTRA_TF", Pou4f1="CONTRA_TF", Foxd1="neutral_TF",
  Neurod2="CONTRA_TF", Fezf1="CONTRA_TF", Atoh7="CONTRA_TF",
  Bcl11a="NUCLEAR_EPI", Bcl11b="NUCLEAR_EPI",
  Kmt2a="NUCLEAR_EPI", Kmt2b="NUCLEAR_EPI",
  Hdac4="NUCLEAR_EPI", Ep300="NUCLEAR_EPI"
)

## ── Load + reclassify ────────────────────────────────────────────────
res_raw <- read.csv("outputs/rna_interactions/zscore_full_results.csv",
                    stringsAsFactors=FALSE)

reclassify <- tribble(
  ~target,   ~new_group,
  "Lhx2",    "IPSI_TF",
  "Unc5c",   "CONTRA",
  "Plxna1",  "CONTRA",      ## Kuwajima 2012: PlexinA1 on cRGC axons
  "Pou4f1",  "CONTRA_TF",
  "Sox4",    "CONTRA_TF",
  "Sox11",   "CONTRA_TF",
  "Isl2",    "CONTRA_TF",
  "Cxcr4",   "SHARED",
  "Igf2",    "CIS_LOCUS"
)

res <- res_raw %>%
  dplyr::left_join(reclassify, by="target") %>%
  dplyr::mutate(
    group = dplyr::coalesce(new_group, group),
    mol_type = dplyr::case_when(
      target %in% NUCLEAR_TF  ~ "NUCLEAR_TF",
      target %in% NUCLEAR_EPI ~ "NUCLEAR_EPI",
      target %in% INTRACELL   ~ "INTRACELLULAR",
      target %in% RBP_CPE     ~ "RBP_CPE",
      group %in% c("NEG_CTRL_TF","IPSI_TF","CONTRA_TF") ~ "NUCLEAR_TF",
      group == "CIS_LOCUS"    ~ "CIS_LOCUS",
      group == "SHARED"       ~ "SHARED",
      group %in% c("IPSI_heatmap","IPSI_new","IPSI") ~ "AXONAL_IPSI",
      group %in% c("CONTRA_heatmap","CONTRA_new","CONTRA") ~ "AXONAL_CONTRA",
      TRUE ~ "OTHER"
    ),
    group_simple = dplyr::case_when(
      mol_type == "AXONAL_IPSI"   ~ "IPSI",
      mol_type == "AXONAL_CONTRA" ~ "CONTRA",
      TRUE ~ mol_type
    ),
    specific  = zscore < -2.0,
    lnc_class = ifelse(query %in% c("H19","Lhx1os"), "lncAI", "lncAC"),
    query = factor(query, levels=c("H19","Lhx1os","Mirg","Meg3","Rian"))
  ) %>%
  dplyr::select(-new_group)

## Load dinucleotide results
dinuc_csv <- "outputs/rna_interactions/zscore_dinucleotide_results.csv"
if (file.exists(dinuc_csv)) {
  dinuc <- read.csv(dinuc_csv, stringsAsFactors=FALSE) %>%
    dplyr::mutate(survives = as.logical(survives)) %>%
    dplyr::mutate(pair=paste0(query,"_",target))
  cat("Dinucleotide results loaded:", nrow(dinuc), "pairs\n")
  cat("Surviving dinuc shuffle:", sum(dinuc$survives), "/", nrow(dinuc), "\n")
} else {
  dinuc <- NULL
  cat("No dinucleotide CSV found\n")
}

res_axon <- dplyr::filter(res, mol_type %in% c("AXONAL_IPSI","AXONAL_CONTRA"))
res_rbp  <- dplyr::filter(res, mol_type == "RBP_CPE")
res_nuc  <- dplyr::filter(res, mol_type %in% c("NUCLEAR_TF","NUCLEAR_EPI"))

cat("Nrp2 in axonal data:", "Nrp2" %in% res_axon$target, "\n")
cat("Plxna1 group:", unique(res$mol_type[res$target=="Plxna1"]), "\n")

row_ann <- data.frame(lnc_class=c("lncAI","lncAI","lncAC","lncAC","lncAC"),
                      row.names=c("H19","Lhx1os","Mirg","Meg3","Rian"))
ann_lnc <- list(lnc_class=c(lncAI=col_ipsi, lncAC=col_contra))
hm_col  <- colorRampPalette(c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#B2182B"))(100)

make_mat <- function(data) {
  data %>%
    dplyr::select(query, target, zscore) %>%
    tidyr::pivot_wider(names_from=target, values_from=zscore, values_fn=mean) %>%
    tibble::column_to_rownames("query") %>%
    as.matrix()
}

## ════════════════════════════════════════════════════════════════════
## FIG 1a: Axonal guidance molecules
## ════════════════════════════════════════════════════════════════════
mat_a    <- make_mat(res_axon)
mat_a_cl <- pmax(mat_a, -5)
col_ann_a <- res_axon %>%
  dplyr::mutate(group=ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA")) %>%
  dplyr::distinct(target, group) %>%
  tibble::column_to_rownames("target")
sig_a <- ifelse(mat_a < -2, "*", "")

pheatmap::pheatmap(
  mat_a_cl, color=hm_col, breaks=seq(-5,3,length.out=101),
  annotation_col=col_ann_a, annotation_row=row_ann,
  annotation_colors=c(ann_lnc, list(group=c(IPSI=col_ipsi, CONTRA=col_contra))),
  cluster_rows=FALSE, cluster_cols=TRUE, border_color=NA,
  fontsize_row=10, fontsize_col=7,
  display_numbers=sig_a, number_color="black", fontsize_number=9,
  main="RNA-RNA interaction Z-score -- axonal guidance & adhesion molecules  (* p<0.05, clamped -5)",
  filename="outputs/rna_interactions/figures/Fig1a_heatmap_axonal.pdf",
  width=max(12,ncol(mat_a)*0.28), height=5
)
cat("Saved: Fig1a_heatmap_axonal.pdf  |  Nrp2:",
    "Nrp2" %in% colnames(mat_a), "\n")


## ════════════════════════════════════════════════════════════════════
## FIG 1b: RBP/CPE translation machinery
## ════════════════════════════════════════════════════════════════════
mat_b <- make_mat(res_rbp)
sig_b <- ifelse(mat_b < -2, "*", "")

rbp_fn <- c(Ago2="RISC/miRNA", Ago3="RISC/miRNA",
            Tnrc6a="RISC/miRNA", Tnrc6b="RISC/miRNA",
            Nova1="Splicing/CPE", Qk="Splicing/CPE",
            Smg1="NMD", Cnot3="Deadenylation")

col_ann_b <- data.frame(
  Function = rbp_fn[colnames(mat_b)],
  row.names = colnames(mat_b),
  stringsAsFactors = FALSE
)
## Fill any missing with "other"
col_ann_b$Function[is.na(col_ann_b$Function)] <- "other"

fn_cols <- c("RISC/miRNA"=col_rbp, "Splicing/CPE"="#1D9E75",
             "NMD"=col_cis, "Deadenylation"="grey60", "other"="grey80")

pheatmap::pheatmap(
  mat_b, color=hm_col, breaks=seq(-8,3,length.out=101),
  annotation_col=col_ann_b, annotation_row=row_ann,
  annotation_colors=c(ann_lnc, list(Function=fn_cols)),
  cluster_rows=FALSE, cluster_cols=TRUE, border_color=NA,
  fontsize_row=10, fontsize_col=9,
  display_numbers=sig_b, number_color="black", fontsize_number=9,
  main="RNA-RNA Z-score -- CPE/translation regulatory machinery  (* p<0.05)\nNote: most hits pending dinucleotide shuffle validation",
  filename="outputs/rna_interactions/figures/Fig1b_heatmap_RBP_CPE.pdf",
  width=max(7,ncol(mat_b)*0.55), height=5
)
cat("Saved: Fig1b_heatmap_RBP_CPE.pdf\n")


## ════════════════════════════════════════════════════════════════════
## FIG 1c: Nuclear TF + epigenetic — FIXED annotation error
## ════════════════════════════════════════════════════════════════════
mat_c <- make_mat(res_nuc)
sig_c <- ifelse(mat_c < -2, "*", "")

## Build annotation carefully — only genes present in matrix
genes_in_mat <- colnames(mat_c)
cat_vec <- dplyr::case_when(
  genes_in_mat %in% NUCLEAR_TF  ~ "Nuclear TF",
  genes_in_mat %in% NUCLEAR_EPI ~ "Epigenetic",
  TRUE ~ "other"
)
tf_vec <- dplyr::coalesce(tf_class_map[genes_in_mat], "unknown")

col_ann_c <- data.frame(
  Category = cat_vec,
  TF_class = tf_vec,
  row.names = genes_in_mat,
  stringsAsFactors = FALSE
)

## Colours — only for levels actually present
cat_levels   <- unique(col_ann_c$Category)
tf_levels    <- unique(col_ann_c$TF_class)

cat_col_all  <- c("Nuclear TF"="grey65", "Epigenetic"="grey30", "other"="grey85")
tf_col_all   <- c(IPSI_TF=col_ipsi, CONTRA_TF=col_contra,
                  neutral_TF="grey70", NUCLEAR_EPI="#444441", unknown="grey90")

pheatmap::pheatmap(
  mat_c, color=hm_col, breaks=seq(-15,3,length.out=101),
  annotation_col=col_ann_c, annotation_row=row_ann,
  annotation_colors=c(ann_lnc,
    list(Category = cat_col_all[cat_levels],
         TF_class = tf_col_all[tf_levels])),
  cluster_rows=FALSE, cluster_cols=TRUE, border_color=NA,
  fontsize_row=10, fontsize_col=8,
  display_numbers=sig_c, number_color="black", fontsize_number=8,
  main="RNA-RNA Z-score -- Nuclear TFs and epigenetic regulators  (* p<0.05)\nSignificant hits explained by 3'UTR length artifact (see Fig 5)",
  filename="outputs/rna_interactions/figures/Fig1c_heatmap_nuclear.pdf",
  width=max(10,ncol(mat_c)*0.45), height=5
)
cat("Saved: Fig1c_heatmap_nuclear.pdf\n")


## ════════════════════════════════════════════════════════════════════
## FIG 2: Method validation
## ════════════════════════════════════════════════════════════════════
panel_a_data <- data.frame(
  pair    = c("H19 x Gda\n(ipsi, exp.validated,\nno interaction)",
              "H19 x Sema5a\n(CONTRA, clean\ncomposition)",
              "H19 x Islr2\n(CONTRA, clean\ncomposition)",
              "Mirg x Plxna1\n(CONTRA, dinuc\nZ=-3.97)"),
  real_dg = c(-51.6,-68.5,-75.0,-57.9),
  bg_mean = c(-52.8,-60.8,-67.0,-44.9),
  bg_std  = c( 3.4,  2.7,  4.0,  2.4),
  zscore  = c(+0.34,-3.52,-2.03,-5.39),
  group   = c("IPSI (neg ctrl)","CONTRA","CONTRA","CONTRA"),
  stringsAsFactors=FALSE
) %>%
  dplyr::mutate(
    bg_lo=bg_mean-bg_std, bg_hi=bg_mean+bg_std,
    pair=factor(pair,levels=rev(pair))
  )

pA <- ggplot(panel_a_data, aes(y=pair)) +
  geom_errorbarh(aes(xmin=bg_lo,xmax=bg_hi,colour=group),
                 height=0.25,linewidth=1,alpha=0.35) +
  geom_point(aes(x=bg_mean,colour=group),shape=124,size=5,alpha=0.4) +
  geom_point(aes(x=real_dg,colour=group,
                 shape=zscore< -2,size=zscore< -2)) +
  geom_segment(aes(x=real_dg,xend=bg_mean,yend=pair,colour=group),
               alpha=0.25,linewidth=0.5) +
  geom_text(aes(x=real_dg-1.5,
                label=ifelse(zscore< -2,
                             paste0("Z=",round(zscore,2),"*"),
                             paste0("Z=",round(zscore,2)))),
            hjust=1,size=3,colour="grey20") +
  scale_colour_manual(values=c("IPSI (neg ctrl)"="grey55",CONTRA=col_contra)) +
  scale_shape_manual(values=c("FALSE"=16,"TRUE"=18),guide="none") +
  scale_size_manual(values=c("FALSE"=3,"TRUE"=5),guide="none") +
  labs(x="delta G (kcal/mol)",y=NULL,colour=NULL,
       title="A. Co-expression != direct RNA-RNA interaction",
       subtitle="Gda validated in ipsi axons by IHC but Z=+0.34") +
  theme_bw(base_size=10) +
  theme(legend.position="bottom",panel.grid.major.y=element_blank(),
        plot.title=element_text(face="bold"))

h19_axon <- dplyr::filter(res,query=="H19",
                           mol_type %in% c("AXONAL_IPSI","AXONAL_CONTRA")) %>%
  dplyr::mutate(group_plot=ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA"))
wt <- wilcox.test(h19_axon$zscore[h19_axon$mol_type=="AXONAL_IPSI"],
                  h19_axon$zscore[h19_axon$mol_type=="AXONAL_CONTRA"],
                  alternative="greater")

pB <- ggplot(h19_axon,aes(x=group_plot,y=zscore,colour=group_plot)) +
  geom_hline(yintercept=0,colour="grey80",linewidth=0.3) +
  geom_hline(yintercept=-2,colour="grey40",linewidth=0.5,linetype="dashed") +
  geom_boxplot(width=0.35,alpha=0.2,outlier.shape=NA,linewidth=0.6) +
  geom_jitter(aes(size=specific),width=0.12,alpha=0.9) +
  ggrepel::geom_text_repel(
    data=dplyr::filter(h19_axon,
                       specific|target %in% c("Gda","Nrp2","Opcml")),
    aes(label=target),size=2.8,colour="grey15",
    box.padding=0.35,max.overlaps=15) +
  annotate("text",x=1.5,y=max(h19_axon$zscore)+0.5,
           label=sprintf("Wilcoxon p=%.3f",wt$p.value),
           size=3.2,colour="grey25") +
  scale_colour_manual(values=c(IPSI=col_ipsi,CONTRA=col_contra)) +
  scale_size_manual(values=c("FALSE"=2,"TRUE"=4),guide="none") +
  scale_x_discrete(labels=c(IPSI="Ipsi\ntargets",CONTRA="Contra\ntargets")) +
  labs(x=NULL,y="Z-score",colour=NULL,
       title="B. H19 prefers CONTRA axon targets",
       subtitle=sprintf("Wilcoxon p=%.3f (guidance molecules only)",wt$p.value)) +
  theme_bw(base_size=10) +
  theme(legend.position="none",panel.grid.major.x=element_blank(),
        plot.title=element_text(face="bold"))

fig2 <- pA+pB+plot_layout(widths=c(1.5,1))+
  plot_annotation(title="RNA-RNA interaction: method validation",
                  theme=theme(plot.title=element_text(size=12,face="bold")))
ggsave("outputs/rna_interactions/figures/Fig2_method_validation.pdf",
       fig2,width=12,height=5)
ggsave("outputs/rna_interactions/figures/Fig2_method_validation.svg",
       fig2,width=12,height=5)
cat("Saved: Fig2_method_validation.pdf\n")


## ════════════════════════════════════════════════════════════════════
## FIG 3: CONTRA preference
## ════════════════════════════════════════════════════════════════════
res_fig3 <- dplyr::filter(res,
  mol_type %in% c("AXONAL_IPSI","AXONAL_CONTRA","RBP_CPE")) %>%
  dplyr::mutate(
    mol_cat = ifelse(mol_type=="RBP_CPE",
                     "CPE/translation machinery",
                     "Guidance/adhesion molecules"),
    grp = ifelse(mol_type=="AXONAL_IPSI","IPSI","CONTRA")
  )

sum3 <- res_fig3 %>%
  dplyr::group_by(query,lnc_class,grp,mol_cat) %>%
  dplyr::summarise(mean_z=mean(zscore),se_z=sd(zscore)/sqrt(dplyr::n()),
                   n_sig=sum(specific),n=dplyr::n(),.groups="drop")

p3 <- ggplot(sum3,aes(x=grp,y=mean_z,colour=query,group=query)) +
  geom_hline(yintercept=0,colour="grey80",linewidth=0.3) +
  geom_hline(yintercept=-2,colour="grey40",linewidth=0.5,linetype="dashed") +
  geom_errorbar(aes(ymin=mean_z-se_z,ymax=mean_z+se_z),
                width=0.12,alpha=0.5,linewidth=0.7) +
  geom_line(alpha=0.7,linewidth=0.9) +
  geom_point(aes(size=n_sig),alpha=0.95) +
  ggrepel::geom_text_repel(
    data=dplyr::filter(sum3,mean_z< -1.2&grp=="CONTRA"),
    aes(label=paste0(query,"\n(",n_sig,"/",n,")")),
    size=2.3,nudge_x=0.18,box.padding=0.3,max.overlaps=8) +
  scale_colour_manual(values=lnc_col) +
  scale_size_continuous(range=c(2,7),breaks=c(1,4,8,12),
                        name="n significant") +
  scale_x_discrete(labels=c(IPSI="Ipsi\ntargets",CONTRA="Contra\ntargets")) +
  facet_grid(mol_cat~lnc_class,
    labeller=labeller(
      lnc_class=c(lncAI="lncAI (ipsi)",lncAC="lncAC (contra)"),
      mol_cat=label_wrap_gen(25))) +
  labs(x=NULL,y="Mean Z-score (+/-SE)",colour="lncRNA",
       title="Axonal lncRNAs preferentially interact with CONTRA program targets",
       subtitle="Pattern consistent across both guidance molecules and CPE machinery") +
  theme_bw(base_size=10) +
  theme(panel.grid.major.x=element_blank(),
        strip.background=element_rect(fill="grey93"),
        strip.text=element_text(face="bold",size=9))
ggsave("outputs/rna_interactions/figures/Fig3_contra_preference.pdf",
       p3,width=11,height=7)
ggsave("outputs/rna_interactions/figures/Fig3_contra_preference.svg",
       p3,width=11,height=7)
cat("Saved: Fig3_contra_preference.pdf\n")


## ════════════════════════════════════════════════════════════════════
## FIG 4: Robust candidates (dinucleotide shuffle results)
## ════════════════════════════════════════════════════════════════════
if (!is.null(dinuc)) {

  ## Add Plxna1 correction to dinuc data
  dinuc <- dinuc %>%
    dplyr::mutate(
      group = dplyr::case_when(
        target=="Plxna1" ~ "CONTRA",
        TRUE ~ group
      ),
      mol_cat = dplyr::case_when(
        target %in% RBP_CPE     ~ "CPE/translation",
        target %in% NUCLEAR_TF  ~ "Nuclear TF",
        target %in% NUCLEAR_EPI ~ "Epigenetic",
        TRUE ~ "Guidance/adhesion"
      )
    )

  ## Note: H19 x Gda survived dinuc (Z=-2.69) — this is the NEG CTRL
  ## It means H19 DOES interact with Gda 3'UTR under dinuc shuffle
  ## BUT Gda is co-expressed (Spearman) and validated in ipsi axons
  ## Interpret as: possible weak interaction, needs wet-lab validation
  ## Mark differently from confident hits

  dinuc_plot <- dinuc %>%
    dplyr::mutate(
      label = dplyr::case_when(
        target=="Gda" ~ paste0(target,"\n(exp. validated,\nneg ctrl)"),
        survives ~ target,
        TRUE ~ ""
      ),
      alpha_val = ifelse(target=="Gda", 0.4, 0.9),
      border_col = ifelse(target=="Gda","grey40","transparent")
    )

  p4 <- ggplot(dinuc_plot,
               aes(x=mono_z, y=dinuc_z, colour=group)) +
    geom_abline(slope=1,intercept=0,linetype="dashed",
                colour="grey60",linewidth=0.5) +
    geom_hline(yintercept=-2,colour="grey40",linewidth=0.4,linetype="dotted") +
    geom_vline(xintercept=-2,colour="grey40",linewidth=0.4,linetype="dotted") +
    geom_point(aes(size=survives,alpha=alpha_val,shape=mol_cat)) +
    ggrepel::geom_text_repel(
      data=dplyr::filter(dinuc_plot,survives|target=="Gda"),
      aes(label=target),size=2.8,colour="grey15",
      box.padding=0.35,max.overlaps=20) +
    scale_colour_manual(values=c(IPSI=col_ipsi,CONTRA=col_contra)) +
    scale_size_manual(values=c("FALSE"=2.5,"TRUE"=4.5),
                      labels=c("FALSE"="artifact","TRUE"="robust"),
                      name="dinuc result") +
    scale_alpha_identity() +
    scale_shape_manual(values=c("Guidance/adhesion"=16,"CPE/translation"=17,
                                 "Nuclear TF"=4,"Epigenetic"=15),
                       name="molecule type") +
    annotate("text",x=-1,y=-12,label="Below diagonal:\ndinuc more stringent\n(composition artifact)",
             size=2.8,colour="grey50",fontface="italic") +
    annotate("text",x=-11,y=-1,label="Above diagonal:\ndinuc less stringent\n(interaction even stronger)",
             size=2.8,colour="grey50",fontface="italic",hjust=0) +
    labs(x="Mononucleotide shuffle Z-score",
         y="Dinucleotide shuffle Z-score",
         title="Validation of significant interactions by dinucleotide shuffle",
         subtitle=paste0(sum(dinuc$survives),"/",nrow(dinuc),
                         " pairs survive more stringent dinucleotide test\n",
                         "Points below diagonal = artefacts removed; ",
                         "points near diagonal = robust interactions")) +
    theme_bw(base_size=11) +
    theme(panel.grid.minor=element_blank())

  ggsave("outputs/rna_interactions/figures/Fig4_mono_vs_dinuc.pdf",
         p4,width=8,height=7)
  ggsave("outputs/rna_interactions/figures/Fig4_mono_vs_dinuc.svg",
         p4,width=8,height=7)
  cat("Saved: Fig4_mono_vs_dinuc.pdf\n")
}


## ════════════════════════════════════════════════════════════════════
## FIG 5: Length artifact — nuclear/epigenetic
## ════════════════════════════════════════════════════════════════════
tf_data <- res_nuc %>%
  dplyr::mutate(tf_len=tf_lens[target],
                tf_class=dplyr::coalesce(tf_class_map[target],"unknown")) %>%
  dplyr::filter(!is.na(tf_len))

r_all <- cor(tf_data$tf_len,tf_data$zscore,use="complete.obs")

p5 <- ggplot(tf_data,aes(x=tf_len,y=zscore)) +
  geom_hline(yintercept=0,colour="grey80",linewidth=0.3) +
  geom_hline(yintercept=-2,colour="grey40",linewidth=0.5,linetype="dashed") +
  geom_smooth(method="lm",colour="grey30",fill="grey88",se=TRUE,linewidth=0.8) +
  geom_point(aes(colour=query,shape=tf_class),size=3,alpha=0.85) +
  ggrepel::geom_text_repel(
    data=dplyr::filter(tf_data,zscore< -5),
    aes(label=paste0(query,"x",target)),
    size=2.5,colour="grey20",box.padding=0.3) +
  scale_colour_manual(values=lnc_col,name="lncRNA") +
  scale_shape_manual(
    values=c(IPSI_TF=16,CONTRA_TF=17,neutral_TF=4,NUCLEAR_EPI=15,unknown=3),
    labels=c(IPSI_TF="Ipsi TF",CONTRA_TF="Contra TF",
             neutral_TF="Neutral TF",NUCLEAR_EPI="Epigenetic",unknown="?"),
    name="Category") +
  annotate("text",x=max(tf_data$tf_len,na.rm=TRUE)*0.65,y=3,
           label=paste0("r = ",round(r_all,2)),
           size=4.5,colour="grey25",fontface="bold") +
  labs(x="3'UTR length (nt)",y="Z-score",
       title="Nuclear regulators: Z-scores driven by 3'UTR length",
       subtitle=paste0("r = ",round(r_all,2),
                       "  -- These proteins absent from growth cones --> excluded from main analysis")) +
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank())
ggsave("outputs/rna_interactions/figures/Fig5_nuclear_artifact.pdf",
       p5,width=9,height=5)
ggsave("outputs/rna_interactions/figures/Fig5_nuclear_artifact.svg",
       p5,width=9,height=5)
cat("Saved: Fig5_nuclear_artifact.pdf\n")


## ════════════════════════════════════════════════════════════════════
## FIG 6: CIS + SHARED special cases
## ════════════════════════════════════════════════════════════════════
special <- dplyr::filter(res,query=="H19",
  mol_type %in% c("CIS_LOCUS","SHARED","AXONAL_IPSI","AXONAL_CONTRA")) %>%
  dplyr::mutate(
    category=dplyr::case_when(
      mol_type=="CIS_LOCUS" ~ "Cis-locus (imprinted)",
      mol_type=="SHARED"    ~ "Shared iRGC+cRGC",
      mol_type=="AXONAL_IPSI"   ~ "Ipsi target",
      mol_type=="AXONAL_CONTRA" ~ "Contra target"
    ),
    category=factor(category,levels=c("Contra target","Ipsi target",
                                       "Cis-locus (imprinted)","Shared iRGC+cRGC")),
    label_text=dplyr::case_when(
      target=="Igf2"  ~ "Igf2 (same imprinted locus)",
      target=="Cxcr4" ~ "Cxcr4 (AU-rich+SHORT)",
      target=="Gda"   ~ "Gda (exp. validated ipsi)",
      specific&zscore< -3 ~ target,
      TRUE ~ ""
    )
  )

p6 <- ggplot(special,aes(x=reorder(target,zscore),y=zscore,fill=category)) +
  geom_col(width=0.7,alpha=0.85) +
  geom_hline(yintercept=-2,linetype="dashed",colour="grey40",linewidth=0.5) +
  geom_hline(yintercept=0,colour="grey70",linewidth=0.3) +
  ggrepel::geom_text_repel(
    data=dplyr::filter(special,label_text!=""),
    aes(label=label_text),size=2.5,colour="grey15",
    box.padding=0.4,max.overlaps=20,direction="x",nudge_x=4) +
  coord_flip() +
  scale_fill_manual(values=c("Ipsi target"=col_ipsi,"Contra target"=col_contra,
                              "Cis-locus (imprinted)"=col_cis,
                              "Shared iRGC+cRGC"=col_shared)) +
  labs(x=NULL,y="Z-score",fill=NULL,
       title="H19 -- all target categories",
       subtitle="Orange: cis-locus (H19/Igf2 imprinted)  |  Purple: Cxcr4 in both cell types") +
  theme_bw(base_size=10) +
  theme(panel.grid.major.y=element_blank(),legend.position="bottom")
ggsave("outputs/rna_interactions/figures/Fig6_special_cases.pdf",p6,width=8,height=10)
ggsave("outputs/rna_interactions/figures/Fig6_special_cases.svg",p6,width=8,height=10)
cat("Saved: Fig6_special_cases.pdf\n")

cat("\nDone. Figures in outputs/rna_interactions/figures/\n")


## ════════════════════════════════════════════════════════════════════
## FIG 4b: Mono vs Dinuc scatter — ALL pairs including Gda
## Three categories visible:
##   - Robust (both significant)
##   - Artifact (mono sig, dinuc not)
##   - Borderline (dinuc sig but mono not — Gda)
## ════════════════════════════════════════════════════════════════════

if (!is.null(dinuc_all) && nrow(dinuc_all) > 0) {

  ## Merge all dinuc results with mono results and molecular categories
  all_results <- dinuc_all %>%
    dplyr::left_join(
      res %>%
        dplyr::distinct(query, target, mol_type, lnc_class, zscore) %>%
        dplyr::rename(mono_z_res=zscore),
      by=c("query","target")
    ) %>%
    dplyr::mutate(
      ## Use mono_z from dinuc CSV (more reliable than joining)
      mol_cat = dplyr::case_when(
        target %in% RBP_CPE     ~ "CPE/translation",
        target %in% NUCLEAR_TF  ~ "Nuclear TF",
        target %in% NUCLEAR_EPI ~ "Epigenetic",
        group   == "CONTRA"     ~ "Guidance (contra)",
        group   == "IPSI"       ~ "Guidance (ipsi)",
        TRUE                    ~ group
      ),
      status = dplyr::case_when(
        target == "Gda"         ~ "Borderline\n(mono n.s., dinuc p<0.05)",
        survives & mono_z < -2  ~ "Robust\n(both methods)",
        !survives & mono_z < -2 ~ "Artifact\n(mono only)",
        TRUE                    ~ "Not significant"
      ),
      status = factor(status, levels=c(
        "Robust\n(both methods)",
        "Borderline\n(mono n.s., dinuc p<0.05)",
        "Artifact\n(mono only)",
        "Not significant"
      )),
      label = dplyr::case_when(
        target == "Gda"             ~ "Gda*",
        survives & dinuc_z < -3.5   ~ paste0(query,"x",target),
        !survives & mono_z < -4     ~ paste0(query,"x",target),
        TRUE                        ~ ""
      )
    )

  status_col <- c(
    "Robust\n(both methods)"              = "#1D9E75",
    "Borderline\n(mono n.s., dinuc p<0.05)" = col_cis,
    "Artifact\n(mono only)"               = "grey60",
    "Not significant"                     = "grey85"
  )

  p4b <- ggplot(all_results, aes(x=mono_z, y=dinuc_z)) +
    ## Reference lines
    geom_abline(slope=1, intercept=0, linetype="dashed",
                colour="grey70", linewidth=0.5) +
    geom_hline(yintercept=-2, colour="grey50", linewidth=0.4,
               linetype="dotted") +
    geom_vline(xintercept=-2, colour="grey50", linewidth=0.4,
               linetype="dotted") +
    ## Points
    geom_point(aes(colour=status, shape=mol_cat, size=status),
               alpha=0.85) +
    ## Labels
    ggrepel::geom_text_repel(
      data=dplyr::filter(all_results, label!=""),
      aes(label=label, colour=status),
      size=2.8, box.padding=0.4, max.overlaps=20,
      fontface=ifelse(dplyr::filter(all_results, label!="")$target=="Gda",
                      "bold", "plain")
    ) +
    ## Gda special annotation
    annotate("text", x=0.34, y=-2.69,
             label="  Gda: co-expressed in ipsi axons\n  mono Z=+0.34, dinuc Z=-2.69\n  borderline — needs wet-lab validation",
             hjust=0, vjust=1, size=2.5, colour=col_cis,
             fontface="italic") +
    scale_colour_manual(values=status_col, name="Interaction status") +
    scale_shape_manual(
      values=c("Guidance (contra)"=16, "Guidance (ipsi)"=17,
               "CPE/translation"=15, "Epigenetic"=18,
               "Nuclear TF"=4, "SHARED"=8, "CIS_LOCUS"=3),
      name="Molecule type"
    ) +
    scale_size_manual(
      values=c("Robust\n(both methods)"=4.5,
               "Borderline\n(mono n.s., dinuc p<0.05)"=4.5,
               "Artifact\n(mono only)"=3,
               "Not significant"=2),
      guide="none"
    ) +
    ## Quadrant labels
    annotate("rect", xmin=-Inf, xmax=-2, ymin=-Inf, ymax=-2,
             fill="#1D9E75", alpha=0.04) +
    annotate("rect", xmin=-2, xmax=Inf, ymin=-2, ymax=Inf,
             fill="grey50", alpha=0.04) +
    annotate("text", x=-9, y=-0.5,
             label="Artifacts removed\nby dinuc shuffle",
             size=2.8, colour="grey55", fontface="italic") +
    annotate("text", x=-9, y=-12,
             label="Robust interactions\n(survive both methods)",
             size=2.8, colour="#1D9E75", fontface="italic") +
    labs(
      x="Mononucleotide shuffle Z-score (N=20)",
      y="Dinucleotide shuffle Z-score (N=50-100)",
      title="RNA-RNA interaction robustness: mononucleotide vs dinucleotide shuffle",
      subtitle=paste0(
        sum(all_results$status=="Robust\n(both methods)"),
        " robust  |  ",
        sum(all_results$status=="Artifact\n(mono only)"),
        " artifacts removed  |  ",
        "1 borderline (Gda)\n",
        "Dinucleotide shuffle preserves AU dinucleotide frequency — ",
        "more stringent test for composition-independent interactions"
      )
    ) +
    theme_bw(base_size=11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.key.size  = unit(0.6, "cm"),
      plot.subtitle    = element_text(size=8.5, colour="grey40")
    )

  ggsave("outputs/rna_interactions/figures/Fig4_mono_vs_dinuc.pdf",
         p4b, width=10, height=8)
  ggsave("outputs/rna_interactions/figures/Fig4_mono_vs_dinuc.svg",
         p4b, width=10, height=8)
  cat("Saved: Fig4_mono_vs_dinuc.pdf\n")
  cat("  Robust:", sum(all_results$status=="Robust\n(both methods)"), "\n")
  cat("  Borderline (Gda):", sum(all_results$target=="Gda"), "\n")
  cat("  Artifacts:", sum(all_results$status=="Artifact\n(mono only)"), "\n")
}


## ════════════════════════════════════════════════════════════════════

## ════════════════════════════════════════════════════════════════════
## FIG 8: Wang et al. 2016 vs Our TRAP-seq comparison
## ════════════════════════════════════════════════════════════════════

wang_data <- tribble(
  ~gene,     ~wang_2016,      ~our_axon,         ~dinuc_z, ~mol_class,    ~notes,
  "EphB1",   "IPSI",          "IPSI",             NA,       "Guidance",    "Canonical ipsi; EphB1-EphrinB2",
  "Nrp2",    "IPSI",          "IPSI",             NA,       "Guidance",    "Neuropilin-2",
  "Boc",     "IPSI",          "IPSI",             NA,       "Guidance",    "Shh co-receptor",
  "Nr2f2",   "IPSI",          "IPSI",             -3.26,    "TF",          "H19 interaction dinuc",
  "Igf2",    "IPSI",          "IPSI",             -7.59,    "Signaling",   "CIS locus with H19",
  "Plxna1",  "CONTRA",        "CONTRA",           -3.97,    "Guidance",    "Sema6D receptor cRGC",
  "Nrcam",   "CONTRA",        "CONTRA",           NA,       "Guidance",    "NrCAM Sema6D co-receptor",
  "Isl2",    "CONTRA",        "CONTRA",           NA,       "TF",          "~30pct cRGCs",
  "Gda",     "not detected",  "IPSI",             -2.69,    "Adhesion",    "Exp validated IHC borderline",
  "Opcml",   "not detected",  "IPSI",             NA,       "Adhesion",    "IgLON adhesion",
  "Tenm1",   "not detected",  "IPSI",             NA,       "Guidance",    "Teneurin",
  "Tenm2",   "not detected",  "IPSI",             -3.70,    "Guidance",    "Rian dinuc",
  "Tenm3",   "not detected",  "IPSI",             -2.40,    "Guidance",    "H19 dinuc",
  "Robo3",   "not detected",  "IPSI",             NA,       "Guidance",    "Slit receptor",
  "Ryk",     "not detected",  "IPSI",             NA,       "Guidance",    "Wnt receptor",
  "Grin2b",  "not detected",  "IPSI",             -3.54,    "Receptor",    "NMDA receptor H19 dinuc",
  "Lhx1",    "not detected",  "IPSI",             -4.20,    "TF",          "ipsi TF H19 dinuc",
  "Rgs7bp",  "not detected",  "CONTRA",           -5.96,    "Effector",    "Exp validated contra H19 Meg3 Mirg",
  "Plxna4",  "not detected",  "CONTRA",           -4.49,    "Guidance",    "H19 dinuc",
  "Islr2",   "not detected",  "CONTRA",           NA,       "Guidance",    "LRR receptor",
  "Sema5a",  "not detected",  "CONTRA",           NA,       "Guidance",    "Semaphorin",
  "Epha4",   "not detected",  "CONTRA",           NA,       "Guidance",    "Eph receptor",
  "Unc5c",   "not detected",  "CONTRA",           NA,       "Guidance",    "Netrin receptor",
  "Robo1",   "not detected",  "CONTRA",           NA,       "Guidance",    "Slit receptor",
  "Robo2",   "not detected",  "CONTRA",           NA,       "Guidance",    "Slit receptor",
  "Slit2",   "not detected",  "CONTRA",           NA,       "Guidance",    "Slit ligand",
  "L1cam",   "not detected",  "CONTRA",           NA,       "Adhesion",    "L1CAM",
  "Cntn2",   "not detected",  "CONTRA",           NA,       "Adhesion",    "Contactin-2",
  "Celsr3",  "not detected",  "CONTRA",           NA,       "Guidance",    "PCP pathway",
  "Zic2",    "IPSI",          "not detected",     NA,       "Nuclear TF",  "Master ipsi TF nuclear",
  "Lhx2",    "IPSI",          "not detected",     NA,       "Nuclear TF",  "ipsi TF weak nuclear",
  "Sox2",    "IPSI",          "not detected",     NA,       "Nuclear TF",  "Progenitor marker nuclear",
  "Atoh7",   "IPSI",          "not detected",     NA,       "Nuclear TF",  "Progenitor TF nuclear",
  "Ccnd2",   "IPSI",          "not detected",     NA,       "Nuclear TF",  "Cyclin D2 nuclear",
  "Igfbp5",  "IPSI",          "not detected",     NA,       "Signaling",   "IGF-binding protein",
  "Slc6a4",  "IPSI",          "not detected",     NA,       "Marker",      "SERT our Cre driver",
  "Igf1",    "CONTRA",        "not detected",     NA,       "Signaling",   "Contra in Wang 2016"
) %>%
  dplyr::mutate(
    category = dplyr::case_when(
      wang_2016 == "IPSI"   & our_axon == "IPSI"   ~ "Consistent IPSI",
      wang_2016 == "CONTRA" & our_axon == "CONTRA"  ~ "Consistent CONTRA",
      wang_2016 == "not detected" & our_axon == "IPSI"   ~ "Novel axonal - IPSI",
      wang_2016 == "not detected" & our_axon == "CONTRA"  ~ "Novel axonal - CONTRA",
      our_axon  == "not detected" ~ "Cell body only",
      TRUE ~ "other"
    ),
    has_interaction = !is.na(dinuc_z),
    wang_2016 = factor(wang_2016,
                       levels = c("IPSI","CONTRA","not detected")),
    our_axon  = factor(our_axon,
                       levels = c("IPSI","CONTRA","not detected"))
  )

cat_col8 <- c(
  "Consistent IPSI"     = col_ipsi,
  "Consistent CONTRA"   = col_contra,
  "Novel axonal - IPSI"  = "#1A7A50",
  "Novel axonal - CONTRA"= "#9B1D4A",
  "Cell body only"      = "grey65"
)

ax_labels    <- c(IPSI="Ipsi\n(our TRAP axon)",
                  CONTRA="Contra\n(our TRAP axon)",
                  "not detected"="Not in\nour data")
wang_labels  <- c(IPSI="Ipsi cell body\n(Wang 2016)",
                  CONTRA="Contra cell body\n(Wang 2016)",
                  "not detected"="Not in\nWang 2016")

n_consistent <- sum(grepl("Consistent", wang_data$category))
n_novel      <- sum(grepl("Novel axonal", wang_data$category))
n_soma       <- sum(wang_data$category == "Cell body only")
n_interact   <- sum(wang_data$has_interaction)

p8 <- ggplot(wang_data, aes(x=our_axon, y=wang_2016, colour=category)) +
  ## Diagonal highlight rectangles
  annotate("rect", xmin=0.5, xmax=1.5, ymin=0.5, ymax=1.5,
           fill=col_ipsi, alpha=0.07) +
  annotate("rect", xmin=1.5, xmax=2.5, ymin=1.5, ymax=2.5,
           fill=col_contra, alpha=0.07) +
  ## Diagonal line
  annotate("segment", x=0.5, xend=2.55, y=0.5, yend=2.55,
           colour="grey55", linewidth=0.7, linetype="dashed") +
  ## Diagonal labels
  annotate("text", x=1.0, y=0.72, label="Consistent\nIPSI",
           size=2.8, colour=col_ipsi, fontface="bold", alpha=0.8) +
  annotate("text", x=2.0, y=1.72, label="Consistent\nCONTRA",
           size=2.8, colour=col_contra, fontface="bold", alpha=0.8) +
  ## Points
  geom_jitter(aes(size=has_interaction, shape=has_interaction),
              width=0.2, height=0.2, alpha=0.9, seed=42) +
  ## Labels for key genes
  ggrepel::geom_text_repel(
    data = dplyr::filter(wang_data,
                         has_interaction |
                         category %in% c("Consistent IPSI",
                                         "Consistent CONTRA") |
                         gene %in% c("Gda","Rgs7bp","Unc5c","Plxna1",
                                     "Zic2","Grin2b","Tenm2")),
    aes(label = gene),
    size = 2.8, colour = "grey15",
    box.padding = 0.4, max.overlaps = 30, seed = 42
  ) +
  scale_colour_manual(values = cat_col8, name = "Category") +
  scale_size_manual(values  = c("FALSE"=3.0, "TRUE"=5.5),
                    labels  = c("FALSE"="expression only",
                                "TRUE"="+ RNA-RNA interaction (dinuc)"),
                    name    = "Evidence") +
  scale_shape_manual(values = c("FALSE"=16, "TRUE"=18), guide="none") +
  scale_x_discrete(labels = ax_labels) +
  scale_y_discrete(labels = wang_labels) +
  labs(
    x = "Our TRAP-seq (axons)",
    y = "Wang et al. 2016 (cell bodies)",
    title = "Axonal TRAP-seq vs cell body microarray (Wang et al. 2016 eNeuro)",
    subtitle = paste0(
      n_consistent, " genes consistent (on diagonal)  |  ",
      n_novel,      " novel axon-specific candidates  |  ",
      n_soma,       " cell body only (nuclear)\n",
      n_interact,   " with significant RNA-RNA interaction (large diamonds)  |  ",
      "* Gda: borderline dinuc Z=-2.69, needs validation"
    )
  ) +
  theme_bw(base_size=11) +
  theme(
    axis.text        = element_text(size=9),
    legend.position  = "right",
    plot.subtitle    = element_text(size=8.5, colour="grey40"),
    panel.grid.major = element_line(colour="grey92")
  )

ggsave("outputs/rna_interactions/figures/Fig8_Wang2016_comparison.pdf",
       p8, width=10, height=8)
ggsave("outputs/rna_interactions/figures/Fig8_Wang2016_comparison.svg",
       p8, width=10, height=8)

## Save comparison table
wang_data %>%
  dplyr::arrange(category, our_axon, wang_2016) %>%
  dplyr::mutate(
    dinuc_z_str = ifelse(is.na(dinuc_z), "n.t.",
                  ifelse(gene=="Gda",
                         paste0(sprintf("%.2f", dinuc_z),"*"),
                         sprintf("%.2f", dinuc_z)))
  ) %>%
  dplyr::select(Gene=gene, Category=category,
                Wang_2016=wang_2016, Our_axon=our_axon,
                Dinuc_Z=dinuc_z_str, Molecule=mol_class, Notes=notes) %>%
  write.csv("outputs/rna_interactions/tables/Table_Wang2016_comparison.csv",
            row.names=FALSE)

cat(sprintf("Saved: Fig8_Wang2016_comparison.pdf\n"))
cat(sprintf("  Consistent: %d | Novel axonal: %d | Soma only: %d | With interaction: %d\n",
            n_consistent, n_novel, n_soma, n_interact))
