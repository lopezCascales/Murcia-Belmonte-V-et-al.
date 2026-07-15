#!/usr/bin/env Rscript
## ─────────────────────────────────────────────────────────────────────────────
## Genome browser tracks — estilo Figuras 2 y 5 del paper
## Reproduce los perfiles de Ribo-RNAseq sobre regiones específicas
##
## Uso: Rscript plot_browser_tracks.R
## Requiere: Gviz, GenomicRanges, rtracklayer, BSgenome.Mmusculus.UCSC.mm39
## ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
  library(GenomicFeatures)
})

OUT_DIR <- "outputs/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ── Colores consistentes con el paper ────────────────────────────────────────
COL_SOMA_CONTRA <- "#C9375A"     ## rojo-rosado
COL_SOMA_IPSI   <- "#1F5B8E"     ## azul oscuro
COL_AXON_CONTRA <- "#E8829A"     ## rosa claro
COL_AXON_IPSI   <- "#3B8BD4"     ## azul claro (cyan en el paper)
COL_UTR_3       <- "#2D7D46"     ## verde (3'UTR boxes)
COL_HIGHLIGHT   <- "#F5F0D8"     ## amarillo muy claro (highlight de fondo)

## ── Función principal: plot de tracks para un gen ────────────────────────────
plot_gene_tracks <- function(
  gene_name,
  chrom, start, end,     ## coordenadas de la región a mostrar
  bam_files,             ## lista nombrada: list(soma_contra=path, axon_ipsi=path, ...)
  gtf_file,              ## path al GTF para anotación de genes
  highlight_start = NULL,## inicio del highlight (3'UTR o región de interés)
  highlight_end   = NULL,
  strand_focus    = "*", ## "+" o "-" para filtrar strand en tracks
  out_pdf         = NULL,
  width_kb        = NULL ## ancho en kb para la barra de escala
) {
  
  message(sprintf("Plotting: %s (%s:%d-%d)", gene_name, chrom, start, end))
  
  ## ── Track de anotación génica desde GTF ──────────────────────────────────
  txdb <- tryCatch(
    makeTxDbFromGFF(gtf_file, format = "GTF", organism = "Mus musculus"),
    error = function(e) NULL
  )
  
  if (!is.null(txdb)) {
    gene_track <- GeneRegionTrack(
      txdb,
      chromosome   = chrom,
      start        = start,
      end          = end,
      name         = "RefSeq\nGenes",
      geneSymbols  = TRUE,
      showId       = TRUE,
      transcriptAnnotation = "symbol",
      col          = "black",
      fill         = "black",
      fontcolor.group = "black",
      fontsize     = 7,
      cex.title    = 0.7,
      background.title = "transparent",
      col.title    = "black",
      rotation.title = 0
    )
  } else {
    ## Fallback: track vacío si no se puede cargar el GTF
    gene_track <- GenomeAxisTrack(name = "Genes")
  }
  
  ## ── Barra de escala ───────────────────────────────────────────────────────
  region_len_kb <- (end - start) / 1000
  scale_unit <- if (is.null(width_kb)) round(region_len_kb / 5) else width_kb
  
  ## ── Tracks de RNA-seq / Ribo-RNAseq ──────────────────────────────────────
  ## Definir colores y etiquetas
  track_specs <- list(
    soma_contra = list(col = COL_SOMA_CONTRA, label = "Soma\nContras"),
    soma_ipsi   = list(col = COL_SOMA_IPSI,   label = "Soma\nIpsis"),
    axon_contra = list(col = COL_AXON_CONTRA,  label = "Axon\nContras"),
    axon_ipsi   = list(col = COL_AXON_IPSI,    label = "Axon\nIpsis")
  )
  
  data_tracks <- list()
  for (type_name in names(bam_files)) {
    bam <- bam_files[[type_name]]
    if (!file.exists(bam)) {
      message("  AVISO: no encontrado ", bam)
      next
    }
    
    spec <- track_specs[[type_name]]
    trk  <- DataTrack(
      range    = bam,
      type     = "histogram",
      name     = spec$label,
      col.histogram = spec$col,
      fill.histogram = spec$col,
      alpha.title = 0.85,
      fontsize = 7,
      cex.title = 0.65,
      background.title = "transparent",
      col.title = spec$col,
      rotation.title = 0,
      chromosome = chrom,
      start = start,
      end   = end
    )
    data_tracks <- c(data_tracks, list(trk))
  }
  
  ## Si no hay BAMs disponibles, crear tracks de ejemplo con datos sintéticos
  if (length(data_tracks) == 0) {
    message("  AVISO: sin BAMs — generando tracks de ejemplo")
    ## En tu entorno real, los BAMs deben estar disponibles
    data_tracks <- lapply(names(track_specs), function(n) {
      spec <- track_specs[[n]]
      DataTrack(
        data = matrix(0, ncol = 10),
        start = seq(start, end, length.out = 10),
        end   = seq(start, end, length.out = 10) + 100,
        chromosome = chrom,
        name  = spec$label,
        col.histogram = spec$col,
        fill.histogram = spec$col,
        fontsize = 7, cex.title = 0.65,
        background.title = "transparent",
        col.title = spec$col,
        rotation.title = 0
      )
    })
  }
  
  ## ── Eje genómico ─────────────────────────────────────────────────────────
  axis_track <- GenomeAxisTrack(
    scale = scale_unit * 1000,
    labelPos = "above",
    fontsize = 7,
    col = "gray40"
  )
  
  ## ── Highlight de la región 3'UTR / interés ────────────────────────────────
  ht_track <- NULL
  if (!is.null(highlight_start) && !is.null(highlight_end)) {
    ht_track <- HighlightTrack(
      trackList = data_tracks,
      start     = highlight_start,
      end       = highlight_end,
      chromosome = chrom,
      col        = "transparent",
      fill       = COL_HIGHLIGHT,
      alpha      = 0.6,
      inBackground = TRUE
    )
  }
  
  ## ── Ensamblar y guardar ───────────────────────────────────────────────────
  all_tracks <- c(
    list(axis_track),
    if (is.null(ht_track)) data_tracks else list(ht_track),
    list(gene_track)
  )
  
  if (is.null(out_pdf)) {
    out_pdf <- file.path(OUT_DIR, sprintf("browser_%s.pdf", gene_name))
  }
  
  pdf(out_pdf, width = 6, height = length(all_tracks) * 0.9 + 0.5)
  tryCatch(
    plotTracks(
      all_tracks,
      from          = start,
      to            = end,
      chromosome    = chrom,
      sizes         = c(0.5, rep(1, length(all_tracks) - 2), 1.2),
      title.width   = 1.3,
      main          = sprintf("%s   [%d]", 
                              gene_name,
                              round((end - start) / 1000)),
      cex.main      = 0.9,
      fontface.main = "italic"
    ),
    error = function(e) message("Error en plotTracks: ", e$message)
  )
  dev.off()
  message("  Guardado: ", out_pdf)
}

## ── Ejemplo: definir regiones de las figuras del paper ───────────────────────
## Estas son las regiones que aparecen en Fig 2 y Fig 5

GTF_FILE <- "gencode.vM38.annotation.gtf"

## BAMs: ajusta paths a tu directorio
BAM_FILES <- list(
  soma_contra = "e15_5_CH1_both_soma.bam",
  soma_ipsi   = "htp55_STA_both_soma.bam",
  axon_contra = "e15_5_CH1_both.txt",   ## o el .bw si tienes bigWig
  axon_ipsi   = "htp55_STA_both.txt"
)

## Si tienes BigWig en lugar de BAM:
## BAM_FILES <- list(
##   soma_contra = "Pou4f2_TRAP_SOMA_merge_ht2.bw",
##   axon_contra = "Pou4f2_TRAP_AXON_merge_ht2.bw",
##   soma_ipsi   = "Slc6a4TRAPsoma_ht2.bw",
##   axon_ipsi   = "Slc6a4_TRAP_AXON_merge_ht2.bw"
## )

## ── Regiones de Fig 2 ─────────────────────────────────────────────────────────
regions_fig2 <- list(
  ## Lrp2: false positive de CPE
  Lrp2   = list(chrom="chr2",  start=69500000, end=69800000,
                 highlight_start=69780000, highlight_end=69800000),
  ## Nbeal1: CONTRA axon
  Nbeal1 = list(chrom="chr2",  start=65800000, end=66200000,
                 highlight_start=66180000, highlight_end=66200000),
  ## Igf2: interacción cis con H19
  Igf2   = list(chrom="chr7",  start=149800000, end=149820000,
                 highlight_start=149810000, highlight_end=149820000)
)

## ── Regiones de Fig 5 (perfiles de lncRNAs clave) ────────────────────────────
regions_fig5 <- list(
  ## H19 + Mir675
  H19    = list(chrom="chr7",  start=142129000, end=142135000),
  ## Lhx1 + Lhx1os
  Lhx1   = list(chrom="chr5",  start=109670000, end=109680000),
  ## Meg3 / Rian cluster
  Meg3   = list(chrom="chr12", start=109540000, end=109600000),
  ## Mirg
  Mirg   = list(chrom="chr12", start=109640000, end=109660000)
)

## ── Generar figuras ───────────────────────────────────────────────────────────
## Descomenta las líneas que necesites:

message("\nGenerando tracks de Fig 2 (CPE + candidatos axonales)...")
for (gname in names(regions_fig2)) {
  r <- regions_fig2[[gname]]
  plot_gene_tracks(
    gene_name       = gname,
    chrom           = r$chrom,
    start           = r$start,
    end             = r$end,
    bam_files       = BAM_FILES,
    gtf_file        = GTF_FILE,
    highlight_start = r$highlight_start,
    highlight_end   = r$highlight_end,
    out_pdf         = file.path(OUT_DIR, sprintf("tracks_fig2_%s.pdf", gname))
  )
}

message("\nGenerando tracks de lncRNAs (Fig 5)...")
for (gname in names(regions_fig5)) {
  r <- regions_fig5[[gname]]
  plot_gene_tracks(
    gene_name = gname,
    chrom     = r$chrom,
    start     = r$start,
    end       = r$end,
    bam_files = BAM_FILES,
    gtf_file  = GTF_FILE,
    out_pdf   = file.path(OUT_DIR, sprintf("tracks_fig5_%s.pdf", gname))
  )
}

message("\nTodos los tracks guardados en: ", OUT_DIR)
message("\nNOTA: Si tienes BigWig (.bw) en lugar de BAM, ajusta BAM_FILES")
message("       con los paths a los .bw disponibles en el directorio")
