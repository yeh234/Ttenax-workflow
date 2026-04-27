#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args) >= 1) {
  args[[1]]
} else {
  "figures/Tt_Figure1_workflow_v7"
}

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

library(grid)

light_fill <- "#FBFBFB"
panel_fill <- "#FFFFFF"
line_col <- "#111111"
note_fill <- "#FFF8EC"
final_fill <- "#F3F8FB"

fmt_int <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

final_read_count <- 3726422
final_base_count <- 11161804683
final_gb <- final_base_count / 1e9
predicted_transcripts <- 33845
any_rna_supported_genes <- 22302
all_three_rna_supported_genes <- 17671
proteomics_supported_genes <- 489

draw_round_box <- function(x, y, width, height, label, fontsize = 14,
                           fontface = "plain", fill = light_fill,
                           lineheight = 1.05, lwd = 2.2, lty = 1) {
  grid.roundrect(
    x = x, y = y, width = width, height = height,
    r = unit(0.018, "snpc"),
    gp = gpar(fill = fill, col = line_col, lwd = lwd, lty = lty)
  )
  grid.text(
    label = label, x = x, y = y,
    gp = gpar(
      fontsize = fontsize,
      fontface = fontface,
      col = line_col,
      lineheight = lineheight
    )
  )
}

draw_titled_box <- function(x, y, width, height, title, body,
                            title_size = 12, body_size = 9,
                            fill = light_fill) {
  grid.roundrect(
    x = x, y = y, width = width, height = height,
    r = unit(0.018, "snpc"),
    gp = gpar(fill = fill, col = line_col, lwd = 2.2)
  )

  left_x <- x - width / 2 + 0.018

  grid.text(
    label = title,
    x = left_x, y = y + height / 2 - 0.02,
    just = c("left", "top"),
    gp = gpar(
      fontsize = title_size,
      fontface = "bold",
      col = line_col,
      lineheight = 1.0
    )
  )

  grid.text(
    label = body,
    x = left_x, y = y + height / 2 - 0.06,
    just = c("left", "top"),
    gp = gpar(fontsize = body_size, col = line_col, lineheight = 1.06)
  )
}

draw_arrow <- function(x0, y0, x1, y1, lwd = 2.2) {
  grid.lines(
    x = unit(c(x0, x1), "npc"),
    y = unit(c(y0, y1), "npc"),
    gp = gpar(col = line_col, lwd = lwd),
    arrow = arrow(length = unit(0.09, "inches"), type = "closed")
  )
}

panel_a_boxes <- list(
  list(
    text = "Genomic DNA\nKAPA cleanup, >1 kb enrichment,\nRNase treatment, QC",
    size = 10.0
  ),
  list(
    text = "PromethION R10.4.1 sequencing\nTwo flow cells: PBE76679 and PBG38932",
    size = 10.0
  ),
  list(
    text = paste0(
      "Dorado SUP basecalling\n",
      "FASTQ + move-aware BAM; ",
      fmt_int(final_read_count), " reads, ",
      sprintf("%.2f", final_gb), " Gb"
    ),
    size = 9.4
  ),
  list(
    text = "Primary assembly strategy comparison\nFlye full-read, asm60x, length-filtered,\nand long-biased runs",
    size = 9.4
  ),
  list(
    text = "Read remapping and depth-based diagnosis\n102.56x single-copy depth proxy;\nestimated genome size 108.83 Mb",
    size = 8.8
  ),
  list(
    text = "Haplotig reduction\n34 haplotig-like contigs removed to build\nan annotation-friendly draft assembly",
    size = 8.8
  ),
  list(
    text = "Archived single-round polish\nCreation of Ttenax.nohap.dorado_polished.fasta\nas a released reference copy",
    size = 8.9
  ),
  list(
    text = "Two-round Dorado accuracy polish\nPrimary-only alignments, round 1 and round 2,\nwith residual-variant comparison",
    size = 8.6
  ),
  list(
    text = "BUSCO / MetaEuk completeness assessment\nBenchmarking of polish1 and polish2 references",
    size = 8.7
  ),
  list(
    text = "Comparative homology mapping\nminiprot alignment of the T. vaginalis proteome\nto the polished T. tenax reference",
    size = 8.4
  ),
  list(
    text = paste0(
      "Structural and functional annotation\n",
      "GALBA/gffread generated ", fmt_int(predicted_transcripts), " transcripts;\n",
      "EggNOG, Tv rescue, InterProScan, evidence merge"
    ),
    size = 7.9
  ),
  list(
    text = paste0(
      "Multi-omics validation\n",
      "DRS, SR/LR cDNA, and proteomics evidence;\n",
      fmt_int(any_rna_supported_genes), " RNA-supported genes, ",
      fmt_int(all_three_rna_supported_genes), " all-three RNA, ",
      fmt_int(proteomics_supported_genes), " proteomics"
    ),
    size = 7.1
  )
)

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(x = 0, y = 0, width = 1, height = 1,
                        just = c("left", "bottom")))

  grid.rect(gp = gpar(fill = panel_fill, col = NA))

  grid.text(
    "Figure 1. Study overview and data generation workflow",
    x = 0.5, y = 0.965,
    gp = gpar(fontsize = 22, fontface = "bold", col = line_col)
  )
  grid.text(
    "Workflow overview for long-read assembly, polishing, and annotation of T. tenax",
    x = 0.5, y = 0.935,
    gp = gpar(fontsize = 11.5, col = line_col)
  )

  grid.text("A", x = 0.045, y = 0.885,
            gp = gpar(fontsize = 20, fontface = "bold", col = line_col))
  grid.text("Sample-to-data and validation overview", x = 0.255, y = 0.885,
            gp = gpar(fontsize = 15, fontface = "bold", col = line_col))

  grid.text("B", x = 0.61, y = 0.885,
            gp = gpar(fontsize = 20, fontface = "bold", col = line_col))
  grid.text("Retained products and downstream annotation", x = 0.80, y = 0.885,
            gp = gpar(fontsize = 15, fontface = "bold", col = line_col))

  panel_a_x <- 0.27
  panel_a_w <- 0.43
  panel_a_h <- 0.045
  panel_a_gap <- 0.0085
  panel_a_y <- seq(0.81, by = -(panel_a_h + panel_a_gap), length.out = length(panel_a_boxes))

  for (i in seq_along(panel_a_boxes)) {
    draw_round_box(
      x = panel_a_x,
      y = panel_a_y[[i]],
      width = panel_a_w,
      height = panel_a_h,
      label = panel_a_boxes[[i]]$text,
      fontsize = panel_a_boxes[[i]]$size,
      fill = light_fill
    )

    if (i < length(panel_a_boxes)) {
      draw_arrow(
        x0 = panel_a_x,
        y0 = panel_a_y[[i]] - panel_a_h / 2,
        x1 = panel_a_x,
        y1 = panel_a_y[[i + 1]] + panel_a_h / 2
      )
    }
  }

  parent_x <- 0.80
  parent_y <- 0.76
  parent_w <- 0.31
  parent_h <- 0.09

  left_prod_x <- 0.69
  left_prod_y <- 0.57
  left_prod_w <- 0.17
  left_prod_h <- 0.19

  archived_x <- 0.89
  archived_y <- 0.64
  archived_w <- 0.16
  archived_h <- 0.145

  polish_x <- 0.89
  polish_y <- 0.46
  polish_w <- 0.16
  polish_h <- 0.12

  final_x <- 0.89
  final_y <- 0.285
  final_w <- 0.16
  final_h <- 0.18

  output_y <- 0.125
  output_w <- 0.31
  output_h <- 0.12

  draw_round_box(
    x = parent_x,
    y = parent_y,
    width = parent_w,
    height = parent_h,
    label = "Selected retained products\nand downstream annotation states",
    fontsize = 11.3,
    fontface = "bold",
    fill = light_fill
  )

  draw_arrow(
    x0 = parent_x - 0.01, y0 = parent_y - parent_h / 2,
    x1 = left_prod_x, y1 = left_prod_y + left_prod_h / 2
  )
  draw_arrow(
    x0 = parent_x + 0.01, y0 = parent_y - parent_h / 2,
    x1 = archived_x, y1 = archived_y + archived_h / 2
  )

  draw_titled_box(
    x = left_prod_x,
    y = left_prod_y,
    width = left_prod_w,
    height = left_prod_h,
    title = "Structural backbone",
    body = paste(
      "flye_len15k/assembly.fasta",
      "101.43 Mb, 212 contigs",
      "N50 924.5 kb",
      "",
      "Primary uses:",
      "- repeat-collapse assessment",
      "- exploratory genome structure",
      "- comparative backbone analyses",
      sep = "\n"
    ),
    title_size = 11.4,
    body_size = 8.3,
    fill = note_fill
  )

  draw_titled_box(
    x = archived_x,
    y = archived_y,
    width = archived_w,
    height = archived_h,
    title = "Archived single-round reference",
    body = paste(
      "Ttenax.nohap.",
      "dorado_polished.fasta.gz",
      "99.19 Mb, 178 contigs",
      "N50 939.5 kb",
      "",
      "released reference state",
      sep = "\n"
    ),
    title_size = 10.2,
    body_size = 7.2,
    fill = light_fill
  )

  draw_titled_box(
    x = polish_x,
    y = polish_y,
    width = polish_w,
    height = polish_h,
    title = "Two-round accuracy polish",
    body = paste(
      "Dorado aligner + primary-only BAM",
      "round 1 -> round 2",
      "98.40% mapped reads",
      "strict INDEL bp 31,295 -> 26,323",
      sep = "\n"
    ),
    title_size = 9.9,
    body_size = 6.8,
    fill = light_fill
  )

  draw_titled_box(
    x = final_x,
    y = final_y,
    width = final_w,
    height = final_h,
    title = "Final local annotation reference",
    body = paste(
      "dorado_polish_gpu1/",
      "draft.polish2.bs6.fasta",
      "99.19 Mb, 178 contigs",
      "",
      "used for:",
      "- miniprot homology mapping",
      "- GALBA + gffread",
      "- EggNOG / Tv rescue /",
      "  InterProScan merge",
      sep = "\n"
    ),
    title_size = 9.8,
    body_size = 6.7,
    fill = final_fill
  )

  draw_arrow(
    x0 = archived_x, y0 = archived_y - archived_h / 2,
    x1 = polish_x, y1 = polish_y + polish_h / 2
  )
  draw_arrow(
    x0 = polish_x, y0 = polish_y - polish_h / 2,
    x1 = final_x, y1 = final_y + final_h / 2
  )

  grid.roundrect(
    x = parent_x, y = output_y,
    width = output_w, height = output_h,
    r = unit(0.018, "snpc"),
    gp = gpar(fill = "#FFFFFF", col = line_col, lwd = 2.2)
  )

  grid.lines(
    x = unit(c(parent_x, parent_x), "npc"),
    y = unit(c(output_y - output_h / 2 + 0.01, output_y + output_h / 2 - 0.01), "npc"),
    gp = gpar(col = line_col, lwd = 1.8)
  )

  draw_arrow(
    x0 = left_prod_x, y0 = left_prod_y - left_prod_h / 2,
    x1 = parent_x - output_w / 4, y1 = output_y + output_h / 2
  )
  draw_arrow(
    x0 = final_x, y0 = final_y - final_h / 2,
    x1 = parent_x + output_w / 4, y1 = output_y + output_h / 2
  )

  grid.text(
    "Structural applications",
    x = parent_x - output_w / 4, y = output_y + 0.045,
    gp = gpar(fontsize = 11.5, fontface = "bold", col = line_col)
  )
  grid.text(
    paste(
      "- repeat-collapse assessment",
      "- exploratory genome structure",
      "- comparative analyses",
      sep = "\n"
    ),
    x = parent_x - output_w / 4, y = output_y - 0.005,
    gp = gpar(fontsize = 8.8, col = line_col, lineheight = 1.08)
  )

  grid.text(
    "Released annotation package",
    x = parent_x + output_w / 4, y = output_y + 0.045,
    gp = gpar(fontsize = 11.5, fontface = "bold", col = line_col)
  )
  grid.text(
    paste(
      "Tt_Final_Annotation.gff3",
      paste0(fmt_int(predicted_transcripts), " proteins and CDS"),
      paste0(fmt_int(predicted_transcripts), "-row annotation table"),
      "Tv rescue and InterPro support",
      sep = "\n"
    ),
    x = parent_x + output_w / 4, y = output_y - 0.003,
    gp = gpar(fontsize = 8.3, col = line_col, lineheight = 1.08)
  )

  grid.text(
    "This schematic highlights retention of an unpolished structural backbone assembly, archival of a",
    x = 0.5, y = 0.035,
    gp = gpar(fontsize = 10.2, col = line_col)
  )
  grid.text(
    "single-round reference, and derivation of a two-round local reference used for multi-layer annotation.",
    x = 0.5, y = 0.014,
    gp = gpar(fontsize = 10.2, col = line_col)
  )

  popViewport()
}

png_file <- paste0(output_prefix, ".png")
pdf_file <- paste0(output_prefix, ".pdf")

png(filename = png_file, width = 16, height = 10.5, units = "in", res = 320, bg = "white")
draw_figure()
dev.off()

pdf(file = pdf_file, width = 16, height = 10.5, bg = "white", useDingbats = FALSE)
draw_figure()
dev.off()

message("Wrote: ", png_file)
message("Wrote: ", pdf_file)
