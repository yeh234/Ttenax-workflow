#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args) >= 1) {
  args[[1]]
} else {
  "figures/Tt_Figure2_assembly_strategy_v2"
}

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

library(grid)

line_col <- "#111111"
muted_col <- "#585858"
axis_col <- "#303030"
bar_col <- "#2E6F9E"
selected_col <- "#D18422"
box_fill <- "#FBFBFB"
note_fill <- "#FFF8EC"

read_fasta_lengths <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)

  lengths <- numeric()
  current <- 0
  repeat {
    lines <- readLines(con, n = 10000, warn = FALSE)
    if (length(lines) == 0) break
    for (line in lines) {
      if (startsWith(line, ">")) {
        if (current > 0) lengths <- c(lengths, current)
        current <- 0
      } else {
        current <- current + nchar(line, type = "chars")
      }
    }
  }
  if (current > 0) lengths <- c(lengths, current)
  lengths
}

read_lengths <- function(path) {
  fai_path <- if (grepl("\\.fai$", path)) path else paste0(path, ".fai")
  if (file.exists(fai_path)) {
    tab <- read.delim(fai_path, header = FALSE, stringsAsFactors = FALSE)
    return(as.numeric(tab[[2]]))
  }
  if (!file.exists(path)) {
    stop("Missing required input: ", path)
  }
  read_fasta_lengths(path)
}

n50 <- function(lengths) {
  sorted <- sort(lengths, decreasing = TRUE)
  sorted[which(cumsum(sorted) >= sum(sorted) / 2)[1]]
}

assembly_stats <- function(label, path, selected = FALSE) {
  lengths <- read_lengths(path)
  data.frame(
    label = label,
    path = path,
    span_mb = sum(lengths) / 1e6,
    n50_kb = n50(lengths) / 1e3,
    contigs = length(lengths),
    selected = selected,
    stringsAsFactors = FALSE
  )
}

fmt1 <- function(x) sprintf("%.1f", x)
fmt2 <- function(x) sprintf("%.2f", x)
fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_pct <- function(x) sprintf("%.1f%%", x)

count_lines <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  length(readLines(path, warn = FALSE))
}

draw_panel_title <- function(letter, title, x, y) {
  grid.text(letter, x = x, y = y, just = c("left", "center"),
            gp = gpar(fontsize = 20, fontface = "bold", col = line_col))
  grid.text(title, x = x + 0.04, y = y, just = c("left", "center"),
            gp = gpar(fontsize = 16, fontface = "bold", col = line_col))
}

draw_axis_ticks <- function(x0, y0, width, ticks, max_value, label_size = 7.5) {
  for (tick in ticks) {
    tx <- x0 + width * tick / max_value
    grid.lines(unit(c(tx, tx), "npc"), unit(c(y0, y0 - 0.006), "npc"),
               gp = gpar(col = axis_col, lwd = 0.8))
    grid.text(fmt_int(tick), x = tx, y = y0 - 0.021,
              gp = gpar(fontsize = label_size, col = line_col))
  }
}

draw_horizontal_bars <- function(x0, y0, width, height, title, values, labels,
                                 selected, axis_title, value_formatter,
                                 tick_step, label_left = TRUE) {
  title_y <- y0 + height + 0.027
  grid.text(title, x = x0 + width / 2, y = title_y,
            gp = gpar(fontsize = 12.5, col = line_col))

  plot_x0 <- if (label_left) x0 + 0.092 else x0 + 0.010
  plot_x1 <- x0 + width - 0.012
  plot_w <- plot_x1 - plot_x0
  axis_y <- y0 + 0.026
  top_pad <- 0.018
  bottom_pad <- 0.065
  bar_area_h <- height - top_pad - bottom_pad
  n <- length(values)
  gap <- bar_area_h / n
  bar_h <- min(0.026, gap * 0.66)

  max_value <- ceiling(max(values) * 1.10 / tick_step) * tick_step
  ticks <- seq(0, max_value, by = tick_step)

  grid.rect(x = x0, y = y0, width = width, height = height,
            just = c("left", "bottom"),
            gp = gpar(fill = "#FFFFFF", col = line_col, lwd = 1.0))
  grid.lines(unit(c(plot_x0, plot_x1), "npc"), unit(c(axis_y, axis_y), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))
  grid.lines(unit(c(plot_x0, plot_x0), "npc"), unit(c(axis_y, y0 + height - top_pad), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))
  draw_axis_ticks(plot_x0, axis_y, plot_w, ticks, max_value)

  grid.text(axis_title, x = plot_x0 + plot_w / 2, y = y0 - 0.018,
            gp = gpar(fontsize = 8.4, col = line_col))

  for (i in seq_along(values)) {
    y <- y0 + height - top_pad - (i - 0.5) * gap
    bar_w <- plot_w * values[[i]] / max_value
    fill <- if (selected[[i]]) selected_col else bar_col

    if (label_left) {
      grid.text(labels[[i]], x = plot_x0 - 0.010, y = y,
                just = c("right", "center"),
                gp = gpar(fontsize = 8.9, col = line_col))
    }

    grid.rect(x = plot_x0, y = y, width = bar_w, height = bar_h,
              just = c("left", "center"),
              gp = gpar(fill = fill, col = NA))

    label_x <- min(plot_x0 + bar_w + 0.006, plot_x1 + 0.002)
    just <- if (label_x > plot_x1 - 0.002) c("right", "center") else c("left", "center")
    grid.text(value_formatter(values[[i]]), x = label_x, y = y,
              just = just,
              gp = gpar(fontsize = 8.2, col = line_col))
  }
}

draw_note_box <- function(x, y, width, height, body, title = "Comparison result:",
                          title_size = 8.6, body_size = 8.0,
                          lineheight = 1.05) {
  grid.roundrect(x = x, y = y, width = width, height = height,
                 r = unit(0.010, "snpc"),
                 gp = gpar(fill = note_fill, col = "#D9B36A", lwd = 0.8))
  grid.text(
    title,
    x = x, y = y + height * 0.30,
    gp = gpar(fontsize = title_size, fontface = "bold", col = line_col)
  )
  grid.text(
    body,
    x = x, y = y - height * 0.10,
    gp = gpar(fontsize = body_size, col = line_col, lineheight = lineheight)
  )
}

draw_round_box <- function(x, y, width, height, title, body,
                           title_size = 11.0, body_size = 8.6,
                           fill = box_fill, lwd = 1.4) {
  grid.roundrect(x = x, y = y, width = width, height = height,
                 r = unit(0.014, "snpc"),
                 gp = gpar(fill = fill, col = line_col, lwd = lwd))
  grid.text(title, x = x, y = y + height / 2 - 0.020,
            gp = gpar(fontsize = title_size, fontface = "bold",
                      col = line_col, lineheight = 0.95))
  grid.text(body, x = x, y = y + height / 2 - 0.055,
            gp = gpar(fontsize = body_size, col = line_col,
                      lineheight = 1.02))
}

draw_arrow <- function(x0, y0, x1, y1, lwd = 1.5) {
  grid.lines(unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
             gp = gpar(col = line_col, lwd = lwd),
             arrow = arrow(length = unit(0.11, "inches"), type = "closed"))
}

draw_vertical_bars <- function(x0, y0, width, height, title, values, labels,
                               y_label, max_value, tick_step, value_formatter) {
  grid.rect(x = x0, y = y0, width = width, height = height,
            just = c("left", "bottom"),
            gp = gpar(fill = "#FFFFFF", col = line_col, lwd = 1.0))
  grid.text(title, x = x0 + width / 2, y = y0 + height + 0.030,
            gp = gpar(fontsize = 12.0, col = line_col))

  left <- x0 + 0.055
  right <- x0 + width - 0.022
  bottom <- y0 + 0.060
  top <- y0 + height - 0.030
  plot_w <- right - left
  plot_h <- top - bottom

  grid.lines(unit(c(left, right), "npc"), unit(c(bottom, bottom), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))
  grid.lines(unit(c(left, left), "npc"), unit(c(bottom, top), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))

  ticks <- seq(0, max_value, by = tick_step)
  for (tick in ticks) {
    ty <- bottom + plot_h * tick / max_value
    grid.lines(unit(c(left - 0.006, left), "npc"), unit(c(ty, ty), "npc"),
               gp = gpar(col = axis_col, lwd = 0.8))
    grid.text(fmt_int(tick), x = left - 0.010, y = ty,
              just = c("right", "center"),
              gp = gpar(fontsize = 7.7, col = line_col))
  }

  grid.text(y_label, x = x0 + 0.018, y = bottom + plot_h / 2, rot = 90,
            gp = gpar(fontsize = 8.8, col = line_col))

  xs <- left + plot_w * c(0.28, 0.72)
  bar_w <- plot_w * 0.36
  fills <- c(bar_col, selected_col)

  for (i in seq_along(values)) {
    bar_h <- plot_h * values[[i]] / max_value
    grid.rect(x = xs[[i]], y = bottom, width = bar_w, height = bar_h,
              just = c("center", "bottom"),
              gp = gpar(fill = fills[[i]], col = NA))
    grid.text(value_formatter(values[[i]]), x = xs[[i]],
              y = min(bottom + bar_h + 0.014, top + 0.004),
              gp = gpar(fontsize = 8.0, col = line_col))
    grid.text(labels[[i]], x = xs[[i]], y = y0 + 0.022, rot = 24,
              gp = gpar(fontsize = 8.2, col = line_col))
  }
}

assemblies <- rbind(
  assembly_stats("Full HQ", "flye_hq/40-polishing/filtered_contigs.fasta.fai"),
  assembly_stats("Full raw", "flye_raw/40-polishing/filtered_contigs.fasta.fai"),
  assembly_stats("asm60x HQ", "flye_asm60x_hq/40-polishing/filtered_contigs.fasta.fai"),
  assembly_stats("asm60x raw", "flye_asm60x_raw/40-polishing/filtered_contigs.fasta.fai"),
  assembly_stats(">=15 kb", "flye_len15k/assembly.fasta", selected = TRUE),
  assembly_stats(">=20 kb", "flye_len20k/40-polishing/filtered_contigs.fasta.fai"),
  assembly_stats("Long-biased", "flye_longbiased_8p22Gb_hq/40-polishing/filtered_contigs.fasta.fai")
)

backbone <- assembly_stats("Structural backbone", "flye_len15k/assembly.fasta")
nohap <- assembly_stats("Haplotig-reduced draft", "flye_len15k/assembly.no_haplotig.fasta")
polish1 <- assembly_stats("Polish round 1", "dorado_polish_gpu1/draft.polish1.bs6.fasta.fai")
polish2 <- assembly_stats("Final annotation reference", "dorado_polish_gpu1/draft.polish2.bs6.fasta.fai")
archived <- assembly_stats("Archived single-round reference", "final_ref/Ttenax.nohap.dorado_polished.fasta.fai")

haplotig_removed <- if (file.exists("haplotig.ids")) length(readLines("haplotig.ids")) else NA_integer_
very_low_contigs <- count_lines("verylow.tsv")
collapsed_repeat_contigs <- count_lines("collapsed.repeat.tsv")
single_copy_depth <- 102.56
estimated_genome_mb <- 11161804683 / single_copy_depth / 1e6

polish_qc <- read.delim("polish_round_qc_summary.tsv", stringsAsFactors = FALSE)
qc_value <- function(metric, column) {
  row <- polish_qc[polish_qc$metric == metric, column, drop = TRUE]
  suppressWarnings(as.numeric(row[[1]]))
}

strict_counts <- c(
  qc_value("strict_indels_q30_dp20", "polish1"),
  qc_value("strict_indels_q30_dp20", "polish2")
)
strict_bp <- c(
  qc_value("strict_indel_bp_q30_dp20", "polish1"),
  qc_value("strict_indel_bp_q30_dp20", "polish2")
)
mapped_pct <- qc_value("primary_reads_mapped_pct", "polish2")

metrics_out <- rbind(
  assemblies,
  backbone,
  nohap,
  archived,
  polish1,
  polish2
)
write.table(metrics_out, paste0(output_prefix, "_source_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

depth_out <- data.frame(
  metric = c(
    "single_copy_depth_x",
    "estimated_genome_size_mb",
    "haplotig_like_contigs",
    "very_low_depth_contigs",
    "strong_high_depth_contigs"
  ),
  value = c(
    single_copy_depth,
    estimated_genome_mb,
    haplotig_removed,
    very_low_contigs,
    collapsed_repeat_contigs
  )
)
write.table(depth_out, paste0(output_prefix, "_depth_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  grid.text("Figure 2. Assembly strategy comparison and selection of final reference assemblies",
            x = 0.5, y = 0.966,
            gp = gpar(fontsize = 22, fontface = "bold", col = line_col))
  grid.text("Data-derived summary of assembly span, contiguity, fragmentation, haplotig reduction, and polishing improvement",
            x = 0.5, y = 0.935,
            gp = gpar(fontsize = 10.5, col = muted_col))

  draw_panel_title("A", "Assembly strategy comparison", 0.035, 0.880)

  draw_horizontal_bars(
    x0 = 0.055, y0 = 0.555, width = 0.285, height = 0.255,
    title = "Span", values = assemblies$span_mb, labels = assemblies$label,
    selected = assemblies$selected, axis_title = "Assembly size (Mb)",
    value_formatter = fmt1, tick_step = 25, label_left = TRUE
  )
  draw_horizontal_bars(
    x0 = 0.382, y0 = 0.555, width = 0.270, height = 0.255,
    title = "Contiguity", values = assemblies$n50_kb, labels = assemblies$label,
    selected = assemblies$selected, axis_title = "N50 (kb)",
    value_formatter = fmt1, tick_step = 250, label_left = FALSE
  )
  draw_horizontal_bars(
    x0 = 0.690, y0 = 0.555, width = 0.255, height = 0.255,
    title = "Fragmentation", values = assemblies$contigs, labels = assemblies$label,
    selected = assemblies$selected, axis_title = "Contigs",
    value_formatter = fmt_int, tick_step = 750, label_left = FALSE
  )

  draw_note_box(
    x = 0.500, y = 0.493, width = 0.720, height = 0.056,
    paste0(
      "The >=15 kb assembly was retained as the structural backbone.\n",
      fmt2(backbone$span_mb), " Mb, ", fmt_int(backbone$contigs),
      " contigs, N50 ", fmt1(backbone$n50_kb), " kb."
    ),
    title_size = 8.6,
    body_size = 8.2,
    lineheight = 1.05
  )

  draw_panel_title("B", "Backbone and annotation-oriented assembly products", 0.035, 0.450)

  b_y <- 0.305
  b_h <- 0.120
  b_w <- 0.177
  x1 <- 0.145
  x2 <- 0.345
  x3 <- 0.545

  draw_round_box(
    x1, b_y, b_w, b_h,
    "Structural backbone",
    paste0(
      "flye_len15k/assembly.fasta\n",
      fmt2(backbone$span_mb), " Mb; ", fmt_int(backbone$contigs), " contigs\n",
      "N50 ", fmt1(backbone$n50_kb), " kb\n",
      "Use: structural comparisons"
    ),
    fill = note_fill
  )
  draw_round_box(
    x2, b_y, b_w, b_h,
    "Annotation draft",
    paste0(
      "assembly.no_haplotig.fasta\n",
      fmt2(nohap$span_mb), " Mb; ", fmt_int(nohap$contigs), " contigs\n",
      "removed ", fmt_int(haplotig_removed), " haplotig-like contigs\n",
      "intermediate draft"
    )
  )
  draw_round_box(
    x3, b_y, b_w, b_h,
    "Final reference",
    paste0(
      "draft.polish2.bs6.fasta\n",
      fmt2(polish2$span_mb), " Mb; ", fmt_int(polish2$contigs), " contigs\n",
      "N50 ", fmt1(polish2$n50_kb), " kb\n",
      "Use: annotation release"
    ),
    fill = "#F3F8FB"
  )

  draw_arrow(x1 + b_w / 2 + 0.005, b_y, x2 - b_w / 2 - 0.005, b_y)
  draw_arrow(x2 + b_w / 2 + 0.005, b_y, x3 - b_w / 2 - 0.005, b_y)

  draw_note_box(
    x = 0.345, y = 0.105, width = 0.505, height = 0.074,
    paste0(
      "Depth screen: ", fmt1(single_copy_depth), "x single-copy depth; ",
      fmt2(estimated_genome_mb), " Mb genome-size estimate.\n",
      fmt_int(haplotig_removed), " haplotig-like and ",
      fmt_int(very_low_contigs), " very-low-depth contigs flagged.\n",
      "Polish2 was retained as the final local working reference."
    ),
    title_size = 8.2,
    body_size = 7.1,
    lineheight = 0.98
  )

  draw_panel_title("C", "ORF-oriented polishing comparison", 0.690, 0.450)

  draw_vertical_bars(
    x0 = 0.675, y0 = 0.185, width = 0.145, height = 0.175,
    title = "Strict INDEL count",
    values = strict_counts,
    labels = c("polish1", "polish2"),
    y_label = "count",
    max_value = 2500,
    tick_step = 500,
    value_formatter = fmt_int
  )

  draw_vertical_bars(
    x0 = 0.845, y0 = 0.185, width = 0.140, height = 0.175,
    title = "Strict INDEL bp",
    values = strict_bp,
    labels = c("polish1", "polish2"),
    y_label = "bp",
    max_value = 35000,
    tick_step = 10000,
    value_formatter = fmt_int
  )

  indel_drop <- strict_counts[[1]] - strict_counts[[2]]
  bp_drop <- strict_bp[[1]] - strict_bp[[2]]
  bp_drop_pct <- bp_drop / strict_bp[[1]] * 100

  draw_note_box(
    x = 0.830, y = 0.095, width = 0.285, height = 0.080,
    paste0(
      "Polish2 retained ", fmt_pct(mapped_pct),
      " mapped reads while\nreducing strict INDEL burden.\n",
      "Strict INDELs -", fmt_int(indel_drop),
      "; strict INDEL bp -", fmt_int(bp_drop), " (", fmt1(bp_drop_pct), "%)."
    ),
    title_size = 8.2,
    body_size = 7.1,
    lineheight = 0.98
  )

  grid.text(
    "Source data: local Flye FASTA/index files, depth-screen TSVs, haplotig IDs, polished FASTA indexes, and polish_round_qc_summary.tsv.",
    x = 0.5, y = 0.033,
    gp = gpar(fontsize = 8.2, col = muted_col)
  )
}

png(paste0(output_prefix, ".png"), width = 16, height = 10.5, units = "in", res = 300)
draw_figure()
dev.off()

pdf(paste0(output_prefix, ".pdf"), width = 16, height = 10.5, useDingbats = FALSE)
draw_figure()
dev.off()

message("Wrote: ", paste0(output_prefix, ".png"))
message("Wrote: ", paste0(output_prefix, ".pdf"))
message("Wrote: ", paste0(output_prefix, "_source_metrics.tsv"))
message("Wrote: ", paste0(output_prefix, "_depth_metrics.tsv"))
