#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args) >= 1) {
  args[[1]]
} else {
  "figures/Tt_Figure3_depth_validation_v1"
}

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

library(grid)

line_col <- "#111111"
axis_col <- "#303030"
muted_col <- "#5C5C5C"
grid_col <- "#E6E6E6"
note_fill <- "#FFF8EC"
panel_fill <- "#FFFFFF"
primary_col <- "#3A7CA5"
haplotig_col <- "#D18422"
verylow_col <- "#B44B3B"
repeat_col <- "#9E6EB0"
collapsed_col <- "#54278F"
final_col <- "#2E8B7C"

required_files <- c(
  "contig_mean_depth.tsv",
  "haplotig_like.tsv",
  "verylow.tsv",
  "repeat_gt1p5C.tsv",
  "collapsed.repeat.tsv",
  "flye_len15k/assembly.fasta",
  "flye_len15k/assembly.no_haplotig.fasta",
  "dorado_polish_gpu1/draft.polish2.bs6.fasta.fai"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files: ", paste(missing_files, collapse = ", "))
}

fmt1 <- function(x) sprintf("%.1f", x)
fmt2 <- function(x) sprintf("%.2f", x)
fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)

read_depth_table <- function(path) {
  tab <- read.delim(
    path,
    header = FALSE,
    col.names = c("contig", "depth_x", "length_bp"),
    stringsAsFactors = FALSE
  )
  tab$depth_x <- as.numeric(tab$depth_x)
  tab$length_bp <- as.numeric(tab$length_bp)
  tab
}

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
  read_fasta_lengths(path)
}

n50 <- function(lengths) {
  sorted <- sort(lengths, decreasing = TRUE)
  sorted[which(cumsum(sorted) >= sum(sorted) / 2)[1]]
}

assembly_stats <- function(label, path) {
  lengths <- read_lengths(path)
  data.frame(
    label = label,
    path = path,
    span_bp = sum(lengths),
    span_mb = sum(lengths) / 1e6,
    n50_bp = n50(lengths),
    n50_kb = n50(lengths) / 1e3,
    contigs = length(lengths),
    stringsAsFactors = FALSE
  )
}

count_lines <- function(path) length(readLines(path, warn = FALSE))

draw_panel_title <- function(letter, title, x, y) {
  grid.text(
    letter,
    x = x,
    y = y,
    just = c("left", "center"),
    gp = gpar(fontsize = 20, fontface = "bold", col = line_col)
  )
  grid.text(
    title,
    x = x + 0.035,
    y = y,
    just = c("left", "center"),
    gp = gpar(fontsize = 15.0, fontface = "bold", col = line_col)
  )
}

draw_panel_frame <- function(x0, y0, width, height) {
  grid.rect(
    x = x0,
    y = y0,
    width = width,
    height = height,
    just = c("left", "bottom"),
    gp = gpar(fill = panel_fill, col = NA, lwd = 1.0)
  )
}

draw_plot_axes <- function(left, bottom, right, top) {
  grid.lines(
    unit(c(left, right), "npc"),
    unit(c(bottom, bottom), "npc"),
    gp = gpar(col = axis_col, lwd = 0.9)
  )
  grid.lines(
    unit(c(left, left), "npc"),
    unit(c(bottom, top), "npc"),
    gp = gpar(col = axis_col, lwd = 0.9)
  )
}

draw_histogram_panel <- function(x0, y0, width, height, depths,
                                 haplotig_range, single_copy_depth,
                                 repeat_threshold, collapsed_threshold,
                                 counts) {
  draw_panel_frame(x0, y0, width, height)

  left <- x0 + 0.060
  right <- x0 + width - 0.025
  bottom <- y0 + 0.060
  top <- y0 + height - 0.040
  plot_w <- right - left
  plot_h <- top - bottom

  x_min <- 0
  x_max <- 600
  breaks <- seq(x_min, x_max, by = 25)
  hist_obj <- hist(depths, breaks = breaks, plot = FALSE, right = FALSE)
  y_max <- ceiling(max(hist_obj$counts) * 1.18 / 20) * 20

  map_x <- function(v) left + plot_w * (v - x_min) / (x_max - x_min)
  map_y <- function(v) bottom + plot_h * v / y_max

  for (tick in seq(0, y_max, by = 10)) {
    ty <- map_y(tick)
    grid.lines(
      unit(c(left, right), "npc"),
      unit(c(ty, ty), "npc"),
      gp = gpar(col = grid_col, lwd = 0.6)
    )
    if (tick %% 20 == 0) {
      grid.lines(
        unit(c(left - 0.006, left), "npc"),
        unit(c(ty, ty), "npc"),
        gp = gpar(col = axis_col, lwd = 0.8)
      )
      grid.text(
        fmt_int(tick),
        x = left - 0.010,
        y = ty,
        just = c("right", "center"),
        gp = gpar(fontsize = 7.5, col = line_col)
      )
    }
  }

  for (tick in seq(0, x_max, by = 100)) {
    tx <- map_x(tick)
    grid.lines(
      unit(c(tx, tx), "npc"),
      unit(c(bottom, bottom - 0.006), "npc"),
      gp = gpar(col = axis_col, lwd = 0.8)
    )
    grid.text(
      fmt_int(tick),
      x = tx,
      y = bottom - 0.020,
      gp = gpar(fontsize = 7.5, col = line_col)
    )
  }

  grid.rect(
    x = map_x(haplotig_range[1]),
    y = bottom,
    width = map_x(haplotig_range[2]) - map_x(haplotig_range[1]),
    height = plot_h,
    just = c("left", "bottom"),
    gp = gpar(fill = adjustcolor(haplotig_col, alpha.f = 0.14), col = NA)
  )
  grid.rect(
    x = map_x(repeat_threshold),
    y = bottom,
    width = map_x(x_max) - map_x(repeat_threshold),
    height = plot_h,
    just = c("left", "bottom"),
    gp = gpar(fill = adjustcolor(repeat_col, alpha.f = 0.12), col = NA)
  )

  draw_plot_axes(left, bottom, right, top)

  for (i in seq_along(hist_obj$counts)) {
    x_left <- map_x(hist_obj$breaks[[i]])
    x_right <- map_x(hist_obj$breaks[[i + 1]])
    bar_h <- map_y(hist_obj$counts[[i]]) - bottom
    grid.rect(
      x = x_left + 0.001,
      y = bottom,
      width = max(x_right - x_left - 0.002, 0.001),
      height = bar_h,
      just = c("left", "bottom"),
      gp = gpar(fill = primary_col, col = NA)
    )
  }

  threshold_lines <- data.frame(
    x = c(single_copy_depth, repeat_threshold, collapsed_threshold),
    col = c(line_col, repeat_col, collapsed_col),
    lty = c("dashed", "dashed", "dotted"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(threshold_lines))) {
    tx <- map_x(threshold_lines$x[[i]])
    grid.lines(
      unit(c(tx, tx), "npc"),
      unit(c(bottom, top), "npc"),
      gp = gpar(col = threshold_lines$col[[i]], lwd = 1.1, lty = threshold_lines$lty[[i]])
    )
  }

  grid.text(
    "Contig mean depth (x)",
    x = left + plot_w / 2,
    y = y0 + 0.018,
    gp = gpar(fontsize = 8.7, col = line_col)
  )
  grid.text(
    "Number of contigs",
    x = x0 + 0.018,
    y = bottom + plot_h / 2,
    rot = 90,
    gp = gpar(fontsize = 8.7, col = line_col)
  )

  grid.text(
    paste0("very-low\nn=", counts$verylow),
    x = map_x(22),
    y = map_y(y_max * 0.82),
    just = c("left", "center"),
    gp = gpar(fontsize = 7.4, col = verylow_col, lineheight = 0.95)
  )
  grid.text(
    paste0("haplotig-like n=", counts$haplotig, "\n",
           fmt1(haplotig_range[1]), "-", fmt1(haplotig_range[2]), "x"),
    x = map_x(mean(haplotig_range)),
    y = map_y(y_max * 0.94),
    gp = gpar(fontsize = 7.2, col = haplotig_col, lineheight = 0.95)
  )
  grid.text(
    paste0("single-copy proxy\n", fmt2(single_copy_depth), "x"),
    x = map_x(single_copy_depth) + 0.007,
    y = map_y(y_max * 0.68),
    just = c("left", "center"),
    gp = gpar(fontsize = 7.4, col = line_col, lineheight = 0.95)
  )
  grid.text(
    paste0(">1.5C repeat-like n=", counts$repeat15, "\n",
           ">2C collapsed n=", counts$collapsed),
    x = map_x(388),
    y = map_y(y_max * 0.82),
    gp = gpar(fontsize = 7.4, col = collapsed_col, lineheight = 0.95)
  )
}

draw_scatter_panel <- function(x0, y0, width, height, depth_df, single_copy_depth) {
  draw_panel_frame(x0, y0, width, height)

  left <- x0 + 0.065
  right <- x0 + width - 0.028
  bottom <- y0 + 0.060
  top <- y0 + height - 0.040
  plot_w <- right - left
  plot_h <- top - bottom

  x_min <- 2.2
  x_max <- 6.3
  y_min <- 0
  y_max <- 600

  map_x <- function(v) left + plot_w * (log10(v) - x_min) / (x_max - x_min)
  map_y <- function(v) bottom + plot_h * (v - y_min) / (y_max - y_min)

  for (tick in seq(0, y_max, by = 100)) {
    ty <- map_y(tick)
    grid.lines(
      unit(c(left, right), "npc"),
      unit(c(ty, ty), "npc"),
      gp = gpar(col = grid_col, lwd = 0.6)
    )
    grid.lines(
      unit(c(left - 0.006, left), "npc"),
      unit(c(ty, ty), "npc"),
      gp = gpar(col = axis_col, lwd = 0.8)
    )
    grid.text(
      fmt_int(tick),
      x = left - 0.010,
      y = ty,
      just = c("right", "center"),
      gp = gpar(fontsize = 7.5, col = line_col)
    )
  }

  x_ticks <- c(1000, 10000, 100000, 1000000)
  x_tick_labels <- c("1 kb", "10 kb", "100 kb", "1 Mb")
  for (i in seq_along(x_ticks)) {
    tx <- map_x(x_ticks[[i]])
    grid.lines(
      unit(c(tx, tx), "npc"),
      unit(c(bottom, bottom - 0.006), "npc"),
      gp = gpar(col = axis_col, lwd = 0.8)
    )
    grid.text(
      x_tick_labels[[i]],
      x = tx,
      y = bottom - 0.020,
      gp = gpar(fontsize = 7.5, col = line_col)
    )
  }

  draw_plot_axes(left, bottom, right, top)

  y_single <- map_y(single_copy_depth)
  grid.lines(
    unit(c(left, right), "npc"),
    unit(c(y_single, y_single), "npc"),
    gp = gpar(col = line_col, lwd = 1.0, lty = "dashed")
  )
  grid.text(
    paste0(fmt2(single_copy_depth), "x"),
    x = right - 0.003,
    y = y_single + 0.012,
    just = c("right", "bottom"),
    gp = gpar(fontsize = 7.5, col = line_col)
  )

  category_levels <- c(
    "primary-like",
    "haplotig-like",
    "very-low",
    ">1.5C repeat-like",
    ">2C collapsed-repeat"
  )
  category_cols <- c(primary_col, haplotig_col, verylow_col, repeat_col, collapsed_col)
  names(category_cols) <- category_levels
  category_sizes <- c(0.115, 0.122, 0.128, 0.136, 0.145)
  names(category_sizes) <- category_levels

  for (category in category_levels) {
    rows <- depth_df[depth_df$category == category, ]
    if (nrow(rows) == 0) next
    grid.points(
      x = unit(map_x(rows$length_bp), "npc"),
      y = unit(map_y(rows$depth_x), "npc"),
      pch = 21,
      size = unit(category_sizes[[category]], "char"),
      gp = gpar(
        col = category_cols[[category]],
        fill = adjustcolor(category_cols[[category]], alpha.f = 0.55),
        lwd = 0.7
      )
    )
  }

  grid.text(
    "Contig length",
    x = left + plot_w / 2,
    y = y0 + 0.018,
    gp = gpar(fontsize = 8.7, col = line_col)
  )
  grid.text(
    "Mean depth (x)",
    x = x0 + 0.020,
    y = bottom + plot_h / 2,
    rot = 90,
    gp = gpar(fontsize = 8.7, col = line_col)
  )

  legend_x <- left + 0.015
  legend_y <- top - 0.014
  legend_items <- c("primary", "haplotig", "very-low", ">1.5C", ">2C")
  legend_keys <- category_levels
  grid.roundrect(
    x = legend_x + 0.030,
    y = legend_y - 0.028,
    width = 0.068,
    height = 0.068,
    r = unit(0.006, "snpc"),
    gp = gpar(fill = adjustcolor("white", alpha.f = 0.92), col = "#B8B8B8", lwd = 0.7)
  )
  for (i in seq_along(legend_items)) {
    lx <- legend_x
    ly <- legend_y - (i - 1) * 0.014
    grid.points(
      x = unit(lx, "npc"),
      y = unit(ly, "npc"),
      pch = 21,
      size = unit(0.078, "char"),
      gp = gpar(
        col = category_cols[[legend_keys[[i]]]],
        fill = adjustcolor(category_cols[[legend_keys[[i]]]], alpha.f = 0.55),
        lwd = 0.7
      )
    )
    grid.text(
      legend_items[[i]],
      x = lx + 0.009,
      y = ly,
      just = c("left", "center"),
      gp = gpar(fontsize = 6.5, col = line_col)
    )
  }

  contig160 <- depth_df[depth_df$contig == "contig_160", ]
  if (nrow(contig160) == 1) {
    px <- map_x(contig160$length_bp)
    py <- map_y(contig160$depth_x)
    label_x <- px + 0.020
    label_y <- min(py + 0.060, top - 0.010)
    grid.lines(
      unit(c(label_x - 0.004, px), "npc"),
      unit(c(label_y - 0.012, py + 0.004), "npc"),
      gp = gpar(col = collapsed_col, lwd = 0.8)
    )
    grid.text(
      paste0("contig_160\n370 kb @ ", fmt1(contig160$depth_x), "x"),
      x = label_x,
      y = label_y,
      just = c("left", "center"),
      gp = gpar(fontsize = 7.3, col = collapsed_col, lineheight = 0.92)
    )
  }
}

draw_reduction_panel <- function(x0, y0, width, height, assemblies,
                                 removed_contigs, span_loss_mb) {
  draw_panel_frame(x0, y0, width, height)

  left <- x0 + 0.105
  right <- x0 + width - 0.025
  bottom <- y0 + 0.130
  top <- y0 + height - 0.020
  plot_w <- right - left
  plot_h <- top - bottom
  x_max <- 110

  map_x <- function(v) left + plot_w * v / x_max

  for (tick in seq(0, x_max, by = 25)) {
    tx <- map_x(tick)
    grid.lines(
      unit(c(tx, tx), "npc"),
      unit(c(bottom, top), "npc"),
      gp = gpar(col = grid_col, lwd = 0.6)
    )
    grid.lines(
      unit(c(tx, tx), "npc"),
      unit(c(bottom, bottom - 0.006), "npc"),
      gp = gpar(col = axis_col, lwd = 0.8)
    )
    grid.text(
      fmt_int(tick),
      x = tx,
      y = bottom - 0.020,
      gp = gpar(fontsize = 7.5, col = line_col)
    )
  }

  draw_plot_axes(left, bottom, right, top)

  bar_y <- seq(top - 0.035, bottom + 0.055, length.out = nrow(assemblies))
  fills <- c(haplotig_col, primary_col, final_col)
  bar_h <- 0.032
  for (i in seq_len(nrow(assemblies))) {
    row <- assemblies[i, ]
    grid.text(
      row$label,
      x = left - 0.012,
      y = bar_y[[i]],
      just = c("right", "center"),
      gp = gpar(fontsize = 8.0, col = line_col)
    )
    grid.rect(
      x = left,
      y = bar_y[[i]],
      width = map_x(row$span_mb) - left,
      height = bar_h,
      just = c("left", "center"),
      gp = gpar(fill = fills[[i]], col = NA)
    )
    grid.text(
      paste0(fmt2(row$span_mb), " Mb / ", fmt_int(row$contigs), " contigs"),
      x = map_x(row$span_mb) - 0.004,
      y = bar_y[[i]],
      just = c("right", "center"),
      gp = gpar(fontsize = 7.7, col = "white", fontface = "bold")
    )
  }

  grid.text(
    "Assembly span (Mb)",
    x = left + plot_w / 2,
    y = y0 + 0.086,
    gp = gpar(fontsize = 8.7, col = line_col)
  )

  grid.roundrect(
    x = left + plot_w / 2,
    y = y0 + 0.044,
    width = 0.270,
    height = 0.050,
    r = unit(0.008, "snpc"),
    gp = gpar(fill = note_fill, col = "#D9B36A", lwd = 0.8)
  )
  grid.text(
    paste0(
      "Haplotig reduction: -", fmt_int(removed_contigs),
      " contigs, -", fmt2(span_loss_mb), " Mb"
    ),
    x = left + plot_w / 2,
    y = y0 + 0.044,
    gp = gpar(fontsize = 7.9, col = line_col)
  )
}

draw_genome_size_panel <- function(x0, y0, width, height, assembly_range,
                                   depth_estimate, flow_center, flow_sd) {
  draw_panel_frame(x0, y0, width, height)

  left <- x0 + 0.105
  right <- x0 + width - 0.030
  bottom <- y0 + 0.095
  top <- y0 + height - 0.040
  plot_w <- right - left
  x_min <- 80
  x_max <- 140
  map_x <- function(v) left + plot_w * (v - x_min) / (x_max - x_min)

  for (tick in seq(x_min, x_max, by = 10)) {
    tx <- map_x(tick)
    grid.lines(
      unit(c(tx, tx), "npc"),
      unit(c(bottom, top), "npc"),
      gp = gpar(col = grid_col, lwd = 0.6)
    )
    grid.lines(
      unit(c(tx, tx), "npc"),
      unit(c(bottom, bottom - 0.006), "npc"),
      gp = gpar(col = axis_col, lwd = 0.8)
    )
    if (tick %% 20 == 0) {
      grid.text(
        fmt_int(tick),
        x = tx,
        y = bottom - 0.020,
        gp = gpar(fontsize = 7.5, col = line_col)
      )
    }
  }

  draw_plot_axes(left, bottom, right, top)

  rows <- data.frame(
    label = c("Assembly-supported span", "Read-depth estimate", "Flow cytometry comparator"),
    y = c(top - 0.038, bottom + (top - bottom) * 0.52, bottom + 0.038),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(rows))) {
    grid.text(
      rows$label[[i]],
      x = left - 0.012,
      y = rows$y[[i]],
      just = c("right", "center"),
      gp = gpar(fontsize = 7.8, col = line_col)
    )
  }

  y1 <- rows$y[[1]]
  grid.lines(
    unit(c(map_x(assembly_range[1]), map_x(assembly_range[2])), "npc"),
    unit(c(y1, y1), "npc"),
    gp = gpar(col = primary_col, lwd = 5, lineend = "round")
  )
  grid.points(
    unit(map_x(assembly_range), "npc"),
    unit(c(y1, y1), "npc"),
    pch = 21,
    size = unit(0.070, "char"),
    gp = gpar(col = primary_col, fill = "white", lwd = 1.0)
  )
  grid.text(
    paste0(fmt2(assembly_range[1]), "-", fmt2(assembly_range[2]), " Mb"),
    x = map_x(assembly_range[2]) + 0.010,
    y = y1,
    just = c("left", "center"),
    gp = gpar(fontsize = 7.7, col = primary_col)
  )

  y2 <- rows$y[[2]]
  grid.points(
    unit(map_x(depth_estimate), "npc"),
    unit(y2, "npc"),
    pch = 21,
    size = unit(0.090, "char"),
    gp = gpar(col = line_col, fill = note_fill, lwd = 1.1)
  )
  grid.text(
    paste0(fmt2(depth_estimate), " Mb"),
    x = map_x(depth_estimate) + 0.011,
    y = y2,
    just = c("left", "center"),
    gp = gpar(fontsize = 7.7, col = line_col)
  )

  y3 <- rows$y[[3]]
  grid.lines(
    unit(c(map_x(flow_center - flow_sd), map_x(flow_center + flow_sd)), "npc"),
    unit(c(y3, y3), "npc"),
    gp = gpar(col = collapsed_col, lwd = 3.0, lineend = "round")
  )
  grid.points(
    unit(map_x(flow_center), "npc"),
    unit(y3, "npc"),
    pch = 21,
    size = unit(0.085, "char"),
    gp = gpar(col = collapsed_col, fill = adjustcolor(collapsed_col, alpha.f = 0.25), lwd = 1.0)
  )
  grid.text(
    paste0(fmt1(flow_center), " +/- ", fmt1(flow_sd), " Mb"),
    x = map_x(flow_center) - 0.011,
    y = y3 - 0.022,
    just = c("right", "center"),
    gp = gpar(fontsize = 7.7, col = collapsed_col)
  )

  grid.text(
    "Genome span (Mb)",
    x = left + plot_w / 2,
    y = y0 + 0.055,
    gp = gpar(fontsize = 8.7, col = line_col)
  )

  grid.roundrect(
    x = x0 + width / 2,
    y = y0 + 0.025,
    width = width - 0.050,
    height = 0.040,
    r = unit(0.008, "snpc"),
    gp = gpar(fill = note_fill, col = "#D9B36A", lwd = 0.8)
  )
  grid.text(
    "Primary span is supported near 100-109 Mb; the remaining difference suggests unresolved or collapsed repeats.",
    x = x0 + width / 2,
    y = y0 + 0.025,
    gp = gpar(fontsize = 7.7, col = line_col)
  )
}

depth_df <- read_depth_table("contig_mean_depth.tsv")
haplotig_df <- read_depth_table("haplotig_like.tsv")
verylow_df <- read_depth_table("verylow.tsv")
repeat15_df <- read_depth_table("repeat_gt1p5C.tsv")
collapsed_df <- read_depth_table("collapsed.repeat.tsv")

single_copy_depth <- 102.56
basecalled_bases <- 11161804683
estimated_genome_mb <- basecalled_bases / single_copy_depth / 1e6
repeat_threshold <- single_copy_depth * 1.5
collapsed_threshold <- single_copy_depth * 2

backbone <- assembly_stats("len15k backbone", "flye_len15k/assembly.fasta")
nohap <- assembly_stats("haplotig-reduced draft", "flye_len15k/assembly.no_haplotig.fasta")
final_ref <- assembly_stats("final annotation reference", "dorado_polish_gpu1/draft.polish2.bs6.fasta.fai")

assemblies <- rbind(backbone, nohap, final_ref)
span_loss_mb <- backbone$span_mb - nohap$span_mb
removed_contigs <- backbone$contigs - nohap$contigs

depth_df$category <- "primary-like"
depth_df$category[depth_df$contig %in% haplotig_df$contig] <- "haplotig-like"
depth_df$category[depth_df$contig %in% verylow_df$contig] <- "very-low"
depth_df$category[depth_df$contig %in% repeat15_df$contig] <- ">1.5C repeat-like"
depth_df$category[depth_df$contig %in% collapsed_df$contig] <- ">2C collapsed-repeat"

counts <- list(
  verylow = count_lines("verylow.tsv"),
  haplotig = count_lines("haplotig_like.tsv"),
  repeat15 = count_lines("repeat_gt1p5C.tsv"),
  collapsed = count_lines("collapsed.repeat.tsv")
)

flow_cytometry_center_mb <- 133
flow_cytometry_sd_mb <- 4
flow_cytometry_source <- "Zubacova et al. 2008, Mol Biochem Parasitol; doi:10.1016/j.molbiopara.2008.06.004"

source_metrics <- rbind(
  backbone,
  nohap,
  final_ref
)
write.table(
  source_metrics,
  paste0(output_prefix, "_source_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

depth_metrics <- data.frame(
  metric = c(
    "single_copy_depth_x",
    "basecalled_bases",
    "estimated_genome_size_mb",
    "haplotig_like_contigs",
    "very_low_depth_contigs",
    "repeat_gt_1p5C_contigs",
    "collapsed_repeat_gt_2C_contigs",
    "haplotig_reduction_span_loss_mb",
    "flow_cytometry_comparator_mb",
    "flow_cytometry_comparator_sd_mb"
  ),
  value = c(
    single_copy_depth,
    basecalled_bases,
    estimated_genome_mb,
    counts$haplotig,
    counts$verylow,
    counts$repeat15,
    counts$collapsed,
    span_loss_mb,
    flow_cytometry_center_mb,
    flow_cytometry_sd_mb
  ),
  source = c(
    "manifest and methods proxy depth",
    "methods final dataset",
    "basecalled_bases / single_copy_depth_x",
    "haplotig_like.tsv",
    "verylow.tsv",
    "repeat_gt1p5C.tsv",
    "collapsed.repeat.tsv",
    "flye_len15k/assembly.fasta vs assembly.no_haplotig.fasta",
    flow_cytometry_source,
    flow_cytometry_source
  ),
  stringsAsFactors = FALSE
)
write.table(
  depth_metrics,
  paste0(output_prefix, "_depth_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  depth_df,
  paste0(output_prefix, "_contig_depth_classes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  grid.text(
    "Figure 3. Read-depth profiling, repeat-collapse diagnosis, and genome-size interpretation",
    x = 0.5,
    y = 0.966,
    gp = gpar(fontsize = 22, fontface = "bold", col = line_col)
  )
  grid.text(
    expression(
      paste(
        "Data-derived validation of the selected ",
        italic("T. tenax"),
        " backbone, haplotig reduction, and residual high-depth repeat signal"
      )
    ),
    x = 0.5,
    y = 0.936,
    gp = gpar(fontsize = 10.8, col = muted_col)
  )

  draw_panel_title("A", "Contig mean-depth distribution", 0.045, 0.888)
  draw_panel_title("B", "Depth versus contig length", 0.535, 0.888)
  draw_panel_title("C", "Effect of haplotig reduction", 0.045, 0.488)
  draw_panel_title("D", "Genome-size interpretation", 0.535, 0.488)

  draw_histogram_panel(
    x0 = 0.055,
    y0 = 0.555,
    width = 0.400,
    height = 0.300,
    depths = depth_df$depth_x,
    haplotig_range = range(haplotig_df$depth_x),
    single_copy_depth = single_copy_depth,
    repeat_threshold = repeat_threshold,
    collapsed_threshold = collapsed_threshold,
    counts = counts
  )

  draw_scatter_panel(
    x0 = 0.545,
    y0 = 0.555,
    width = 0.400,
    height = 0.300,
    depth_df = depth_df,
    single_copy_depth = single_copy_depth
  )

  draw_reduction_panel(
    x0 = 0.055,
    y0 = 0.155,
    width = 0.400,
    height = 0.300,
    assemblies = assemblies,
    removed_contigs = removed_contigs,
    span_loss_mb = span_loss_mb
  )

  draw_genome_size_panel(
    x0 = 0.545,
    y0 = 0.155,
    width = 0.400,
    height = 0.300,
    assembly_range = range(c(backbone$span_mb, nohap$span_mb, final_ref$span_mb)),
    depth_estimate = estimated_genome_mb,
    flow_center = flow_cytometry_center_mb,
    flow_sd = flow_cytometry_sd_mb
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
message("Wrote: ", paste0(output_prefix, "_contig_depth_classes.tsv"))
