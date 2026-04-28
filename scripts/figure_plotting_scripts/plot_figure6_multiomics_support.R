#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args) >= 1) {
  args[[1]]
} else {
  "figures/Tt_Figure6_multiomics_support_v1"
}

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

library(grid)

line_col <- "#111111"
axis_col <- "#303030"
muted_col <- "#5C5C5C"
grid_col <- "#E7E7E7"
panel_fill <- "#FFFFFF"
drs_col <- "#2E8B7C"
sr_col <- "#3A7CA5"
lr_col <- "#7A5AA6"
rna_col <- "#4F88B6"
any_col <- "#2F7F63"
all_col <- "#6C6C6C"
prot_col <- "#D18422"
class_col <- "#4A90A4"
none_col <- "#B8B8B8"

paths <- list(
  matrix = "Express/multiomics_gene_evidence_matrix.tsv",
  proteomics = "Express/proteomics_gene_level.tsv",
  drs_flagstat = "Express/qc/Tt_DRS_merged.hac.genome.flagstat.txt",
  lr_flagstat = "Express/LRcDNA_qc/Tt_cDNA_merged.sup.Q10.min200.max50k.genome.flagstat.txt",
  sr_star_log = "Express/SR_aln/Tt_HVFMLDSX2.Log.final.out"
)

missing_paths <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing_paths) > 0) {
  stop("Missing required input files: ", paste(missing_paths, collapse = ", "))
}

fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_pct1 <- function(x) sprintf("%.1f%%", x)
fmt_count_pct <- function(n, total) paste0(fmt_int(n), "\n(", fmt_pct1(n / total * 100), ")")
yes <- function(x) x == "yes"

wrap_text <- function(x, width = 34) {
  paste(strwrap(x, width = width), collapse = "\n")
}

parse_flagstat_rate <- function(path, prefer_primary = FALSE) {
  lines <- readLines(path, warn = FALSE)
  extract_first_int <- function(line) as.numeric(sub("^\\s*([0-9]+).*", "\\1", line))
  extract_pct <- function(line) {
    pct <- sub(".*\\(([0-9.]+)%.*", "\\1", line)
    suppressWarnings(as.numeric(pct))
  }

  if (prefer_primary && any(grepl("primary mapped", lines, fixed = TRUE))) {
    primary_line <- grep(" primary$", lines, value = TRUE)
    primary_mapped_line <- grep("primary mapped", lines, value = TRUE)[1]
    total <- extract_first_int(primary_line[1])
    mapped <- extract_first_int(primary_mapped_line)
    pct <- extract_pct(primary_mapped_line)
  } else {
    total_line <- grep(" in total ", lines, value = TRUE)[1]
    mapped_line <- grep(" mapped \\(", lines, value = TRUE)[1]
    total <- extract_first_int(total_line)
    mapped <- extract_first_int(mapped_line)
    pct <- extract_pct(mapped_line)
  }

  list(mapped = mapped, total = total, percent = pct)
}

parse_star_unique_rate <- function(path) {
  lines <- readLines(path, warn = FALSE)
  value_for <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE, fixed = TRUE)[1]
    sub(".*\\|\\s*", "", line)
  }
  input <- as.numeric(value_for("Number of input reads"))
  unique <- as.numeric(value_for("Uniquely mapped reads number"))
  pct <- as.numeric(sub("%", "", value_for("Uniquely mapped reads %"), fixed = TRUE))
  list(mapped = unique, total = input, percent = pct)
}

draw_panel_title <- function(letter, title, x, y) {
  grid.text(letter, x = x, y = y, just = c("left", "center"),
            gp = gpar(fontsize = 20, fontface = "bold", col = line_col))
  grid.text(title, x = x + 0.035, y = y, just = c("left", "center"),
            gp = gpar(fontsize = 15, fontface = "bold", col = line_col))
}

draw_panel_frame <- function(x0, y0, width, height) {
  grid.rect(x = x0, y = y0, width = width, height = height,
            just = c("left", "bottom"),
            gp = gpar(fill = panel_fill, col = line_col, lwd = 1.0))
}

draw_axes <- function(left, bottom, right, top) {
  grid.lines(unit(c(left, right), "npc"), unit(c(bottom, bottom), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))
  grid.lines(unit(c(left, left), "npc"), unit(c(bottom, top), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))
}

draw_vertical_bars <- function(x0, y0, width, height, df, y_max, y_step,
                               y_label, label_mode = c("percent", "count_pct")) {
  label_mode <- match.arg(label_mode)
  draw_panel_frame(x0, y0, width, height)

  left <- x0 + 0.070
  right <- x0 + width - 0.030
  bottom <- y0 + 0.075
  top <- y0 + height - 0.055
  plot_w <- right - left
  plot_h <- top - bottom
  map_y <- function(v) bottom + plot_h * v / y_max

  for (tick in seq(0, y_max, by = y_step)) {
    ty <- map_y(tick)
    grid.lines(unit(c(left, right), "npc"), unit(c(ty, ty), "npc"),
               gp = gpar(col = grid_col, lwd = 0.6))
    grid.lines(unit(c(left - 0.006, left), "npc"), unit(c(ty, ty), "npc"),
               gp = gpar(col = axis_col, lwd = 0.8))
    grid.text(fmt_int(tick), x = left - 0.010, y = ty,
              just = c("right", "center"),
              gp = gpar(fontsize = 8.0, col = line_col))
  }
  draw_axes(left, bottom, right, top)

  centers <- left + plot_w * (seq_len(nrow(df)) - 0.5) / nrow(df)
  bar_w <- plot_w / nrow(df) * 0.58
  for (i in seq_len(nrow(df))) {
    value <- df$value[[i]]
    grid.rect(x = centers[[i]], y = bottom, width = bar_w,
              height = map_y(value) - bottom,
              just = c("center", "bottom"),
              gp = gpar(fill = df$color[[i]], col = NA))
    label <- if (label_mode == "percent") {
      fmt_pct1(value)
    } else {
      paste0(fmt_int(value), "\n(", fmt_pct1(value / df$total[[i]] * 100), ")")
    }
    grid.text(label, x = centers[[i]], y = map_y(value) + 0.018,
              gp = gpar(fontsize = 8.6, col = line_col, lineheight = 0.95))
    grid.text(df$label[[i]], x = centers[[i]], y = bottom - 0.037,
              gp = gpar(fontsize = 8.7, col = line_col, lineheight = 0.92))
  }

  grid.text(y_label, x = left - 0.050, y = bottom + plot_h / 2,
            rot = 90, gp = gpar(fontsize = 9.5, col = line_col))
}

draw_horizontal_bars <- function(x0, y0, width, height, df, x_max, x_step,
                                 x_label, label_total = NULL, label_font = 8.4) {
  draw_panel_frame(x0, y0, width, height)

  left <- x0 + 0.145
  right <- x0 + width - 0.030
  bottom <- y0 + 0.060
  top <- y0 + height - 0.045
  plot_w <- right - left
  plot_h <- top - bottom
  n <- nrow(df)
  map_x <- function(v) left + plot_w * v / x_max

  for (tick in seq(0, x_max, by = x_step)) {
    tx <- map_x(tick)
    grid.lines(unit(c(tx, tx), "npc"), unit(c(bottom, top), "npc"),
               gp = gpar(col = grid_col, lwd = 0.6))
    grid.text(fmt_int(tick), x = tx, y = bottom - 0.030,
              gp = gpar(fontsize = 7.6, col = line_col))
  }
  grid.lines(unit(c(left, right), "npc"), unit(c(bottom, bottom), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))
  grid.lines(unit(c(left, left), "npc"), unit(c(bottom, top), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))

  row_h <- plot_h / n
  for (i in seq_len(n)) {
    y <- top - row_h * (i - 0.5)
    value <- df$value[[i]]
    fill <- df$color[[i]]
    grid.rect(x = left, y = y, width = map_x(value) - left,
              height = row_h * 0.55,
              just = c("left", "center"),
              gp = gpar(fill = fill, col = NA))
    grid.text(df$label[[i]], x = left - 0.013, y = y,
              just = c("right", "center"),
              gp = gpar(fontsize = 8.2, col = line_col, lineheight = 0.92))
    pct_text <- if (!is.null(label_total)) {
      paste0(" (", fmt_pct1(value / label_total * 100), ")")
    } else {
      ""
    }
    grid.text(paste0(fmt_int(value), pct_text),
              x = min(map_x(value) + 0.007, right - 0.006), y = y,
              just = c("left", "center"),
              gp = gpar(fontsize = label_font, col = line_col))
  }

  grid.text(x_label, x = left + plot_w / 2, y = y0 + 0.020,
            gp = gpar(fontsize = 8.6, col = line_col))
}

draw_inset_note <- function(x, y, width, height, title, body) {
  grid.roundrect(x = x, y = y, width = width, height = height,
                 r = unit(0.010, "snpc"),
                 gp = gpar(fill = "#FFF8EC", col = line_col, lwd = 0.9))
  grid.text(title, x = x, y = y + height / 2 - 0.026,
            gp = gpar(fontsize = 9.2, fontface = "bold", col = line_col))
  grid.text(wrap_text(body, 58), x = x, y = y - 0.008,
            gp = gpar(fontsize = 7.3, col = line_col, lineheight = 0.95))
}

matrix <- read.delim(paths$matrix, stringsAsFactors = FALSE, quote = "", comment.char = "")
proteomics <- read.delim(paths$proteomics, stringsAsFactors = FALSE, quote = "", comment.char = "")

total_genes <- nrow(matrix)

drs_rate <- parse_flagstat_rate(paths$drs_flagstat, prefer_primary = FALSE)
sr_rate <- parse_star_unique_rate(paths$sr_star_log)
lr_rate <- parse_flagstat_rate(paths$lr_flagstat, prefer_primary = TRUE)

mapping_df <- data.frame(
  label = c("DRS", "SR cDNA", "LR cDNA"),
  value = c(drs_rate$percent, sr_rate$percent, lr_rate$percent),
  mapped = c(drs_rate$mapped, sr_rate$mapped, lr_rate$mapped),
  total = c(drs_rate$total, sr_rate$total, lr_rate$total),
  color = c(drs_col, sr_col, lr_col),
  source = c(paths$drs_flagstat, paths$sr_star_log, paths$lr_flagstat),
  stringsAsFactors = FALSE
)

rna_df <- data.frame(
  label = c("DRS", "SR cDNA", "LR cDNA", "Any RNA", "All 3 RNA"),
  value = c(
    sum(yes(matrix$DRS_support)),
    sum(yes(matrix$SR_cDNA_support)),
    sum(yes(matrix$LR_cDNA_support)),
    sum(yes(matrix$Any_RNA_support)),
    sum(yes(matrix$All_three_RNA_support))
  ),
  total = total_genes,
  color = c(drs_col, sr_col, lr_col, any_col, all_col),
  stringsAsFactors = FALSE
)

class_order <- c(
  "no_detected_omics_support",
  "single_RNA_supported",
  "multi_RNA_supported",
  "three_RNA_supported",
  "RNA_and_proteomics_supported",
  "three_RNA_and_proteomics_supported",
  "three_RNA_and_high_confidence_proteomics_supported"
)
class_labels <- c(
  "No detected\nomics support",
  "Single RNA",
  "Two RNA",
  "Three RNA only",
  "RNA + proteomics",
  "Three RNA +\nproteomics",
  "Three RNA +\n>=2 unique peptides"
)
class_counts <- table(factor(matrix$Evidence_class, levels = class_order))
class_df <- data.frame(
  label = class_labels,
  class = class_order,
  value = as.integer(class_counts),
  color = c(none_col, rna_col, rna_col, rna_col, prot_col, prot_col, prot_col),
  stringsAsFactors = FALSE
)

proteomics_df <- data.frame(
  label = c(
    ">=1 unique\npeptide",
    ">=2 unique\npeptides",
    ">=3 unique\npeptides",
    ">=5 unique\npeptides",
    "Three RNA +\nproteomics",
    "Three RNA +\n>=2 unique peptides"
  ),
  value = c(
    sum(yes(matrix$Proteomics_support)),
    sum(yes(matrix$High_confidence_proteomics_2_unique_peptides)),
    sum(yes(matrix$Strong_proteomics_3_unique_peptides)),
    sum(yes(matrix$Very_strong_proteomics_5_unique_peptides)),
    sum(yes(matrix$All_three_RNA_support) & yes(matrix$Proteomics_support)),
    sum(yes(matrix$All_three_RNA_support) & yes(matrix$High_confidence_proteomics_2_unique_peptides))
  ),
  color = c(prot_col, prot_col, prot_col, prot_col, "#8C6D31", "#8C6D31"),
  stringsAsFactors = FALSE
)

source_data <- rbind(
  data.frame(panel = "A", metric = "mapping_rate_percent", label = mapping_df$label,
             value = mapping_df$value, total = mapping_df$total, percent = mapping_df$value,
             source = mapping_df$source, stringsAsFactors = FALSE),
  data.frame(panel = "B", metric = "RNA_supported_genes", label = rna_df$label,
             value = rna_df$value, total = rna_df$total,
             percent = rna_df$value / rna_df$total * 100,
             source = paths$matrix, stringsAsFactors = FALSE),
  data.frame(panel = "C", metric = "mutually_exclusive_evidence_class",
             label = gsub("\n", " ", class_df$label),
             value = class_df$value, total = total_genes,
             percent = class_df$value / total_genes * 100,
             source = paths$matrix, stringsAsFactors = FALSE),
  data.frame(panel = "D", metric = "proteomics_gene_support",
             label = gsub("\n", " ", proteomics_df$label),
             value = proteomics_df$value, total = total_genes,
             percent = proteomics_df$value / total_genes * 100,
             source = paths$matrix, stringsAsFactors = FALSE)
)
write.table(source_data, paste0(output_prefix, "_source_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

proteomics_summary <- data.frame(
  metric = c(
    "proteomics_gene_level_rows",
    "proteomics_supported_genes_ge1_unique_peptide",
    "high_confidence_ge2_unique_peptides",
    "strong_ge3_unique_peptides",
    "very_strong_ge5_unique_peptides",
    "rna_and_proteomics_supported_genes",
    "three_rna_and_proteomics_supported_genes",
    "three_rna_and_ge2_unique_peptide_supported_genes"
  ),
  value = c(
    nrow(proteomics),
    proteomics_df$value[[1]],
    proteomics_df$value[[2]],
    proteomics_df$value[[3]],
    proteomics_df$value[[4]],
    sum(yes(matrix$Any_RNA_support) & yes(matrix$Proteomics_support)),
    proteomics_df$value[[5]],
    proteomics_df$value[[6]]
  ),
  source = c(paths$proteomics, rep(paths$matrix, 7)),
  stringsAsFactors = FALSE
)
write.table(proteomics_summary, paste0(output_prefix, "_proteomics_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  grid.text(expression(paste("Figure 6. Multi-omics support for the final ", italic("T. tenax"), " gene annotation")),
            x = 0.5, y = 0.966,
            gp = gpar(fontsize = 22, fontface = "bold", col = line_col))
  grid.text(paste0(
    "RNA and proteomics evidence across ", fmt_int(total_genes),
    " predicted genes; counts are annotation-support evidence, not differential expression"
  ),
  x = 0.5, y = 0.936,
  gp = gpar(fontsize = 10.4, col = muted_col))

  draw_panel_title("A", "RNA dataset alignment to final genome", 0.045, 0.888)
  draw_panel_title("B", "Gene support by RNA evidence", 0.535, 0.888)
  draw_panel_title("C", "Multi-omics evidence classes", 0.045, 0.488)
  draw_panel_title("D", "Proteomics-supported gene evidence", 0.535, 0.488)

  draw_vertical_bars(0.055, 0.555, 0.400, 0.300, mapping_df,
                     y_max = 105, y_step = 20,
                     y_label = "Mapping rate (%)", label_mode = "percent")
  draw_vertical_bars(0.545, 0.555, 0.400, 0.300, rna_df,
                     y_max = 25000, y_step = 5000,
                     y_label = "Supported genes", label_mode = "count_pct")
  draw_horizontal_bars(0.055, 0.155, 0.400, 0.300, class_df,
                       x_max = 20000, x_step = 5000,
                       x_label = "Number of genes", label_total = total_genes,
                       label_font = 8.0)
  draw_horizontal_bars(0.545, 0.155, 0.400, 0.300, proteomics_df,
                       x_max = 540, x_step = 100,
                       x_label = "Number of genes", label_total = NULL,
                       label_font = 8.2)

  grid.text(
    "DRS, direct RNA sequencing; SR cDNA, short-read cDNA; LR cDNA, Nanopore long-read cDNA. Primary count mode: featureCounts -s 0.",
    x = 0.5, y = 0.035,
    gp = gpar(fontsize = 8.0, col = muted_col)
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
message("Wrote: ", paste0(output_prefix, "_source_data.tsv"))
message("Wrote: ", paste0(output_prefix, "_proteomics_summary.tsv"))
