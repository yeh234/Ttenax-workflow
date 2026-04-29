#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args) >= 1) {
  args[[1]]
} else {
  "figures/Tt_Figure4_annotation_completeness_v1"
}

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

library(grid)

line_col <- "#111111"
muted_col <- "#5C5C5C"
axis_col <- "#303030"
grid_col <- "#E6E6E6"
panel_fill <- "#FFFFFF"
note_fill <- "#FFF8EC"
struct_col <- "#3A7CA5"
eggnog_col <- "#D18422"
tv_col <- "#2E8B7C"
ips_col <- "#7A5AA6"
missing_col <- "#C7C7C7"
frag_col <- "#D9B36A"
none_col <- "#B65C4A"

paths <- list(
  galba_gff = "dorado_polish_gpu1/galba_out/galba.gff3",
  proteins = "dorado_polish_gpu1/Tt_predicted_proteins.aa.fasta",
  cds = "dorado_polish_gpu1/Tt_predicted_cds.nt.fasta",
  all_ids = "dorado_polish_gpu1/eggnogHit/all_ids.txt",
  eggnog_ids = "dorado_polish_gpu1/eggnogHit/annotated_ids.txt",
  unannotated_ids = "dorado_polish_gpu1/eggnogHit/unannotated_ids.txt",
  tv_blast_all = "dorado_polish_gpu1/Tt_vs_Tv_annotation.txt",
  tv_rescue = "dorado_polish_gpu1/eggnogHit/unannotated_vs_Tv_with_desc.txt",
  tv_mapping_gff = "final_Tt/Tt_alignmentTrichDBv68.gff",
  tv_reference_proteins = "dorado_polish_gpu1/galba_out/prot_seq_all.fa",
  busco_miniprot = "dorado_polish_gpu1/busco_odb12_p2/short_summary.specific.eukaryota_odb12.busco_odb12_p2.txt",
  busco_metaeuk = "dorado_polish_gpu1/busco_odb12_p2_metaeuk/short_summary.specific.eukaryota_odb12.busco_odb12_p2_metaeuk.txt",
  interpro_tsv = "dorado_polish_gpu1/interproscan_results/Tt_predicted_proteins.aa.fasta.tsv",
  final_table = "final_Tt/Tt_Comprehensive_Annotation_Table.tsv"
)

missing_paths <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing_paths) > 0) {
  stop("Missing required input files: ", paste(missing_paths, collapse = ", "))
}

fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt1 <- function(x) sprintf("%.1f", x)
fmt_pct <- function(n, d) sprintf("%.1f%%", n / d * 100)
present <- function(x) !is.na(x) & nzchar(x) & x != "-"

stream_count_fasta <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  total <- 0L
  repeat {
    lines <- readLines(con, n = 20000, warn = FALSE)
    if (length(lines) == 0) break
    total <- total + sum(startsWith(lines, ">"))
  }
  total
}

stream_count_gff_feature <- function(path, feature) {
  con <- file(path, "rt")
  on.exit(close(con), add = TRUE)
  total <- 0L
  repeat {
    lines <- readLines(con, n = 20000, warn = FALSE)
    if (length(lines) == 0) break
    lines <- lines[!startsWith(lines, "#")]
    if (length(lines) == 0) next
    fields <- strsplit(lines, "\t", fixed = TRUE)
    total <- total + sum(vapply(fields, function(x) length(x) >= 3 && x[[3]] == feature, logical(1)))
  }
  total
}

stream_unique_first_col <- function(path) {
  con <- file(path, "rt")
  on.exit(close(con), add = TRUE)
  ids <- character()
  repeat {
    lines <- readLines(con, n = 2000, warn = FALSE)
    if (length(lines) == 0) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (length(lines) == 0) next
    ids <- unique(c(ids, sub("\t.*$", "", lines)))
  }
  ids
}

read_busco_summary <- function(path, label) {
  lines <- readLines(path, warn = FALSE)
  summary_line <- grep("C:[0-9.]+%", lines, value = TRUE)[1]
  pct <- function(prefix) {
    as.numeric(sub(paste0(".*", prefix, ":([0-9.]+)%.*"), "\\1", summary_line))
  }
  count_for <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE)[1]
    as.integer(sub("^\\s*([0-9]+).*", "\\1", line))
  }
  internal_line <- grep("internal stop codons", lines, value = TRUE)
  internal_stop <- if (length(internal_line) == 0) {
    NA_integer_
  } else {
    as.integer(sub(".*of which ([0-9]+).*", "\\1", internal_line[[1]]))
  }
  data.frame(
    workflow = label,
    complete_pct = pct("C"),
    single_pct = pct("S"),
    duplicated_pct = pct("D"),
    fragmented_pct = pct("F"),
    missing_pct = pct("M"),
    complete_count = count_for("Complete BUSCOs \\(C\\)"),
    single_count = count_for("Complete and single-copy BUSCOs"),
    duplicated_count = count_for("Complete and duplicated BUSCOs"),
    fragmented_count = count_for("Fragmented BUSCOs"),
    missing_count = count_for("Missing BUSCOs"),
    total_buscos = count_for("Total BUSCO groups searched"),
    internal_stop_complete = internal_stop,
    stringsAsFactors = FALSE
  )
}

read_blast_unique_queries <- function(path, normalize_dot = FALSE) {
  tab <- read.delim(path, header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  ids <- tab[[1]]
  if (normalize_dot) ids <- gsub("\\.", "", ids)
  list(rows = nrow(tab), unique_queries = length(unique(ids)), ids = unique(ids))
}

stream_tv_mapping <- function(path) {
  con <- file(path, "rt")
  on.exit(close(con), add = TRUE)
  mrna_records <- 0L
  targets <- character()
  high_identity_records <- 0L
  high_identity_targets <- character()
  repeat {
    lines <- readLines(con, n = 20000, warn = FALSE)
    if (length(lines) == 0) break
    lines <- lines[grepl("\tmRNA\t", lines, fixed = TRUE)]
    if (length(lines) == 0) next
    mrna_records <- mrna_records + length(lines)
    target <- sub(".*Target=([^; ]+).*", "\\1", lines)
    identity_txt <- sub(".*Identity=([0-9.]+).*", "\\1", lines)
    identity <- suppressWarnings(as.numeric(identity_txt))
    targets <- unique(c(targets, target))
    high <- !is.na(identity) & identity >= 0.6
    high_identity_records <- high_identity_records + sum(high)
    high_identity_targets <- unique(c(high_identity_targets, target[high]))
  }
  list(
    mrna_records = mrna_records,
    unique_targets = length(unique(targets)),
    high_identity_records = high_identity_records,
    high_identity_unique_targets = length(unique(high_identity_targets))
  )
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
            gp = gpar(fill = panel_fill, col = NA, lwd = 1.0))
}

draw_arrow <- function(x0, y0, x1, y1, lwd = 1.4) {
  grid.lines(unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
             gp = gpar(col = line_col, lwd = lwd),
             arrow = arrow(length = unit(0.10, "inches"), type = "closed"))
}

draw_pipeline_panel <- function(x0, y0, width, height, metrics) {
  draw_panel_frame(x0, y0, width, height)
  boxes <- list(
    list(
      title = "GALBA structural prediction",
      detail = paste0(fmt_int(metrics$galba_genes), " genes; ", fmt_int(metrics$galba_mrna), " mRNAs")
    ),
    list(
      title = "gffread extraction",
      detail = paste0(fmt_int(metrics$proteins), " proteins and CDS sequences")
    ),
    list(
      title = "EggNOG-mapper",
      detail = paste0(fmt_int(metrics$eggnog_annotated), " annotated proteins (",
                      fmt_pct(metrics$eggnog_annotated, metrics$proteins), ")")
    ),
    list(
      title = "Tv BLAST rescue",
      detail = paste0(fmt_int(metrics$tv_rescue_unique), " / ",
                      fmt_int(metrics$unannotated), " EggNOG-unannotated hit")
    ),
    list(
      title = "InterProScan domain support",
      detail = paste0(fmt_int(metrics$interproscan_unique), " proteins with records")
    ),
    list(
      title = "Evidence merge",
      detail = "EggNOG + Tv rescue + InterPro + GO/KEGG"
    ),
    list(
      title = "Final master table",
      detail = paste0(fmt_int(metrics$final_rows), " transcript-level rows")
    )
  )

  box_w <- width - 0.060
  box_h <- 0.0225
  gap <- 0.0160
  top_y <- y0 + height - 0.030
  cx <- x0 + width / 2
  ys <- top_y - seq(0, by = box_h + gap, length.out = length(boxes))

  for (i in seq_along(boxes)) {
    fill <- if (i %in% c(1, 2, 7)) "#F3F8FB" else "#FBFBFB"
    if (i == 4) fill <- "#F3FAF6"
    if (i == 5) fill <- "#F7F2FA"
    grid.roundrect(x = cx, y = ys[[i]], width = box_w, height = box_h,
                   r = unit(0.010, "snpc"),
                   gp = gpar(fill = fill, col = line_col, lwd = 1.2))
    grid.text(bquote(bold(.(boxes[[i]]$title)) ~ .(boxes[[i]]$detail)),
              x = cx, y = ys[[i]],
              gp = gpar(fontsize = 7.45, col = line_col))
    if (i < length(boxes)) {
      draw_arrow(cx, ys[[i]] - box_h / 2, cx, ys[[i + 1]] + box_h / 2)
    }
  }
}

draw_horizontal_bar_panel <- function(x0, y0, width, height, counts_df, total) {
  draw_panel_frame(x0, y0, width, height)
  left <- x0 + 0.120
  right <- x0 + width - 0.055
  bottom <- y0 + 0.055
  top <- y0 + height - 0.035
  plot_w <- right - left
  x_max <- 36000
  map_x <- function(v) left + plot_w * v / x_max

  for (tick in seq(0, x_max, by = 10000)) {
    tx <- map_x(tick)
    grid.lines(unit(c(tx, tx), "npc"), unit(c(bottom, top), "npc"),
               gp = gpar(col = grid_col, lwd = 0.6))
    grid.lines(unit(c(tx, tx), "npc"), unit(c(bottom, bottom - 0.006), "npc"),
               gp = gpar(col = axis_col, lwd = 0.8))
    grid.text(fmt_int(tick), x = tx, y = bottom - 0.020,
              gp = gpar(fontsize = 7.2, col = line_col))
  }
  grid.lines(unit(c(left, right), "npc"), unit(c(bottom, bottom), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))
  grid.lines(unit(c(left, left), "npc"), unit(c(bottom, top), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))

  ys <- seq(top - 0.022, bottom + 0.030, length.out = nrow(counts_df))
  bar_h <- min(0.025, (top - bottom) / nrow(counts_df) * 0.54)
  for (i in seq_len(nrow(counts_df))) {
    row <- counts_df[i, ]
    grid.text(row$label, x = left - 0.010, y = ys[[i]],
              just = c("right", "center"),
              gp = gpar(fontsize = 7.3, col = line_col))
    grid.rect(x = left, y = ys[[i]], width = map_x(row$count) - left, height = bar_h,
              just = c("left", "center"),
              gp = gpar(fill = row$color, col = NA))
    label <- paste0(fmt_int(row$count), " (", fmt_pct(row$count, total), ")")
    label_x <- min(map_x(row$count) + 0.006, right - 0.002)
    just <- if (label_x > right - 0.004) c("right", "center") else c("left", "center")
    grid.text(label, x = label_x, y = ys[[i]], just = just,
              gp = gpar(fontsize = 7.2, col = line_col))
  }
  grid.text("Protein or feature count", x = left + plot_w / 2, y = y0 + 0.018,
            gp = gpar(fontsize = 8.4, col = line_col))
}

draw_stacked_bar <- function(x0, y, width, height, values, colors, labels,
                             label_prefix = NULL) {
  start <- x0
  for (i in seq_along(values)) {
    w <- width * values[[i]] / 100
    grid.rect(x = start, y = y, width = w, height = height,
              just = c("left", "center"),
              gp = gpar(fill = colors[[i]], col = "white", lwd = 0.5))
    if (values[[i]] >= 8) {
      grid.text(labels[[i]], x = start + w / 2, y = y,
                gp = gpar(fontsize = 7.2, col = "white", fontface = "bold"))
    }
    start <- start + w
  }
  if (!is.null(label_prefix)) {
    grid.text(label_prefix, x = x0 - 0.010, y = y, just = c("right", "center"),
              gp = gpar(fontsize = 7.7, col = line_col))
  }
}

draw_completeness_panel <- function(x0, y0, width, height, busco_df,
                                    tv_mapping, tv_reference_total) {
  draw_panel_frame(x0, y0, width, height)
  left <- x0 + 0.105
  right <- x0 + width - 0.055
  bar_w <- right - left

  grid.text("BUSCO eukaryota_odb12 summary", x = x0 + width / 2, y = y0 + height - 0.027,
            gp = gpar(fontsize = 9.4, fontface = "bold", col = line_col))

  y_top <- y0 + height - 0.070
  bar_h <- 0.030
  for (i in seq_len(nrow(busco_df))) {
    row <- busco_df[i, ]
    y <- y_top - (i - 1) * 0.052
    draw_stacked_bar(
      left, y, bar_w, bar_h,
      c(row$complete_pct, row$fragmented_pct, row$missing_pct),
      c(tv_col, frag_col, missing_col),
      c(paste0("C ", fmt1(row$complete_pct), "%"),
        paste0("F ", fmt1(row$fragmented_pct), "%"),
        paste0("M ", fmt1(row$missing_pct), "%")),
      row$workflow
    )
  }

  grid.points(unit(c(left + 0.010, left + 0.115, left + 0.220), "npc"),
              unit(rep(y0 + height - 0.149, 3), "npc"),
              pch = 15, size = unit(0.10, "char"),
              gp = gpar(col = c(tv_col, frag_col, missing_col)))
  grid.text("Complete", x = left + 0.020, y = y0 + height - 0.149, just = c("left", "center"),
            gp = gpar(fontsize = 6.9, col = line_col))
  grid.text("Fragmented", x = left + 0.125, y = y0 + height - 0.149, just = c("left", "center"),
            gp = gpar(fontsize = 6.9, col = line_col))
  grid.text("Missing", x = left + 0.230, y = y0 + height - 0.149, just = c("left", "center"),
            gp = gpar(fontsize = 6.9, col = line_col))

  miniprot_stop <- busco_df$internal_stop_complete[busco_df$workflow == "Miniprot"]
  grid.text(paste0("Miniprot complete BUSCOs with internal stops: ",
                   fmt_int(miniprot_stop), " / ",
                   fmt_int(busco_df$complete_count[busco_df$workflow == "Miniprot"])),
            x = x0 + width / 2, y = y0 + height - 0.178,
            gp = gpar(fontsize = 7.3, col = muted_col))

  grid.text("Tv proteome-to-genome mapping", x = x0 + width / 2, y = y0 + 0.105,
            gp = gpar(fontsize = 9.4, fontface = "bold", col = line_col))
  map_y <- y0 + 0.066
  grid.rect(x = left, y = map_y, width = bar_w, height = 0.030,
            just = c("left", "center"),
            gp = gpar(fill = "#F1F1F1", col = "#BFBFBF", lwd = 0.6))
  grid.rect(x = left, y = map_y,
            width = bar_w * tv_mapping$unique_targets / tv_reference_total,
            height = 0.030, just = c("left", "center"),
            gp = gpar(fill = struct_col, col = NA))
  grid.text(
    paste0(fmt_int(tv_mapping$unique_targets), " / ", fmt_int(tv_reference_total),
           " unique Tv proteins (",
           fmt_pct(tv_mapping$unique_targets, tv_reference_total), ")"),
    x = left + bar_w / 2, y = map_y,
    gp = gpar(fontsize = 7.5, col = "white", fontface = "bold")
  )
  grid.text(
    paste0(fmt_int(tv_mapping$mrna_records), " mRNA alignment records; ",
           fmt_int(tv_mapping$high_identity_records), " alignments with Identity >= 0.6"),
    x = x0 + width / 2, y = y0 + 0.027,
    gp = gpar(fontsize = 7.1, col = muted_col)
  )
}

draw_naming_rescue_panel <- function(x0, y0, width, height, naming_df, rescue_df, total) {
  draw_panel_frame(x0, y0, width, height)
  grid.text("Final product naming sources", x = x0 + width * 0.225, y = y0 + height - 0.030,
            gp = gpar(fontsize = 9.0, fontface = "bold", col = line_col))
  grid.text("EggNOG-unannotated fraction", x = x0 + width * 0.735, y = y0 + height - 0.030,
            gp = gpar(fontsize = 9.0, fontface = "bold", col = line_col))

  stack_x <- x0 + 0.035
  stack_y <- y0 + height - 0.083
  stack_w <- width * 0.385
  stack_h <- 0.040
  start <- stack_x
  for (i in seq_len(nrow(naming_df))) {
    seg_w <- stack_w * naming_df$count[[i]] / total
    grid.rect(x = start, y = stack_y, width = seg_w, height = stack_h,
              just = c("left", "center"),
              gp = gpar(fill = naming_df$color[[i]], col = "white", lwd = 0.5))
    if (seg_w > 0.045) {
      grid.text(paste0(naming_df$short[[i]], "\n", fmt_int(naming_df$count[[i]])),
                x = start + seg_w / 2, y = stack_y,
                gp = gpar(fontsize = 6.8, col = "white", fontface = "bold", lineheight = 0.9))
    }
    start <- start + seg_w
  }

  legend_y <- y0 + height - 0.143
  for (i in seq_len(nrow(naming_df))) {
    lx <- stack_x + ((i - 1) %% 2) * 0.155
    ly <- legend_y - floor((i - 1) / 2) * 0.024
    grid.points(unit(lx, "npc"), unit(ly, "npc"), pch = 15, size = unit(0.10, "char"),
                gp = gpar(col = naming_df$color[[i]]))
    grid.text(paste0(naming_df$source[[i]], ": ", fmt_int(naming_df$count[[i]])),
              x = lx + 0.010, y = ly, just = c("left", "center"),
              gp = gpar(fontsize = 7.0, col = line_col))
  }

  note_x <- x0 + width * 0.225
  grid.roundrect(x = note_x, y = y0 + 0.057, width = width * 0.38, height = 0.054,
                 r = unit(0.008, "snpc"),
                 gp = gpar(fill = note_fill, col = "#D9B36A", lwd = 0.8))
  grid.text("BUSCO is a supplemental marker set;\nannotation support is primarily evidence-layer based.",
            x = note_x, y = y0 + 0.057,
            gp = gpar(fontsize = 7.2, col = line_col, lineheight = 0.95))

  bar_left <- x0 + width * 0.610
  bar_right <- x0 + width - 0.065
  bar_w <- bar_right - bar_left
  ys <- seq(y0 + height - 0.078, y0 + 0.075, length.out = nrow(rescue_df))
  for (i in seq_len(nrow(rescue_df))) {
    grid.text(rescue_df$label[[i]], x = bar_left - 0.012, y = ys[[i]],
              just = c("right", "center"),
              gp = gpar(fontsize = 7.1, col = line_col))
    grid.rect(x = bar_left, y = ys[[i]], width = bar_w, height = 0.022,
              just = c("left", "center"),
              gp = gpar(fill = "#F1F1F1", col = "#BFBFBF", lwd = 0.5))
    grid.rect(x = bar_left, y = ys[[i]],
              width = bar_w * rescue_df$count[[i]] / rescue_df$total[[i]],
              height = 0.022, just = c("left", "center"),
              gp = gpar(fill = rescue_df$color[[i]], col = NA))
    grid.text(paste0(fmt_int(rescue_df$count[[i]]), " (",
                     fmt_pct(rescue_df$count[[i]], rescue_df$total[[i]]), ")"),
              x = min(bar_left + bar_w * rescue_df$count[[i]] / rescue_df$total[[i]] + 0.006,
                      bar_right - 0.002),
              y = ys[[i]],
              just = c("left", "center"),
              gp = gpar(fontsize = 7.0, col = line_col))
  }
}

galba_genes <- stream_count_gff_feature(paths$galba_gff, "gene")
galba_mrna <- stream_count_gff_feature(paths$galba_gff, "mRNA")
protein_count <- stream_count_fasta(paths$proteins)
cds_count <- stream_count_fasta(paths$cds)
all_ids <- readLines(paths$all_ids, warn = FALSE)
eggnog_ids <- readLines(paths$eggnog_ids, warn = FALSE)
unannotated_ids <- readLines(paths$unannotated_ids, warn = FALSE)
tv_all <- read_blast_unique_queries(paths$tv_blast_all, normalize_dot = TRUE)
tv_rescue <- read_blast_unique_queries(paths$tv_rescue, normalize_dot = FALSE)
tv_mapping <- stream_tv_mapping(paths$tv_mapping_gff)
tv_reference_total <- stream_count_fasta(paths$tv_reference_proteins)
busco <- rbind(
  read_busco_summary(paths$busco_miniprot, "Miniprot"),
  read_busco_summary(paths$busco_metaeuk, "MetaEuk")
)
interpro_ids <- stream_unique_first_col(paths$interpro_tsv)

final_table <- read.delim(paths$final_table, stringsAsFactors = FALSE, check.names = FALSE)
final_rows <- nrow(final_table)
go_count <- sum(present(final_table$GO_Terms))
final_interpro_count <- sum(present(final_table$InterPro_IDs))
kegg_count <- sum(present(final_table$KEGG_Pathways))
naming_counts <- table(final_table$Naming_Source)
naming_df <- data.frame(
  source = c("Tv_Rescue", "EggNOG", "InterProScan", "None"),
  short = c("Tv", "EggNOG", "IPS", "None"),
  count = as.integer(naming_counts[c("Tv_Rescue", "EggNOG", "InterProScan", "None")]),
  color = c(tv_col, eggnog_col, ips_col, none_col),
  stringsAsFactors = FALSE
)

unannotated_with_ips <- length(intersect(unannotated_ids, interpro_ids))
unannotated_with_tv <- length(intersect(unannotated_ids, tv_rescue$ids))
unannotated_with_either <- length(union(intersect(unannotated_ids, interpro_ids),
                                       intersect(unannotated_ids, tv_rescue$ids)))
unannotated_without_either <- length(unannotated_ids) - unannotated_with_either

metrics <- list(
  galba_genes = galba_genes,
  galba_mrna = galba_mrna,
  proteins = protein_count,
  cds = cds_count,
  all_ids = length(all_ids),
  eggnog_annotated = length(eggnog_ids),
  unannotated = length(unannotated_ids),
  tv_all_rows = tv_all$rows,
  tv_all_unique = tv_all$unique_queries,
  tv_rescue_rows = tv_rescue$rows,
  tv_rescue_unique = tv_rescue$unique_queries,
  interproscan_unique = length(interpro_ids),
  final_rows = final_rows,
  go_count = go_count,
  final_interpro_count = final_interpro_count,
  kegg_count = kegg_count,
  unannotated_with_ips = unannotated_with_ips,
  unannotated_with_tv = unannotated_with_tv,
  unannotated_with_either = unannotated_with_either,
  unannotated_without_either = unannotated_without_either,
  tv_reference_total = tv_reference_total
)

evidence_counts <- data.frame(
  label = c(
    "GALBA gene features",
    "Extracted proteins",
    "InterProScan-supported proteins",
    "EggNOG annotated",
    "GO terms in final table",
    "InterPro IDs in final table",
    "KEGG pathways in final table"
  ),
  count = c(
    galba_genes,
    protein_count,
    length(interpro_ids),
    length(eggnog_ids),
    go_count,
    final_interpro_count,
    kegg_count
  ),
  color = c(struct_col, struct_col, ips_col, eggnog_col, "#6F9FBF", ips_col, "#8CA35C"),
  stringsAsFactors = FALSE
)

rescue_df <- data.frame(
  label = c("Tv rescue hit", "InterProScan record", "Either rescue layer", "No rescue support"),
  count = c(unannotated_with_tv, unannotated_with_ips, unannotated_with_either, unannotated_without_either),
  total = rep(length(unannotated_ids), 4),
  color = c(tv_col, ips_col, struct_col, none_col),
  stringsAsFactors = FALSE
)

metrics_out <- data.frame(
  metric = c(
    "galba_gene_features",
    "galba_mrna_features",
    "predicted_protein_sequences",
    "predicted_cds_sequences",
    "eggnog_annotated_proteins",
    "eggnog_unannotated_proteins",
    "full_tv_blast_rows",
    "full_tv_blast_unique_ttenax_queries",
    "tv_rescue_rows_for_eggnog_unannotated",
    "tv_rescue_unique_queries_for_eggnog_unannotated",
    "interproscan_unique_proteins",
    "final_master_table_rows",
    "final_table_go_terms",
    "final_table_interpro_ids",
    "final_table_kegg_pathways",
    "unannotated_with_tv_rescue",
    "unannotated_with_interproscan",
    "unannotated_with_tv_or_interproscan",
    "unannotated_without_tv_or_interproscan",
    "tv_mapping_mrna_alignment_records",
    "tv_mapping_unique_tv_targets",
    "tv_mapping_reference_tv_proteins",
    "tv_mapping_high_identity_records_identity_ge_0p6",
    "tv_mapping_high_identity_unique_tv_targets"
  ),
  value = c(
    galba_genes,
    galba_mrna,
    protein_count,
    cds_count,
    length(eggnog_ids),
    length(unannotated_ids),
    tv_all$rows,
    tv_all$unique_queries,
    tv_rescue$rows,
    tv_rescue$unique_queries,
    length(interpro_ids),
    final_rows,
    go_count,
    final_interpro_count,
    kegg_count,
    unannotated_with_tv,
    unannotated_with_ips,
    unannotated_with_either,
    unannotated_without_either,
    tv_mapping$mrna_records,
    tv_mapping$unique_targets,
    tv_reference_total,
    tv_mapping$high_identity_records,
    tv_mapping$high_identity_unique_targets
  ),
  source = c(
    paths$galba_gff,
    paths$galba_gff,
    paths$proteins,
    paths$cds,
    paths$eggnog_ids,
    paths$unannotated_ids,
    paths$tv_blast_all,
    paths$tv_blast_all,
    paths$tv_rescue,
    paths$tv_rescue,
    paths$interpro_tsv,
    paths$final_table,
    paths$final_table,
    paths$final_table,
    paths$final_table,
    paths$tv_rescue,
    paths$interpro_tsv,
    paste(paths$tv_rescue, paths$interpro_tsv, sep = "; "),
    paste(paths$tv_rescue, paths$interpro_tsv, sep = "; "),
    paths$tv_mapping_gff,
    paths$tv_mapping_gff,
    paths$tv_reference_proteins,
    paths$tv_mapping_gff,
    paths$tv_mapping_gff
  ),
  stringsAsFactors = FALSE
)
write.table(metrics_out, paste0(output_prefix, "_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(evidence_counts[, c("label", "count")],
            paste0(output_prefix, "_evidence_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(busco, paste0(output_prefix, "_busco_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(naming_df[, c("source", "count")],
            paste0(output_prefix, "_naming_sources.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rescue_df[, c("label", "count", "total")],
            paste0(output_prefix, "_eggnog_unannotated_rescue.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  grid.text("Figure 4. Annotation workflow and summary of annotation completeness",
            x = 0.5, y = 0.966,
            gp = gpar(fontsize = 22, fontface = "bold", col = line_col))
  grid.text("Data-derived structural prediction, functional evidence layers, BUSCO context, and lineage-aware coding-space support",
            x = 0.5, y = 0.936,
            gp = gpar(fontsize = 10.6, col = muted_col))

  draw_panel_title("A", "Annotation pipeline", 0.045, 0.888)
  draw_panel_title("B", "Gene, protein, and evidence counts", 0.535, 0.888)
  draw_panel_title("C", "Completeness assessment", 0.045, 0.488)
  draw_panel_title("D", "Final evidence integration", 0.535, 0.488)

  draw_pipeline_panel(0.055, 0.555, 0.400, 0.300, metrics)
  draw_horizontal_bar_panel(0.545, 0.555, 0.400, 0.300, evidence_counts, protein_count)
  draw_completeness_panel(0.055, 0.155, 0.400, 0.300, busco, tv_mapping, tv_reference_total)
  draw_naming_rescue_panel(0.545, 0.155, 0.400, 0.300, naming_df, rescue_df, protein_count)

}

png(paste0(output_prefix, ".png"), width = 16, height = 10.5, units = "in", res = 300)
draw_figure()
dev.off()

pdf(paste0(output_prefix, ".pdf"), width = 16, height = 10.5, useDingbats = FALSE)
draw_figure()
dev.off()

message("Wrote: ", paste0(output_prefix, ".png"))
message("Wrote: ", paste0(output_prefix, ".pdf"))
message("Wrote: ", paste0(output_prefix, "_metrics.tsv"))
message("Wrote: ", paste0(output_prefix, "_evidence_counts.tsv"))
message("Wrote: ", paste0(output_prefix, "_busco_metrics.tsv"))
message("Wrote: ", paste0(output_prefix, "_naming_sources.tsv"))
message("Wrote: ", paste0(output_prefix, "_eggnog_unannotated_rescue.tsv"))
