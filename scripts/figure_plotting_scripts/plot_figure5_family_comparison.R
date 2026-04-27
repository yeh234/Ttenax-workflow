#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args) >= 1) {
  args[[1]]
} else {
  "figures/Tt_Figure5_family_comparison_v1"
}

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

library(grid)

line_col <- "#111111"
axis_col <- "#303030"
muted_col <- "#5C5C5C"
grid_col <- "#E6E6E6"
panel_fill <- "#FFFFFF"
note_fill <- "#FFF8EC"
tv_col <- "#2E8B7C"
ips_col <- "#7A5AA6"
lrr_col <- "#3A7CA5"
warn_col <- "#B65C4A"
empty_col <- "#D6D6D6"

paths <- list(
  full_tv_blast = "dorado_polish_gpu1/Tt_vs_Tv_annotation.txt",
  tv_rescue = "dorado_polish_gpu1/eggnogHit/unannotated_vs_Tv_with_desc.txt",
  interpro_tsv = "dorado_polish_gpu1/interproscan_results/Tt_predicted_proteins.aa.fasta.tsv"
)

missing_paths <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing_paths) > 0) {
  stop("Missing required input files: ", paste(missing_paths, collapse = ", "))
}

families <- c("BspA", "LRR non-BspA", "Armadillo", "Kinase")

fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_pct <- function(n, d) sprintf("%.1f%%", n / d * 100)

wrap_text <- function(x, width = 42) {
  paste(strwrap(x, width = width), collapse = "\n")
}

new_set <- function() new.env(hash = TRUE, parent = emptyenv())
add_to_set <- function(env, id) {
  if (!is.na(id) && nzchar(id)) assign(id, TRUE, envir = env)
}
set_keys <- function(env) ls(env, all.names = TRUE)

family_hits <- function(desc) {
  list(
    BspA = grepl("bspa", desc, fixed = TRUE),
    `LRR non-BspA` = grepl("leucine-rich repeat", desc, fixed = TRUE) &
      !grepl("bspa", desc, fixed = TRUE),
    Armadillo = grepl("armadillo", desc, fixed = TRUE) |
      grepl("arm repeat", desc, fixed = TRUE),
    Kinase = grepl("kinase", desc, fixed = TRUE)
  )
}

read_blast_family_counts <- function(path, source_label) {
  tab <- read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  if (ncol(tab) < 13) {
    stop("Expected at least 13 BLAST columns in ", path)
  }
  query_id <- gsub("\\.", "", tab[[1]])
  desc <- tolower(tab[[13]])
  hit_list <- family_hits(desc)

  do.call(rbind, lapply(families, function(family) {
    idx <- hit_list[[family]]
    data.frame(
      source = source_label,
      family = family,
      hit_rows = sum(idx),
      unique_queries = length(unique(query_id[idx])),
      stringsAsFactors = FALSE
    )
  }))
}

read_interpro_families <- function(path) {
  armadillo <- new_set()
  bspa <- new_set()
  lrr <- new_set()
  kinase <- new_set()
  signalp <- new_set()
  tm_or_phobius <- new_set()
  all_ids <- new_set()

  con <- file(path, "rt")
  on.exit(close(con), add = TRUE)

  repeat {
    lines <- readLines(con, n = 1000, warn = FALSE)
    if (length(lines) == 0) break
    fields <- strsplit(lines, "\t", fixed = TRUE)

    for (row in fields) {
      if (length(row) < 6) next
      protein_id <- row[[1]]
      analysis <- row[[4]]
      signature_acc <- row[[5]]
      description <- tolower(row[[6]])

      add_to_set(all_ids, protein_id)

      if (grepl("SignalP", analysis, fixed = TRUE)) add_to_set(signalp, protein_id)
      if (grepl("TMHMM", analysis, fixed = TRUE) ||
          grepl("Phobius", analysis, fixed = TRUE)) {
        add_to_set(tm_or_phobius, protein_id)
      }

      if (grepl("armadillo", description, fixed = TRUE) ||
          grepl("PF00514", signature_acc, fixed = TRUE)) {
        add_to_set(armadillo, protein_id)
      }
      if (grepl("bspa", description, fixed = TRUE) ||
          grepl("PF13516", signature_acc, fixed = TRUE)) {
        add_to_set(bspa, protein_id)
      }
      if (grepl("leucine-rich repeat", description, fixed = TRUE) ||
          grepl("PF00560", signature_acc, fixed = TRUE)) {
        add_to_set(lrr, protein_id)
      }
      if (grepl("kinase", description, fixed = TRUE) &&
          grepl("protein", description, fixed = TRUE)) {
        add_to_set(kinase, protein_id)
      }
    }
  }

  arm_ids <- set_keys(armadillo)
  bspa_ids <- set_keys(bspa)
  lrr_ids <- set_keys(lrr)
  kinase_ids <- set_keys(kinase)
  other_lrr_ids <- setdiff(lrr_ids, bspa_ids)
  signalp_ids <- set_keys(signalp)
  tm_or_phobius_ids <- set_keys(tm_or_phobius)

  list(
    unique_interpro_ids = length(set_keys(all_ids)),
    domain_counts = data.frame(
      family = families,
      interproscan_domain_proteins = c(
        length(bspa_ids),
        length(other_lrr_ids),
        length(arm_ids),
        length(kinase_ids)
      ),
      stringsAsFactors = FALSE
    ),
    signal = data.frame(
      family = c("Armadillo", "BspA"),
      signalp_positive = c(
        length(intersect(arm_ids, signalp_ids)),
        length(intersect(bspa_ids, signalp_ids))
      ),
      total = c(length(arm_ids), length(bspa_ids)),
      tmhmm_or_phobius_positive = c(
        length(intersect(arm_ids, tm_or_phobius_ids)),
        length(intersect(bspa_ids, tm_or_phobius_ids))
      ),
      stringsAsFactors = FALSE
    )
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
            gp = gpar(fill = panel_fill, col = line_col, lwd = 1.0))
}

draw_axes <- function(left, bottom, right, top) {
  grid.lines(unit(c(left, right), "npc"), unit(c(bottom, bottom), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))
  grid.lines(unit(c(left, left), "npc"), unit(c(bottom, top), "npc"),
             gp = gpar(col = axis_col, lwd = 0.9))
}

draw_arrow <- function(x0, y0, x1, y1, lwd = 1.4) {
  grid.lines(unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
             gp = gpar(col = line_col, lwd = lwd),
             arrow = arrow(length = unit(0.10, "inches"), type = "closed"))
}

draw_label_box <- function(x, y, width, height, title, body, fill = "#FBFBFB",
                           title_col = line_col) {
  grid.roundrect(x = x, y = y, width = width, height = height,
                 r = unit(0.010, "snpc"),
                 gp = gpar(fill = fill, col = line_col, lwd = 1.0))
  grid.text(title, x = x, y = y + height / 2 - 0.022,
            gp = gpar(fontsize = 10.6, fontface = "bold", col = title_col))
  grid.text(body, x = x, y = y - 0.014,
            gp = gpar(fontsize = 8.0, col = line_col, lineheight = 0.94))
}

draw_family_bar_panel <- function(x0, y0, width, height, counts) {
  draw_panel_frame(x0, y0, width, height)

  left <- x0 + 0.068
  right <- x0 + width - 0.030
  bottom <- y0 + 0.075
  top <- y0 + height - 0.060
  plot_w <- right - left
  plot_h <- top - bottom
  y_max <- ceiling(max(counts$interproscan_domain_proteins,
                       counts$tv_blast_text_unique_queries) * 1.20 / 200) * 200
  y_max <- max(y_max, 1600)
  map_y <- function(v) bottom + plot_h * v / y_max

  for (tick in seq(0, y_max, by = 400)) {
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

  group_centers <- left + plot_w * (seq_along(families) - 0.5) / length(families)
  bar_w <- plot_w / length(families) * 0.25

  for (i in seq_along(families)) {
    tv_value <- counts$tv_blast_text_unique_queries[[i]]
    ips_value <- counts$interproscan_domain_proteins[[i]]
    x_tv <- group_centers[[i]] - bar_w * 0.58
    x_ips <- group_centers[[i]] + bar_w * 0.58

    grid.rect(x = x_tv, y = bottom, width = bar_w, height = map_y(tv_value) - bottom,
              just = c("center", "bottom"),
              gp = gpar(fill = tv_col, col = NA))
    grid.rect(x = x_ips, y = bottom, width = bar_w, height = map_y(ips_value) - bottom,
              just = c("center", "bottom"),
              gp = gpar(fill = ips_col, col = NA))

    grid.text(fmt_int(tv_value), x = x_tv, y = map_y(tv_value) + 0.014,
              gp = gpar(fontsize = 8.0, col = tv_col, fontface = "bold"))
    grid.text(fmt_int(ips_value), x = x_ips, y = map_y(ips_value) + 0.014,
              gp = gpar(fontsize = 8.0, col = ips_col, fontface = "bold"))
    grid.text(gsub(" ", "\n", families[[i]], fixed = TRUE),
              x = group_centers[[i]], y = bottom - 0.038,
              gp = gpar(fontsize = 8.2, col = line_col, lineheight = 0.90))
  }

  grid.text("Protein count", x = left - 0.047, y = bottom + plot_h / 2,
            rot = 90, gp = gpar(fontsize = 8.8, col = line_col))

  legend_x <- left + 0.010
  legend_y <- top + 0.025
  grid.rect(x = legend_x, y = legend_y, width = 0.014, height = 0.014,
            gp = gpar(fill = tv_col, col = NA))
  grid.text("Tv BLAST text labels", x = legend_x + 0.020, y = legend_y,
            just = c("left", "center"), gp = gpar(fontsize = 8.8, col = line_col))
  grid.rect(x = legend_x + 0.170, y = legend_y, width = 0.014, height = 0.014,
            gp = gpar(fill = ips_col, col = NA))
  grid.text("InterProScan domain calls", x = legend_x + 0.190, y = legend_y,
            just = c("left", "center"), gp = gpar(fontsize = 8.8, col = line_col))

}

draw_reinterpretation_panel <- function(x0, y0, width, height, counts) {
  draw_panel_frame(x0, y0, width, height)

  tv_lrr <- sum(counts$tv_blast_text_unique_queries[
    counts$family %in% c("BspA", "LRR non-BspA")
  ])
  ips_lrr <- sum(counts$interproscan_domain_proteins[
    counts$family %in% c("BspA", "LRR non-BspA")
  ])
  tv_arm <- counts$tv_blast_text_unique_queries[counts$family == "Armadillo"]
  ips_arm <- counts$interproscan_domain_proteins[counts$family == "Armadillo"]
  tv_kin <- counts$tv_blast_text_unique_queries[counts$family == "Kinase"]
  ips_kin <- counts$interproscan_domain_proteins[counts$family == "Kinase"]

  box_w <- width - 0.060
  box_h <- 0.078
  cx <- x0 + width / 2
  y_top <- y0 + height - 0.065
  y_bot <- y0 + 0.120

  draw_label_box(
    cx, y_top, box_w, box_h,
    "Text-label view",
    paste0(
      "Tv BLAST descriptions\n",
      "BspA + non-BspA LRR: ", fmt_int(tv_lrr), "\n",
      "Armadillo: ", fmt_int(tv_arm), "; Kinase: ", fmt_int(tv_kin)
    ),
    fill = "#F3FAF6",
    title_col = tv_col
  )
  draw_arrow(cx, y_top - box_h / 2 - 0.006, cx, y_bot + box_h / 2 + 0.006)
  grid.text("reinterpretation", x = cx + 0.082, y = (y_top + y_bot) / 2,
            gp = gpar(fontsize = 8.0, col = muted_col))
  draw_label_box(
    cx, y_bot, box_w, box_h,
    "Domain-architecture view",
    paste0(
      "InterProScan structural calls\n",
      "BspA + non-BspA LRR: ", fmt_int(ips_lrr), "\n",
      "Armadillo: ", fmt_int(ips_arm), "; Kinase: ", fmt_int(ips_kin)
    ),
    fill = "#F7F2FA",
    title_col = ips_col
  )

  grid.text(
    wrap_text(
      "The domain layer reveals a much broader LRR/BspA-related repertoire than text-label transfer alone. Armadillo behaves in the opposite direction: many Tv text labels but few InterPro ARM-domain calls.",
      48
    ),
    x = cx, y = y0 + 0.033,
    gp = gpar(fontsize = 7.4, col = muted_col, lineheight = 0.94)
  )
}

draw_signal_panel <- function(x0, y0, width, height, signal_df) {
  draw_panel_frame(x0, y0, width, height)

  left <- x0 + 0.080
  right <- x0 + width * 0.480
  bottom <- y0 + 0.090
  top <- y0 + height - 0.075
  plot_w <- right - left
  y_positions <- c(top - 0.045, bottom + 0.040)
  names(y_positions) <- signal_df$family

  grid.text("Canonical SignalP calls among domain-supported proteins",
            x = left, y = y0 + height - 0.038,
            just = c("left", "center"),
            gp = gpar(fontsize = 11.0, fontface = "bold", col = line_col))

  for (tick in seq(0, 100, by = 25)) {
    tx <- left + plot_w * tick / 100
    grid.lines(unit(c(tx, tx), "npc"), unit(c(bottom, top), "npc"),
               gp = gpar(col = grid_col, lwd = 0.6))
    grid.text(paste0(tick, "%"), x = tx, y = bottom - 0.035,
              gp = gpar(fontsize = 7.8, col = line_col))
  }
  grid.lines(unit(c(left, right), "npc"), unit(c(bottom, bottom), "npc"),
             gp = gpar(col = axis_col, lwd = 0.8))

  for (i in seq_len(nrow(signal_df))) {
    family <- signal_df$family[[i]]
    pct <- signal_df$signalp_positive[[i]] / signal_df$total[[i]] * 100
    y <- y_positions[[family]]
    grid.text(family, x = left - 0.018, y = y,
              just = c("right", "center"),
              gp = gpar(fontsize = 9.0, fontface = "bold", col = line_col))
    grid.roundrect(x = left, y = y, width = plot_w, height = 0.026,
                   r = unit(0.010, "snpc"),
                   just = c("left", "center"),
                   gp = gpar(fill = empty_col, col = NA))
    if (pct > 0) {
      grid.roundrect(x = left, y = y, width = plot_w * pct / 100, height = 0.026,
                     r = unit(0.010, "snpc"),
                     just = c("left", "center"),
                     gp = gpar(fill = lrr_col, col = NA))
    }
    grid.points(x = left + plot_w * pct / 100, y = y, pch = 21, size = unit(0.050, "in"),
                gp = gpar(fill = warn_col, col = line_col, lwd = 0.6))
    grid.text(
      paste0(
        fmt_int(signal_df$signalp_positive[[i]]), " / ",
        fmt_int(signal_df$total[[i]]), " (",
        fmt_pct(signal_df$signalp_positive[[i]], signal_df$total[[i]]), ")"
      ),
      x = right - 0.012, y = y,
      just = c("right", "center"),
      gp = gpar(fontsize = 8.8, col = line_col)
    )
  }

  caution_x <- x0 + width * 0.735
  caution_y <- y0 + height / 2
  caution_w <- width * 0.405
  caution_h <- height - 0.090
  grid.roundrect(x = caution_x, y = caution_y, width = caution_w, height = caution_h,
                 r = unit(0.010, "snpc"),
                 gp = gpar(fill = note_fill, col = line_col, lwd = 1.0))
  grid.text("Interpretation guardrail",
            x = caution_x, y = caution_y + caution_h / 2 - 0.035,
            gp = gpar(fontsize = 11.0, fontface = "bold", col = line_col))
  grid.text(
    wrap_text(
      "Armadillo proteins lack canonical signal peptides in this run, supporting a mostly non-secretory interpretation. For BspA/LRR proteins, absence of SignalP should remain cautious because incomplete N-termini, divergent signal peptides, or non-canonical trafficking can reduce prediction sensitivity.",
      58
    ),
    x = caution_x, y = caution_y - 0.008,
    gp = gpar(fontsize = 8.4, col = line_col, lineheight = 0.96)
  )
}

tv_full_counts <- read_blast_family_counts(paths$full_tv_blast, "Tv BLAST all hits")
tv_rescue_counts <- read_blast_family_counts(paths$tv_rescue, "Tv rescue subset")
interpro <- read_interpro_families(paths$interpro_tsv)

family_counts <- merge(
  tv_full_counts[, c("family", "hit_rows", "unique_queries")],
  tv_rescue_counts[, c("family", "hit_rows", "unique_queries")],
  by = "family",
  suffixes = c("_full_tv_blast", "_tv_rescue"),
  sort = FALSE
)
family_counts <- merge(family_counts, interpro$domain_counts, by = "family", sort = FALSE)
family_counts <- family_counts[match(families, family_counts$family), ]
names(family_counts) <- c(
  "family",
  "tv_blast_text_hit_rows",
  "tv_blast_text_unique_queries",
  "tv_rescue_text_hit_rows",
  "tv_rescue_text_unique_queries",
  "interproscan_domain_proteins"
)

signal_df <- interpro$signal
signal_df$signalp_percent <- signal_df$signalp_positive / signal_df$total * 100

blast_table <- read.delim(paths$full_tv_blast, header = FALSE, stringsAsFactors = FALSE,
                          quote = "", comment.char = "")
rescue_table <- read.delim(paths$tv_rescue, header = FALSE, stringsAsFactors = FALSE,
                           quote = "", comment.char = "")

metrics_out <- data.frame(
  metric = c(
    "full_tv_blast_rows",
    "full_tv_blast_unique_ttenax_queries",
    "tv_rescue_rows_for_eggnog_unannotated",
    "tv_rescue_unique_ttenax_queries",
    "interproscan_unique_proteins",
    "tv_blast_text_bspa_plus_lrr_unique_queries",
    "interproscan_bspa_plus_lrr_domain_proteins",
    "tv_blast_text_armadillo_unique_queries",
    "interproscan_armadillo_domain_proteins",
    "tv_blast_text_kinase_unique_queries",
    "interproscan_kinase_domain_proteins"
  ),
  value = c(
    nrow(blast_table),
    length(unique(gsub("\\.", "", blast_table[[1]]))),
    nrow(rescue_table),
    length(unique(gsub("\\.", "", rescue_table[[1]]))),
    interpro$unique_interpro_ids,
    sum(family_counts$tv_blast_text_unique_queries[
      family_counts$family %in% c("BspA", "LRR non-BspA")
    ]),
    sum(family_counts$interproscan_domain_proteins[
      family_counts$family %in% c("BspA", "LRR non-BspA")
    ]),
    family_counts$tv_blast_text_unique_queries[family_counts$family == "Armadillo"],
    family_counts$interproscan_domain_proteins[family_counts$family == "Armadillo"],
    family_counts$tv_blast_text_unique_queries[family_counts$family == "Kinase"],
    family_counts$interproscan_domain_proteins[family_counts$family == "Kinase"]
  ),
  source = c(
    paths$full_tv_blast,
    paths$full_tv_blast,
    paths$tv_rescue,
    paths$tv_rescue,
    paths$interpro_tsv,
    paths$full_tv_blast,
    paths$interpro_tsv,
    paths$full_tv_blast,
    paths$interpro_tsv,
    paths$full_tv_blast,
    paths$interpro_tsv
  ),
  stringsAsFactors = FALSE
)

write.table(family_counts, paste0(output_prefix, "_family_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(signal_df, paste0(output_prefix, "_signal_peptide.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metrics_out, paste0(output_prefix, "_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  grid.text("Figure 5. Lineage-aware versus domain-based annotation of expanded gene families",
            x = 0.5, y = 0.966,
            gp = gpar(fontsize = 22, fontface = "bold", col = line_col))
  grid.text("Data-derived comparison of Tv BLAST text labels, InterProScan domain calls, and SignalP localization context",
            x = 0.5, y = 0.936,
            gp = gpar(fontsize = 10.6, col = muted_col))

  draw_panel_title("A", "Family comparison", 0.045, 0.888)
  draw_panel_title("B", "Conceptual reinterpretation", 0.648, 0.888)
  draw_panel_title("C", "Signal peptide context", 0.045, 0.488)

  draw_family_bar_panel(0.055, 0.545, 0.555, 0.305, family_counts)
  draw_reinterpretation_panel(0.655, 0.545, 0.290, 0.305, family_counts)
  draw_signal_panel(0.055, 0.150, 0.890, 0.300, signal_df)

  grid.text(
    "Source data: Tt_vs_Tv_annotation.txt, EggNOG-unannotated Tv rescue BLAST table, and InterProScan TSV parsed with the same family rules as parse_ips.py.",
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
message("Wrote: ", paste0(output_prefix, "_family_counts.tsv"))
message("Wrote: ", paste0(output_prefix, "_signal_peptide.tsv"))
message("Wrote: ", paste0(output_prefix, "_metrics.tsv"))
