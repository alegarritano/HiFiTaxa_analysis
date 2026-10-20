#!/usr/bin/env Rscript
# =====================================================================
# Figure 1: the HiFiTaxa pipeline, drawn at its final width (180 mm) so the
# font sizes in the file are the printed sizes (titles 7.5 pt, body 6 pt).
# Pure ggplot2 (geom_shape from ggforce for rounded boxes); no data input.
# Writes fig1_pipeline.png (300 dpi) and fig1_pipeline.pdf to
# HIFITAXA_OUT/figures (default: <this folder>/regenerated/figures).
# =====================================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggforce); library(grid) })
.SCRIPT_DIR <- local({ a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd()) })
OUT <- Sys.getenv("HIFITAXA_OUT", unset = file.path(.SCRIPT_DIR, "regenerated"))
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)

# ---- palette (colour blind safe, muted) ----
C <- list(input = "#2b2d42", qc = "#4a4e69", panel16 = "#eef3fb", panel16b = "#9db4d6", panelits = "#fdf4e7", panelitsb = "#e0b98a",
          prep = "#ffffff", asv = "#f3eefa", asvb = "#8e6cc1", clf = "#6f4fb0", clfb = "#4d3380", prof = "#c94c6a", profb = "#a23552",
          out = "#1b9e77", outb = "#13704f", ref = "#f4f4f4", refb = "#8c8c8c", text = "#1a1a1a", arrow = "#555555")
W <- 180; H <- 114
boxes <- list(); texts <- list(); segs <- list()
box <- function(x1, y1, x2, y2, fill, col, lwd = 0.35, lty = 1) boxes[[length(boxes) + 1]] <<- data.frame(x = c(x1, x2, x2, x1), y = c(y1, y1, y2, y2), fill = fill, col = col, lwd = lwd, lty = lty, id = length(boxes) + 1)
txt <- function(x, y, label, size = 5, face = "plain", col = C$text, hjust = 0.5, vjust = 0.5, lh = 0.95) texts[[length(texts) + 1]] <<- data.frame(x = x, y = y, label = label, size = size, face = face, col = col, hjust = hjust, vjust = vjust, lh = lh)
seg <- function(x1, y1, x2, y2, arrow = TRUE, col = C$arrow, lty = 1, lwd = 0.4) segs[[length(segs) + 1]] <<- data.frame(x = x1, y = y1, xend = x2, yend = y2, arrow = arrow, col = col, lty = lty, lwd = lwd)
PT <- 25.4 / 72   # mm per point
# a box with a bold title line and a body of up to six 5 pt lines, laid out from the top edge
labelled_box <- function(x1, y1, x2, y2, title, body = NULL, fill, col, tcol = C$text, bcol = C$text, tsize = 7.5, bsize = 5, lwd = 0.35) {
  box(x1, y1, x2, y2, fill, col, lwd = lwd); xm <- (x1 + x2) / 2
  if (is.null(body)) { txt(xm, (y1 + y2) / 2, title, size = tsize, face = "bold", col = tcol); return(invisible()) }
  ttop <- y2 - 1.6; txt(xm, ttop, title, size = tsize, face = "bold", col = tcol, vjust = 1)
  txt(xm, ttop - tsize * PT - 0.8, body, size = bsize, col = bcol, vjust = 1, lh = 1.0) }

# ---- front end (left column) ----
labelled_box(4, 88, 30, 100, "PacBio HiFi reads", "demultiplexed FASTQ files\nand a sample sheet", fill = C$input, col = C$input, tcol = "white", bcol = "white")
labelled_box(4, 68, 30, 82, "Quality filter", "seqkit: mean per base error\nprobability at most 1%\n(average quality Q20)", fill = "white", col = C$qc)
dm <- data.frame(x = c(17, 27, 17, 7), y = c(60, 51, 42, 51))
txt(17, 52.8, "Marker?", size = 7.5, face = "bold"); txt(17, 49.0, "16S rRNA gene\nor ITS", size = 5, lh = 1.0)
seg(17, 88, 17, 82.5); seg(17, 68, 17, 60.5)
seg(27, 51, 33, 51, arrow = FALSE); seg(33, 51, 33, 91, arrow = FALSE); seg(33, 91, 40, 91)
seg(33, 51, 33, 39, arrow = FALSE); seg(33, 39, 40, 39)

# ---- reference lane behind the classifier column ----
box(100, 1, 134, 113, fill = C$ref, col = C$refb, lwd = 0.35, lty = 2)
txt(117, 112.3, "References built at launch", size = 6.5, face = "bold", vjust = 1)
txt(117, 109.3, "one build per source, shared by all classifiers", size = 5, vjust = 1, col = "#444444")
txt(117, 1.8, "GTDB: SSU sequences of the current release\nUNITE: the supplied release file", size = 5, col = "#444444", vjust = 0, lh = 1.0)

# ---- 16S panel ----
box(36, 58, 148, 105, fill = C$panel16, col = C$panel16b, lwd = 0.4)
txt(38.5, 103.6, "16S rRNA gene (bacteria and archaea): GTDB SSU, r232 here", size = 7, face = "bold", hjust = 0, vjust = 1, col = "#2c4a7a")
labelled_box(40, 84, 66, 97, "Primer trimming", "cutadapt, linked adapter;\n27F and 1492R by default,\nboth primers required", fill = C$prep, col = C$qc)
labelled_box(70, 79, 98, 99, "ASV inference", "DADA2 denoise-ccs (QIIME 2)\n1,000 to 1,600 bases\nat most two expected errors\nchimera removal\nASV filter: 5 reads in\n1 sample by default", fill = C$asv, col = C$asvb)
labelled_box(102, 88, 132, 100, "BLCA", "Bayesian LCA on BLAST hits,\nper rank bootstrap\n16S rRNA gene default", fill = C$clf, col = C$clfb, tcol = "white", bcol = "white")
labelled_box(102, 73, 132, 85, "Naive Bayes", "genus by assignTaxonomy,\nspecies by exact match", fill = C$clf, col = C$clfb, tcol = "white", bcol = "white")
labelled_box(70, 61, 98, 71, "Length filter", "1,000 to 1,600 bases;\nuntrimmed reads", fill = C$prep, col = C$qc)
labelled_box(102, 60, 132, 72, "Emu", "read level profiler\n(expectation maximisation)", fill = C$prof, col = C$profb, tcol = "white", bcol = "white")
txt(53, 82.2, "ASV classifiers", size = 5, col = "#555555"); txt(53, 72.6, "read level profiler", size = 5, col = "#555555")
seg(37, 91, 37, 66, arrow = FALSE); seg(37, 66, 70, 66)
seg(66, 91, 70, 91); seg(98, 94, 102, 94); seg(98, 79, 102, 79); seg(98, 66, 102, 66)

# ---- ITS panel ----
box(36, 6, 148, 53, fill = C$panelits, col = C$panelitsb, lwd = 0.4)
txt(38.5, 51.6, "Fungal ITS: UNITE general release for Fungi, 19.02.2025 here", size = 7, face = "bold", hjust = 0, vjust = 1, col = "#7a4a12")
labelled_box(40, 32, 66, 46, "ITS extraction", "ITSxRust extracts the full\nITS region; no primer\ntrimming", fill = C$prep, col = C$qc)
labelled_box(70, 28, 98, 48, "ASV inference", "DADA2 denoise-ccs (QIIME 2)\n200 to 1,300 bases\nat most two expected errors\nchimera removal\nASV filter: 5 reads in\n1 sample by default", fill = C$asv, col = C$asvb)
labelled_box(102, 36, 132, 48, "BLCA", "Bayesian LCA on BLAST hits,\nper rank bootstrap", fill = C$clf, col = C$clfb, tcol = "white", bcol = "white")
labelled_box(102, 21, 132, 33, "Naive Bayes", "single step assignTaxonomy\nover seven ranks", fill = C$clf, col = C$clfb, tcol = "white", bcol = "white")
labelled_box(70, 9, 98, 19, "Length filter", "200 to 1,300 bases;\nextracted ITS reads", fill = C$prep, col = C$qc)
labelled_box(102, 8, 132, 20, "EMITS", "read level profiler\n(expectation maximisation)\nITS default", fill = C$prof, col = C$profb, tcol = "white", bcol = "white")
txt(53, 30.2, "ASV classifiers", size = 5, col = "#555555"); txt(53, 20.6, "read level profiler", size = 5, col = "#555555")
seg(66, 39, 67.5, 39, arrow = FALSE); seg(67.5, 39, 70, 39); seg(67.5, 39, 67.5, 14, arrow = FALSE); seg(67.5, 14, 70, 14)
seg(98, 42, 102, 42); seg(98, 27, 102, 27); seg(98, 14, 102, 14)

# ---- outputs ----
for (y in c(94, 79, 66, 42, 27, 14)) seg(132, y, 141, y, arrow = FALSE)
seg(141, 14, 141, 94, arrow = FALSE); seg(141, 56, 152, 56)
labelled_box(152, 40, 178, 72, "Outputs", "taxonomy table\nper classifier\n\nread tracking table\n\nASV by sample\nfeature table\n\nexecution report\nand parameters.txt", fill = "white", col = C$outb, tcol = C$outb, bsize = 5.5)

B <- do.call(rbind, boxes); Tt <- do.call(rbind, texts); Sg <- do.call(rbind, segs)
g <- ggplot() +
  geom_shape(data = B, aes(x, y, group = id), fill = B$fill, colour = NA, radius = unit(1.6, "mm")) +
  geom_shape(data = B, aes(x, y, group = id), fill = NA, colour = B$col, linewidth = B$lwd, linetype = B$lty, radius = unit(1.6, "mm")) +
  geom_polygon(data = dm, aes(x, y), fill = "white", colour = C$qc, linewidth = 0.35) +
  geom_segment(data = Sg[!Sg$arrow, ], aes(x = x, y = y, xend = xend, yend = yend), colour = Sg$col[!Sg$arrow], linewidth = Sg$lwd[!Sg$arrow], linetype = Sg$lty[!Sg$arrow], lineend = "round") +
  geom_segment(data = Sg[Sg$arrow, ], aes(x = x, y = y, xend = xend, yend = yend), colour = Sg$col[Sg$arrow], linewidth = Sg$lwd[Sg$arrow], linetype = Sg$lty[Sg$arrow], lineend = "round", linejoin = "mitre", arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_text(data = Tt, aes(x, y, label = label), size = Tt$size / .pt, fontface = Tt$face, colour = Tt$col, hjust = Tt$hjust, vjust = Tt$vjust, lineheight = Tt$lh) +
  coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE, clip = "off") + theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(1, 1, 1, 1, "mm"))
png_out <- file.path(OUT, "figures", "fig1_pipeline.png"); pdf_out <- file.path(OUT, "figures", "fig1_pipeline.pdf")
ggsave(png_out, g, width = W, height = H, units = "mm", dpi = 300, device = ragg::agg_png, bg = "white")
ggsave(pdf_out, g, width = W, height = H, units = "mm", device = pdf, useDingbats = FALSE, bg = "white")
cat("Wrote", png_out, "and", pdf_out, "\n")
