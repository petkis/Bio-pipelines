library(tidyverse)
library(grid)

# ----------------------------
# INPUT CHECK
# ----------------------------
if (!file.exists("combined_results.csv")) {
  stop("ERROR: combined_results.csv not found!")
}

dat <- read.csv("combined_results.csv")

if (nrow(dat) == 0) stop("ERROR: CSV is empty!")

required_cols <- c("tool", "chrom", "REF", "chromType", "Length", "Percent")
missing_cols <- setdiff(required_cols, colnames(dat))

if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

cat("CSV OK | Rows:", nrow(dat), "\n")

# ----------------------------
# PREP
# ----------------------------
mydat <- dat %>%
  mutate(
    Length_kb = Length / 1000,
    REF = factor(REF, levels = c("Masked", "NoTelo", "Telo")),
    chromType = factor(chromType, levels = c("Masked", "NoTelo", "Telo")),
    chrom = factor(chrom),
    tool = factor(tool)
  ) %>%
  filter(!chrom %in% c("chrM", "chrMT"))

# ----------------------------
# LABELS
# ----------------------------
ref_labels <- c(
  Masked = "Masked reference",
  NoTelo = "Reference without telomeres",
  Telo   = "Telomeric reference"
)

chromtype_cols <- c(
  Masked = "#F8766D",
  NoTelo = "#00BA38",
  Telo   = "#619CFF"
)

chromtype_shapes <- c(
  Masked = 16,
  NoTelo = 17,
  Telo   = 15
)

chromtype_alpha <- c(
  Masked = 1.0,
  NoTelo = 0.7,
  Telo   = 0.7
)

chromtype_labels <- c(
  Masked = "Masked: telomeric sequence masked in the reference",
  NoTelo = "NoTelo: subtelomeres/telomeres removed from the reference",
  Telo   = "Telo: complete reference incl. subtelomeres and telomeres"
)

# ----------------------------
# THEME
# ----------------------------
base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

readable_theme <- theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  )

# ----------------------------
# CHROM ORDER
# ----------------------------
chr_order_key <- function(ch) {
  ch <- as.character(ch)
  ch <- sub("^chr", "", ch)
  ch <- sub("_.*$", "", ch)

  if (grepl("^[0-9]+$", ch)) return(as.integer(ch))
  if (toupper(ch) == "X") return(23L)
  if (toupper(ch) == "Y") return(24L)
  return(999L)
}

pretty_chr_name <- function(ch) {
  ch <- as.character(ch)
  ch <- sub("^chr", "", ch)
  ch <- sub("_.*$", "", ch)
  paste0("chromosome", ch)
}

# ----------------------------
# FUNCTION: FULL PDF PER TOOL
# ----------------------------
make_tool_pdf <- function(tool_name) {

  cat("Processing tool:", tool_name, "\n")

  tool_data <- mydat %>% filter(tool == tool_name)

  if (nrow(tool_data) == 0) {
    cat("No data for tool:", tool_name, "\n")
    return(NULL)
  }

  chrom_levels <- unique(tool_data$chrom)
  chrom_levels <- chrom_levels[order(vapply(chrom_levels, chr_order_key, integer(1)))]

  out_pdf <- paste0("All_chromosomes_", tool_name, ".pdf")
  pdf(out_pdf, width = 13, height = 7)

  for (chr in chrom_levels) {

    df_chr <- tool_data %>% filter(chrom == chr)

    if (nrow(df_chr) == 0) {
      cat("Skipping empty:", chr, "\n")
      next
    }

    if (tool_name == "minimap2") {
      tool_label <- "Minimap2"
    } else if (tool_name == "winnowmap") {
      tool_label <- "Winnowmap"
    } else {
      tool_label <- tool_name
    }

    total_reads <- max(df_chr$N, na.rm = TRUE)
    if (!is.finite(total_reads)) {
      total_reads <- 10000
    }

    default_desc <- paste0(
      "Reads were simulated with wgsim.\n",
      tool_label, " was used for alignment against each reference setup.\n",
      "Accuracy percentage = reads aligned to the correct chromosome / total simulated reads × 100."
    )

    p <- ggplot(df_chr, aes(x = Length_kb, y = Percent)) +
      geom_line(aes(color = chromType, group = chromType), alpha = 0.5) +
      geom_point(aes(color = chromType, shape = chromType, alpha = chromType), size = 2.5) +
      scale_alpha_manual(values = chromtype_alpha, guide = "none") +
      facet_wrap(~ REF, ncol = 3, labeller = labeller(REF = ref_labels)) +
      scale_y_continuous(
        limits = c(0, 100),
        breaks = seq(0, 100, by = 25),
        labels = function(x) paste0(x, "%")
      ) +
      scale_color_manual(
        name = "Reads setup",
        values = chromtype_cols,
        labels = chromtype_labels
      ) +
      scale_shape_manual(values = chromtype_shapes, guide = "none") +
      labs(
        title = paste0(pretty_chr_name(chr), " (", tool_name, ")"),
        x = "Read length (kb)",
        y = "Accuracy percentage (%)"
      ) +
      base_theme

    grid.newpage()

    pushViewport(
      viewport(
        layout = grid.layout(
          nrow = 2,
          heights = unit(c(0.84, 0.16), "npc")
        )
      )
    )

    print(p, vp = viewport(layout.pos.row = 1))

    grid.text(
      default_desc,
      vp = viewport(layout.pos.row = 2),
      x = 0.5, y = 0.5,
      just = c("center", "center"),
      gp = gpar(fontsize = 9)
    )
  }

  dev.off()
  cat("Saved:", out_pdf, "\n")
}

# ----------------------------
# CREATE TWO PDFS
# ----------------------------
make_tool_pdf("minimap2")
make_tool_pdf("winnowmap")
