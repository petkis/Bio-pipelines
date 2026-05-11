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
    chromType = factor(chromType),
    chrom = factor(chrom),
    tool = factor(tool)
  ) %>%
  filter(!chrom %in% c("chrM", "chrMT"))

# ----------------------------
# LEGEND
# ----------------------------
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
  NoTelo = 0.6,
  Telo   = 0.5
)

chromtype_labels <- c(
  Masked = "Masked: telomeric sequence masked",
  NoTelo = "NoTelo: telomeres removed",
  Telo   = "Telo: telomeres included"
)

# ----------------------------
# THEME
# ----------------------------
base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

# ----------------------------
# CHROM ORDER
# ----------------------------
chr_order_key <- function(ch) {
  ch <- as.character(ch)
  ch <- sub("^chr", "", ch)
  if (grepl("^[0-9]+$", ch)) return(as.integer(ch))
  if (toupper(ch) == "X") return(23L)
  if (toupper(ch) == "Y") return(24L)
  return(999L)
}

# ----------------------------
# MAIN LOOP (PER TOOL)
# ----------------------------
for (tool_name in unique(mydat$tool)) {

  cat("Processing tool:", tool_name, "\n")

  tool_data <- mydat %>% filter(tool == tool_name)

  default_desc <- paste0(
    "Reads were simulated with wgsim (10,000 reads per length category).\n",
    tool_name, " was used for alignment.\n",
    "The reported percentage equals:\n",
    "( number of simulated reads aligned to the correct chromosome / 10,000 ) × 100."
  )

  if (nrow(tool_data) == 0) {
    cat("No data for tool:", tool_name, "\n")
    next
  }

  chrom_levels <- unique(tool_data$chrom)
  chrom_levels <- chrom_levels[order(vapply(chrom_levels, chr_order_key, integer(1)))]

  out_pdf <- paste0("All_chromosomes_", tool_name, ".pdf")
  pdf(out_pdf, width = 13, height = 7)

  for (chr in chrom_levels) {

    df_chr <- tool_data %>% filter(chrom == chr)

    # 🔥 CRITICAL FIX (prevents your crash)
    if (nrow(df_chr) == 0) {
      cat("Skipping empty:", chr, "\n")
      next
    }

    chr_char <- as.character(chr)
    chr_pretty <- toupper(sub("^chr", "", chr_char))
    title_txt <- paste0("chromosome", chr_pretty)

    p <- ggplot(df_chr, aes(x = Length_kb, y = Percent)) +
      geom_line(aes(color = chromType, group = chromType), alpha = 0.5) +
      geom_point(aes(color = chromType, shape = chromType, alpha = chromType), size = 2.5) +
      scale_alpha_manual(values = chromtype_alpha, guide = "none") +
      facet_wrap(~ REF, ncol = 3) +
      scale_y_continuous(limits = c(0, 100)) +
      scale_color_manual(
        name = "Reads setup",
        values = chromtype_cols,
        labels = chromtype_labels
      ) +
      scale_shape_manual(values = chromtype_shapes, guide = "none") +
      labs(
        title = paste0(title_txt, " (", tool_name, ")"),
        x = "Read length (kb)",
        y = "Percent"
      ) +
      base_theme

    # ----------------------------
    # GRID LAYOUT (PLOT + DESCRIPTION)
    # ----------------------------
    grid.newpage()

    pushViewport(
      viewport(
        layout = grid.layout(
          nrow = 2,
          heights = unit(c(0.84, 0.16), "npc")
        )
      )
    )

    # Plot
    print(p, vp = viewport(layout.pos.row = 1))

    # Description
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
