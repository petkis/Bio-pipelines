#!/usr/bin/env Rscript

# =========================
# AUTO-INSTALL DEPENDENCIES
# =========================
required_packages <- c("ggplot2", "dplyr", "stringr", "purrr", "readr", "patchwork")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos='http://cran.us.r-project.org')

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(patchwork) # Required for stitching independent plots
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript telo_histograms.R <base_dir>")
}

base_dir <- args[1]

# =========================
# VERSION MAP
# =========================

version_map <- c(
  "0.9.1" = "0.9.1 (2024 early)",
  "0.9.6" = "0.9.6 (2024 late)",
  "1.0.0" = "1.0.0 (2025-05)",
  "1.2.0" = "1.2.0 (2025-10)",
  "1.3.0" = "1.3.0 (2025-11)",
  "1.4.0" = "1.4.0 (2026-02)",
  "6.5.7" = "6.5.7 (2023)"
)

# =========================
# HELPERS
# =========================

model_family_from_model <- function(model) {
  case_when(
    str_detect(str_to_lower(model), "sup") ~ "SUP",
    str_detect(str_to_lower(model), "hac") ~ "HAC",
    TRUE ~ NA_character_
  )
}

# =========================
# LOAD AND PARSE
# =========================

files <- list.files(base_dir, pattern = "telomere_lengths.tsv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No files found.")

parse_file <- function(f) {
  df <- tryCatch({ read_tsv(f, show_col_types = FALSE) }, error = function(e) return(NULL))
  if (is.null(df) || !"length" %in% names(df)) return(NULL)

  df$length <- suppressWarnings(as.numeric(df$length))
  df <- df[!is.na(df$length) & df$length > 0, ]
  if (nrow(df) == 0) return(NULL)

  path <- str_split(f, "/")[[1]]
  folder <- path[which(str_detect(path, "dorado_|guppy_"))][1]
  if (is.na(folder)) return(NULL)

  parts <- str_split(folder, "_")[[1]]
  basecaller <- parts[1]
  version_raw <- str_remove(parts[2], "^v")
  version <- ifelse(version_raw %in% names(version_map), version_map[version_raw], version_raw)

  tool <- case_when(
    str_detect(f, "/ncrf/") ~ "NCRF",
    str_detect(f, "/telobp/") ~ "TeloBP",
    TRUE ~ "Unknown"
  )

  df$basecaller <- basecaller
  df$version <- version
  df$model_family <- model_family_from_model(paste(parts[3:length(parts)], collapse="_"))
  df$tool <- tool
  # New header format: Tool - Basecaller - Version
  df$display_label <- paste0(tool, " - ", str_to_title(basecaller), " - ", version)
  
  return(df)
}

df_all <- map_dfr(files, parse_file) %>%
  filter(basecaller %in% c("dorado", "guppy"), model_family %in% c("SUP", "HAC"))

# Ensure chronological order
ordered_levels <- unique(df_all$display_label)
df_all$display_label <- factor(df_all$display_label, levels = ordered_levels)

# =========================
# PLOTTING LOGIC
# =========================

create_individual_plot <- function(sub_data, label_name) {
  ggplot(sub_data, aes(x = length)) +
    geom_histogram(bins = 80, color = "black", fill = "steelblue") +
    theme_minimal(base_size = 11) +
    labs(
      title = label_name,
      x = "Telomere length (bp)",
      y = "Count"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      # Ensure axes are drawn on every plot
      axis.line = element_line(color = "black")
    )
}

# =========================
# GENERATE OUTPUTS
# =========================

out_root <- "histograms"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Group by Tool and Family for separate PDFs
groupings <- df_all %>% distinct(tool, model_family)

for (i in 1:nrow(groupings)) {
  curr_tool <- groupings$tool[i]
  curr_fam  <- groupings$model_family[i]
  
  plot_data <- df_all %>% filter(tool == curr_tool, model_family == curr_fam)
  if (nrow(plot_data) == 0) next
  
  # Create a list of separate plots for each version
  unique_labels <- levels(droplevels(plot_data$display_label))
  plot_list <- map(unique_labels, function(lbl) {
    create_individual_plot(filter(plot_data, display_label == lbl), lbl)
  })

  # Use patchwork to combine plots
  # ncol = 2 forces 2 columns
  # Every plot keeps its own axis labels because they are separate ggplot objects
  combined_plot <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = paste(curr_tool, "Analysis -", curr_fam, "Models"),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

  # Calculate size
  n_rows <- ceiling(length(plot_list) / 2)
  
  tool_dir <- file.path(out_root, str_to_lower(curr_tool))
  dir.create(tool_dir, showWarnings = FALSE)
  
  out_path <- file.path(tool_dir, paste0("hist_", curr_tool, "_", curr_fam, ".pdf"))
  
  ggsave(out_path, combined_plot, width = 12, height = max(5, n_rows * 4.5), units = "in", dpi = 300)
}

# =========================
# SUMMARY
# =========================

summary_stats <- df_all %>%
  group_by(tool, basecaller, version, model_family) %>%
  summarise(n = n(), median = median(length), .groups = "drop")

write_tsv(summary_stats, "histogram_summary_stats.tsv")
message("Files generated in 'histograms' folder. Each plot has full X/Y axes.")
