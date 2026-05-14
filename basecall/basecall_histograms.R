#!/usr/bin/env Rscript

# Create histogram PDFs for telomere length distributions.
# The script searches for telomere_lengths.tsv files, parses basecaller/version/model/tool
# information from folder names, creates grouped histogram PDFs, and saves summary statistics.

# Install required packages if they are missing.
required_packages <- c("ggplot2", "dplyr", "stringr", "purrr", "readr", "patchwork")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos='http://cran.us.r-project.org')

# Load packages without printing startup messages.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(patchwork)
})

# Read input directory from the command line.
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript telo_histograms.R <base_dir>")
}

base_dir <- args[1]

# Map raw basecaller versions to readable labels used in plot titles.
version_map <- c(
  "0.9.1" = "0.9.1 (2024 early)",
  "0.9.6" = "0.9.6 (2024 late)",
  "1.0.0" = "1.0.0 (2025-05)",
  "1.2.0" = "1.2.0 (2025-10)",
  "1.3.0" = "1.3.0 (2025-11)",
  "1.4.0" = "1.4.0 (2026-02)",
  "6.5.7" = "6.5.7 (2023)"
)

# Classify the model into SUP or HAC based on the model name.
model_family_from_model <- function(model) {
  case_when(
    str_detect(str_to_lower(model), "sup") ~ "SUP",
    str_detect(str_to_lower(model), "hac") ~ "HAC",
    TRUE ~ NA_character_
  )
}

# Find all telomere length files in the input directory.
files <- list.files(base_dir, pattern = "telomere_lengths.tsv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No files found.")

# Read one telomere_lengths.tsv file and add metadata parsed from its path.
parse_file <- function(f) {
  # Load the file safely and skip it if reading fails.
  df <- tryCatch({ read_tsv(f, show_col_types = FALSE) }, error = function(e) return(NULL))

  # Skip files that could not be read or do not contain a length column.
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

# Load all files and keep only Dorado/Guppy SUP and HAC data.
df_all <- map_dfr(files, parse_file) %>%
  filter(basecaller %in% c("dorado", "guppy"), model_family %in% c("SUP", "HAC"))

# Keep display labels in the order they appear in the loaded data.
ordered_levels <- unique(df_all$display_label)
df_all$display_label <- factor(df_all$display_label, levels = ordered_levels)

# Create one histogram for one version/basecaller/tool group.
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

# Create output folder for histogram PDFs.
out_root <- "histograms"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Create separate PDFs for each tool and model family combination.
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

   # Combine individual histograms into one multi-panel PDF.
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

# Save basic summary statistics for each tool/basecaller/version/model family group.
summary_stats <- df_all %>%
  group_by(tool, basecaller, version, model_family) %>%
  summarise(n = n(), median = median(length), .groups = "drop")

write_tsv(summary_stats, "histogram_summary_stats.tsv")
message("Files generated in 'histograms' folder. Each plot has full X/Y axes.")
