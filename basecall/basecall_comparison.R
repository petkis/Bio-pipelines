#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- args[1]

# =========================
# VERSION MAPS
# =========================

version_map <- c(
  "0.9.1" = "0.9.1 (2024 early)",
  "0.9.6" = "0.9.6 (2024 late)",
  "1.0.0" = "1.0.0 (2025-05)",
  "1.2.0" = "1.2.0 (2025-10)",
  "1.3.0" = "1.3.0 (2025-11)",
  "1.4.0" = "1.4.0 (2026-02)",
  "6.0.6" = "6.0.6",
  "6.3.8" = "6.3.8",
  "6.5.7" = "6.5.7 (2023)"
)

# Short labels used only for plotting.
# This keeps the PDF x-axis readable.
version_short_map <- c(
  "0.9.1" = "0.9.1",
  "0.9.6" = "0.9.6",
  "1.0.0" = "1.0.0",
  "1.2.0" = "1.2.0",
  "1.3.0" = "1.3.0",
  "1.4.0" = "1.4.0",
  "6.0.6" = "6.0.6",
  "6.3.8" = "6.3.8",
  "6.5.7" = "6.5.7"
)

version_order <- unname(version_map)
version_short_order <- unname(version_short_map)

# =========================
# LOAD FILES
# =========================

files <- list.files(
  base_dir,
  pattern = "telomere_lengths.tsv",
  recursive = TRUE,
  full.names = TRUE
)

# =========================
# PARSE FUNCTION
# =========================

parse_file <- function(f) {

  df <- tryCatch({
    read_tsv(f, show_col_types = FALSE)
  }, error = function(e) return(NULL))

  if (is.null(df)) return(NULL)

  df$length <- suppressWarnings(as.numeric(df$length))
  df <- df[!is.na(df$length) & df$length > 0, ]

  if (nrow(df) == 0) return(NULL)

  path <- str_split(f, "/")[[1]]
  folder <- path[which(str_detect(path, "dorado_|guppy_"))]

  if (length(folder) == 0) return(NULL)

  parts <- str_split(folder, "_")[[1]]

  basecaller <- parts[1]
  version_raw <- str_remove(parts[2], "v")

  version <- ifelse(
    version_raw %in% names(version_map),
    version_map[version_raw],
    version_raw
  )

  version_short <- ifelse(
    version_raw %in% names(version_short_map),
    version_short_map[version_raw],
    version_raw
  )

  model <- paste(parts[3:length(parts)], collapse = "_")
  tool <- if (str_detect(f, "/ncrf/")) "ncrf" else "telobp"

  df$basecaller <- basecaller
  df$version <- version
  df$version_short <- version_short
  df$model <- model
  df$tool <- tool

  return(df)
}

# =========================
# MERGE
# =========================

df_all <- map_dfr(files, parse_file)

if (nrow(df_all) == 0) stop("No data loaded")

# =========================
# ADD SHORT MODEL NAME
# =========================

df_all <- df_all %>%
  mutate(
    model_short = case_when(
      str_detect(model, "fast") ~ "fast",
      str_detect(model, "hac")  ~ "hac",
      str_detect(model, "sup")  ~ "sup",
      TRUE ~ model
    )
  )

# =========================
# FACTOR ORDERING
# =========================

df_all <- df_all %>%
  mutate(
    basecaller = factor(basecaller, levels = c("dorado", "guppy")),
    version = factor(version, levels = version_order),
    version_short = factor(version_short, levels = version_short_order),
    tool = factor(tool, levels = c("telobp", "ncrf")),
    model_short = factor(model_short, levels = c("fast", "hac", "sup"))
  )

# =========================
# ADD GROUP COLUMNS
# =========================

df_all <- df_all %>%
  mutate(
    group = case_when(
      basecaller == "guppy"  & tool == "telobp" ~ "guppy - telobp",
      basecaller == "guppy"  & tool == "ncrf"   ~ "guppy - ncrf",
      basecaller == "dorado" & tool == "telobp" ~ "dorado - telobp",
      basecaller == "dorado" & tool == "ncrf"   ~ "dorado - ncrf"
    ),
    group_model = paste(basecaller, tool, model_short, sep = " - ")
  )

# =========================
# PLOTTING FUNCTIONS
# =========================

plot_violin <- function(
    data,
    title,
    filename,
    fill_var = "model",
    force_width = NULL,
    force_height = 7.5,
    y_max = NULL
) {

  if (nrow(data) == 0) {
    warning(paste("Skipping plot because no data:", filename))
    return(NULL)
  }

  # For version plots, use the shorter version labels.
  x_var <- if (fill_var == "version") "version_short" else fill_var

  p <- ggplot(data, aes(x = .data[[x_var]], y = length, fill = basecaller)) +
    geom_violin(
      trim = TRUE,
      alpha = 0.75,
      color = "black",
      linewidth = 0.35
    ) +
    geom_boxplot(
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.65,
      linewidth = 0.35
    ) +
    scale_fill_manual(
      values = c(
        "dorado" = "#F08A80",
        "guppy"  = "#41B6C4"
      ),
      drop = FALSE
    ) +
    scale_x_discrete(
      guide = guide_axis(n.dodge = 2)
    ) +
    labs(
      title = title,
      x = ifelse(fill_var == "version", "Version", "Category"),
      y = "Telomere length (bp)",
      fill = "Basecaller"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 15, margin = margin(t = 10)),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 13),
      axis.text.y = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      strip.text = element_text(size = 13, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 15, 10, 15)
    )

  # Version plots look cleaner split by basecaller.
  if (fill_var == "version" && length(unique(na.omit(data$basecaller))) > 1) {
    p <- p + facet_wrap(~ basecaller, scales = "free_x", nrow = 1)
  }

  if (!is.null(y_max)) {
    p <- p + coord_cartesian(ylim = c(0, y_max))
  } else {
    p <- p + coord_cartesian(ylim = c(0, NA))
  }

  n_groups <- length(unique(data[[x_var]]))
  plot_width <- if (!is.null(force_width)) force_width else max(10, n_groups * 1.3)

  ggsave(
    filename,
    p,
    width = plot_width,
    height = force_height,
    units = "in",
    dpi = 300
  )
}

plot_violin_log <- function(
    data,
    title,
    filename,
    fill_var = "model",
    force_width = NULL,
    force_height = 7.5
) {

  if (nrow(data) == 0) {
    warning(paste("Skipping plot because no data:", filename))
    return(NULL)
  }

  # For version plots, use the shorter version labels.
  x_var <- if (fill_var == "version") "version_short" else fill_var

  p <- ggplot(data, aes(x = .data[[x_var]], y = length, fill = basecaller)) +
    geom_violin(
      trim = TRUE,
      alpha = 0.75,
      color = "black",
      linewidth = 0.35
    ) +
    geom_boxplot(
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.65,
      linewidth = 0.35
    ) +
    scale_y_log10() +
    scale_fill_manual(
      values = c(
        "dorado" = "#F08A80",
        "guppy"  = "#41B6C4"
      ),
      drop = FALSE
    ) +
    scale_x_discrete(
      guide = guide_axis(n.dodge = 2)
    ) +
    labs(
      title = paste0(title, " (log10 scale)"),
      x = ifelse(fill_var == "version", "Version", "Category"),
      y = "Telomere length (log10 bp)",
      fill = "Basecaller"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 15, margin = margin(t = 10)),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 13),
      axis.text.y = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      strip.text = element_text(size = 13, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 15, 10, 15)
    )

  # Version plots look cleaner split by basecaller.
  if (fill_var == "version" && length(unique(na.omit(data$basecaller))) > 1) {
    p <- p + facet_wrap(~ basecaller, scales = "free_x", nrow = 1)
  }

  n_groups <- length(unique(data[[x_var]]))
  plot_width <- if (!is.null(force_width)) force_width else max(10, n_groups * 1.3)

  ggsave(
    filename,
    p,
    width = plot_width,
    height = force_height,
    units = "in",
    dpi = 300
  )
}

# =========================
# 1. SUP MODELS -> VERSION COMPARISON
# =========================

df_sup_telobp <- df_all %>%
  filter(str_detect(model, "sup"), tool == "telobp")

df_sup_ncrf <- df_all %>%
  filter(str_detect(model, "sup"), tool == "ncrf")

plot_violin(
  df_sup_telobp,
  "Version comparison: SUP models - TeloBP",
  "compare_versions_sup_telobp_violin.pdf",
  "version",
  12,
  7.5
)

plot_violin_log(
  df_sup_telobp,
  "Version comparison: SUP models - TeloBP",
  "compare_versions_sup_telobp_violin_log.pdf",
  "version",
  12,
  7.5
)

plot_violin(
  df_sup_ncrf,
  "Version comparison: SUP models - NCRF",
  "compare_versions_sup_ncrf_violin.pdf",
  "version",
  12,
  7.5
)

plot_violin_log(
  df_sup_ncrf,
  "Version comparison: SUP models - NCRF",
  "compare_versions_sup_ncrf_violin_log.pdf",
  "version",
  12,
  7.5
)

# =========================
# 1b. HAC MODELS -> VERSION COMPARISON
# =========================

df_hac_telobp <- df_all %>%
  filter(str_detect(model, "hac"), tool == "telobp")

df_hac_ncrf <- df_all %>%
  filter(str_detect(model, "hac"), tool == "ncrf")

plot_violin(
  df_hac_telobp,
  "Version comparison: HAC models - TeloBP",
  "compare_versions_hac_telobp_violin.pdf",
  "version",
  12,
  7.5
)

plot_violin_log(
  df_hac_telobp,
  "Version comparison: HAC models - TeloBP",
  "compare_versions_hac_telobp_violin_log.pdf",
  "version",
  12,
  7.5
)

plot_violin(
  df_hac_ncrf,
  "Version comparison: HAC models - NCRF",
  "compare_versions_hac_ncrf_violin.pdf",
  "version",
  12,
  7.5
)

plot_violin_log(
  df_hac_ncrf,
  "Version comparison: HAC models - NCRF",
  "compare_versions_hac_ncrf_violin_log.pdf",
  "version",
  12,
  7.5
)

# =========================
# 2. NCRF vs TELOBP for GUPPY
# =========================

df_g65 <- df_all %>%
  filter(
    basecaller == "guppy",
    version == "6.5.7 (2023)",
    str_detect(model, "sup")
  )

plot_violin(
  df_g65,
  "Guppy 6.5.7 SUP - NCRF vs TeloBP",
  "guppy65_tool_compare.pdf",
  "tool"
)

plot_violin_log(
  df_g65,
  "Guppy 6.5.7 SUP - NCRF vs TeloBP",
  "guppy65_tool_compare_log.pdf",
  "tool"
)

# =========================
# 3. NCRF vs TELOBP for DORADO
# =========================

df_d14 <- df_all %>%
  filter(
    basecaller == "dorado",
    version == "1.4.0 (2026-02)",
    str_detect(model, "sup")
  )

plot_violin(
  df_d14,
  "Dorado 1.4.0 SUP - NCRF vs TeloBP",
  "dorado14_tool_compare.pdf",
  "tool"
)

plot_violin_log(
  df_d14,
  "Dorado 1.4.0 SUP - NCRF vs TeloBP",
  "dorado14_tool_compare_log.pdf",
  "tool"
)

# =========================
# 4. DORADO MODEL COMPARISON
# =========================

df_d14_all <- df_all %>%
  filter(
    basecaller == "dorado",
    version == "1.4.0 (2026-02)"
  )

plot_violin(
  df_d14_all %>% filter(tool == "telobp"),
  "Dorado 1.4.0 models: TeloBP",
  "dorado14_models_telobp.pdf",
  "model_short"
)

plot_violin_log(
  df_d14_all %>% filter(tool == "telobp"),
  "Dorado 1.4.0 models: TeloBP",
  "dorado14_models_telobp_log.pdf",
  "model_short"
)

plot_violin(
  df_d14_all %>% filter(tool == "ncrf"),
  "Dorado 1.4.0 models: NCRF",
  "dorado14_models_ncrf.pdf",
  "model_short"
)

plot_violin_log(
  df_d14_all %>% filter(tool == "ncrf"),
  "Dorado 1.4.0 models: NCRF",
  "dorado14_models_ncrf_log.pdf",
  "model_short"
)

plot_violin(
  df_d14_all,
  "Dorado 1.4.0 models: combined",
  "dorado14_models_combined.pdf",
  "group_model"
)

plot_violin_log(
  df_d14_all,
  "Dorado 1.4.0 models: combined",
  "dorado14_models_combined_log.pdf",
  "group_model"
)

# =========================
# 5. BASECALLER COMPARISON
# =========================

df_compare <- df_all %>%
  filter(
    (basecaller == "guppy" & version == "6.5.7 (2023)") |
      (basecaller == "dorado" & version == "1.4.0 (2026-02)")
  ) %>%
  filter(str_detect(model, "sup"))

plot_violin(
  df_compare,
  "Guppy vs Dorado: combined",
  "compare_basecaller_combined.pdf",
  "group"
)

plot_violin_log(
  df_compare,
  "Guppy vs Dorado: combined",
  "compare_basecaller_combined_log.pdf",
  "group"
)

# =========================
# SAVE DATA
# =========================

write_tsv(df_all, "summary_stats.tsv")
