#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(readr)

args <- commandArgs(trailingOnly=TRUE)
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
# LOAD FILES
# =========================

files <- list.files(
  base_dir,
  pattern="telomere_lengths.tsv",
  recursive=TRUE,
  full.names=TRUE
)

# =========================
# PARSE FUNCTION
# =========================

parse_file <- function(f) {

  df <- tryCatch({
    read_tsv(f, show_col_types = FALSE)
  }, error=function(e) return(NULL))

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

  version <- ifelse(version_raw %in% names(version_map),
                    version_map[version_raw],
                    version_raw)

  model <- paste(parts[3:length(parts)], collapse="_")
  tool <- if (str_detect(f, "/ncrf/")) "ncrf" else "telobp"

  df$basecaller <- basecaller
  df$version <- version
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
# ADD GROUP COLUMN
# =========================

df_all <- df_all %>%
  mutate(
    group = case_when(
      basecaller=="guppy"  & tool=="telobp" ~ "guppy - telobp",
      basecaller=="guppy"  & tool=="ncrf"   ~ "guppy - ncrf",
      basecaller=="dorado" & tool=="telobp" ~ "dorado - telobp",
      basecaller=="dorado" & tool=="ncrf"   ~ "dorado - ncrf"
    )
  )

# =========================
# PLOTTING FUNCTIONS
# =========================

plot_violin <- function(data, title, filename, fill_var="model", force_width=NULL) {

  p <- ggplot(data, aes_string(x=fill_var, y="length", fill="basecaller")) +
    geom_violin(trim=FALSE, alpha=0.7, color="black") +
    geom_boxplot(width=0.12, outlier.size=0.4, alpha=0.5) +
    coord_cartesian(ylim = c(0, NA)) +
    scale_fill_manual(values=c(
      "dorado" = "#F08A80",
      "guppy"  = "#41B6C4"
    )) +
    theme_minimal(base_size = 14) +
    labs(
      title=title,
      x="Category",
      y="Telomere length (bp)",
      fill="basecaller"
    ) +
    theme(
      axis.text.x = element_text(angle=45, hjust=1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )

  n_groups <- length(unique(data[[fill_var]]))
  plot_width <- if (!is.null(force_width)) force_width else max(10, n_groups * 1.5)

  ggsave(filename, p, width=plot_width, height=6, units="in", dpi=300)
}

plot_violin_log <- function(data, title, filename, fill_var="model", force_width=NULL) {

  p <- ggplot(data, aes_string(x=fill_var, y="length", fill="basecaller")) +
    geom_violin(trim=FALSE, alpha=0.7, color="black") +
    geom_boxplot(width=0.12, outlier.size=0.4, alpha=0.5) +
    scale_y_log10() +
    scale_fill_manual(values=c(
      "dorado" = "#F08A80",
      "guppy"  = "#41B6C4"
    )) +
    theme_minimal(base_size = 14) +
    labs(
      title=paste0(title, " (log10 scale)"),
      x="Category",
      y="Telomere length (log10 scale)",
      fill="basecaller"
    ) +
    theme(
      axis.text.x = element_text(angle=45, hjust=1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )

  n_groups <- length(unique(data[[fill_var]]))
  plot_width <- if (!is.null(force_width)) force_width else max(10, n_groups * 1.5)

  ggsave(filename, p, width=plot_width, height=6, units="in", dpi=300)
}

# =========================
# 1. SUP MODELS (SPLIT BY TOOL)
# =========================

df_sup_telobp <- df_all %>%
  filter(str_detect(model, "sup"), tool=="telobp")

df_sup_ncrf <- df_all %>%
  filter(str_detect(model, "sup"), tool=="ncrf")

plot_violin(df_sup_telobp,
            "Version comparison (SUP models - TeloBP)",
            "compare_versions_sup_telobp_violin.pdf",
            "version", 22)

plot_violin_log(df_sup_telobp,
                "Version comparison (SUP models - TeloBP)",
                "compare_versions_sup_telobp_violin_log.pdf",
                "version", 22)

plot_violin(df_sup_ncrf,
            "Version comparison (SUP models - NCRF)",
            "compare_versions_sup_ncrf_violin.pdf",
            "version", 22)

plot_violin_log(df_sup_ncrf,
                "Version comparison (SUP models - NCRF)",
                "compare_versions_sup_ncrf_violin_log.pdf",
                "version", 22)

# =========================
# 1b. HAC MODELS (SPLIT + LOG ADDED)
# =========================

df_hac_telobp <- df_all %>%
  filter(str_detect(model, "hac"), tool=="telobp")

df_hac_ncrf <- df_all %>%
  filter(str_detect(model, "hac"), tool=="ncrf")

# TELOBP
plot_violin(df_hac_telobp,
            "Version comparison (HAC models - TeloBP)",
            "compare_versions_hac_telobp_violin.pdf",
            "version", 22)

plot_violin_log(df_hac_telobp,
                "Version comparison (HAC models - TeloBP)",
                "compare_versions_hac_telobp_violin_log.pdf",
                "version", 22)

# NCRF
plot_violin(df_hac_ncrf,
            "Version comparison (HAC models - NCRF)",
            "compare_versions_hac_ncrf_violin.pdf",
            "version", 22)

plot_violin_log(df_hac_ncrf,
                "Version comparison (HAC models - NCRF)",
                "compare_versions_hac_ncrf_violin_log.pdf",
                "version", 22)

# =========================
# 2. NCRF vs TELOBP (GUPPY)
# =========================

df_g65 <- df_all %>%
  filter(basecaller=="guppy",
         version=="6.5.7 (2023)",
         str_detect(model, "sup"))

plot_violin(df_g65, "Guppy 6.5.7 SUP - NCRF vs TeloBP",
            "guppy65_tool_compare.pdf", "tool")

plot_violin_log(df_g65, "Guppy 6.5.7 SUP - NCRF vs TeloBP",
                "guppy65_tool_compare_log.pdf", "tool")

# =========================
# 3. NCRF vs TELOBP (DORADO)
# =========================

df_d14 <- df_all %>%
  filter(basecaller=="dorado",
         version=="1.4.0 (2026-02)",
         str_detect(model, "sup"))

plot_violin(df_d14, "Dorado 1.4.0 SUP - NCRF vs TeloBP",
            "dorado14_tool_compare.pdf", "tool")

plot_violin_log(df_d14, "Dorado 1.4.0 SUP - NCRF vs TeloBP",
                "dorado14_tool_compare_log.pdf", "tool")

# =========================
# 4. DORADO MODELS
# =========================

df_d14_all <- df_all %>%
  filter(basecaller=="dorado",
         version=="1.4.0 (2026-02)")

plot_violin(df_d14_all %>% filter(tool=="telobp"),
            "Dorado 1.4.0 models (TeloBP)",
            "dorado14_models_telobp.pdf")

plot_violin(df_d14_all %>% filter(tool=="ncrf"),
            "Dorado 1.4.0 models (NCRF)",
            "dorado14_models_ncrf.pdf")

plot_violin(df_d14_all,
            "Dorado 1.4.0 models (combined)",
            "dorado14_models_combined.pdf",
            "group")

plot_violin_log(df_d14_all,
                "Dorado 1.4.0 models (combined)",
                "dorado14_models_combined_log.pdf",
                "group")

# =========================
# 5. BASECALLER COMPARISON
# =========================

df_compare <- df_all %>%
  filter(
    (basecaller=="guppy" & version=="6.5.7 (2023)") |
    (basecaller=="dorado" & version=="1.4.0 (2026-02)")
  ) %>%
  filter(str_detect(model, "sup"))

plot_violin(df_compare, "Guppy vs Dorado (combined)",
            "compare_basecaller_combined.pdf",
            "group")

plot_violin_log(df_compare, "Guppy vs Dorado (combined)",
                "compare_basecaller_combined_log.pdf",
                "group")

# =========================
# SAVE DATA
# =========================

write_tsv(df_all, "summary_stats.tsv")
