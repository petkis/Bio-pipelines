#!/usr/bin/env Rscript

library(ggplot2)

input_file <- "telomere_lengths.tsv"
output_pdf <- "graph.pdf"

df <- read.table(input_file, header = TRUE, sep = "\t", fill = TRUE)

df$length <- suppressWarnings(as.numeric(df$length))
df <- df[!is.na(df$length) & df$length > 0, ]

if (nrow(df) == 0) {
  stop("No valid telomere lengths")
}

# Metadata
basecaller <- Sys.getenv("BASECALLER")
version <- Sys.getenv("VERSION")
model <- Sys.getenv("MODEL")
tool <- Sys.getenv("GRAPH_TOOL")

meta <- paste(
  "Basecaller:", basecaller,
  "| Version:", version,
  "| Model:", model,
  "| Tool:", tool
)

ggplot(df, aes(x = length)) +
  geom_histogram(bins = 50,
                 fill = "steelblue",
                 color = "black",
                 alpha = 0.7) +
  labs(
    title = "Telomere length distribution",
    subtitle = meta,
    x = "Telomere length (bp)",
    y = "Count"
  ) +
  theme_minimal()

ggsave(output_pdf)
