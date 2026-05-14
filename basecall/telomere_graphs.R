#!/usr/bin/env Rscript

library(ggplot2)

# Create a histogram PDF from one telomere_lengths.tsv file.
# The script reads telomere lengths, removes invalid values,
# adds run metadata from environment variables, and saves graph.pdf.

# Define input and output files.
input_file <- "telomere_lengths.tsv"
output_pdf <- "graph.pdf"

# Read telomere length table.
df <- read.table(input_file, header = TRUE, sep = "\t", fill = TRUE)

# Convert the length column to numeric and keep only positive values.
df$length <- suppressWarnings(as.numeric(df$length))
df <- df[!is.na(df$length) & df$length > 0, ]

# Stop if no valid telomere lengths remain after filtering.
if (nrow(df) == 0) {
  stop("No valid telomere lengths")
}

# Read metadata from environment variables.
basecaller <- Sys.getenv("BASECALLER")
version <- Sys.getenv("VERSION")
model <- Sys.getenv("MODEL")
tool <- Sys.getenv("GRAPH_TOOL")

# Combine metadata into one subtitle string for the plot.
meta <- paste(
  "Basecaller:", basecaller,
  "| Version:", version,
  "| Model:", model,
  "| Tool:", tool
)

# Create histogram of telomere length distribution.
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

# Save the last created plot to a PDF file.
ggsave(output_pdf)