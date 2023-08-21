# ---- libraries attaching ----
library(MSstatsTMT)
library(tidyverse)
library(ggrepel)
# ---- data attaching ----
proteinGroups <- read.table(
  "data/default/proteinGroups.txt", sep="\t", header=TRUE
)

evidence <- read.table(
  "data/default/evidence.txt", sep="\t", header=TRUE
)

annotation <- read.table(
  "data/default/annotation.tsv", sep="\t", header=TRUE
)
# ---- MSstats input ----
input_mq <- MaxQtoMSstatsTMTFormat(
  evidence, 
  proteinGroups, 
  annotation,
  use_log_file=T,
  append=T,
  verbose=T,
  log_file_path = "logs/MSstats.log"
)

input_mq_df <- as.data.frame(input_mq)

quant_msstats <- proteinSummarization(
  data = input_mq_df,
  method = "msstats",
  global_norm=T,
  reference_norm=T,
  remove_norm_channel=T,
  remove_empty_channel=T,
  use_log_file=T,
  append=T,
  verbose=T,
  log_file_path = "logs/MSstats.log",
  msstats_log_path ="logs/MSstats.log"
)
# ---- samples comparison ----
# Check the conditions in the protein level data
sample_types <- levels(quant_msstats$ProteinLevelData$Condition)
# Only compare condition 0.125 and 1
cells_comparison_matrix <- matrix(c(-1, 1, 0, 0), nrow = 1)
# Set the names of each row
row.names(cells_comparison_matrix) <- "cells_DMSO_vs_PlB"
# Set the column names
colnames(cells_comparison_matrix) <- sample_types

cells_comparison_matrix

cells_comparison <- groupComparisonTMT(
  data = quant_msstats,
  contrast.matrix = cells_comparison_matrix,
  adj.method = "BH",
  moderated=T,
  use_log_file=T,
  append=T,
  verbose=T,
  log_file_path = "logs/MSstats.log"
)

head(cells_comparison$ComparisonResult)

cells_comparison_results <- cells_comparison$ComparisonResult

cells_comparison_results <- left_join(
  cells_comparison_results,
  select(evidence, Proteins, Gene.names),
  by = c("Protein" = "Proteins")
) %>% distinct()

cells_comparison_results <- cells_comparison_results %>% 
  dplyr::mutate(
    color = case_when(
      pvalue < 0.01 & log2FC > 0 ~ "#488f31",
      pvalue < 0.01 & log2FC < 0 ~ "#de425b",
      pvalue < 0.05 & log2FC > 0 ~ "#a8c162",
      pvalue < 0.05 & log2FC < 0 ~ "#f9a160",
      TRUE ~ "grey"
    ),
    label = case_when(
      pvalue < 0.05 ~ paste(Gene.names, Protein, sep = "\n")
    )
  )

cells_comparison_results %>% 
  ggplot(aes(
    x = log2FC,
    y = -log10(pvalue),
    color = color,
    label = label
  )) +
  geom_point() +
  scale_colour_identity() +
  geom_text_repel(size = 3, force_pull = 5)
