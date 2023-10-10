# ---- libraries loading ----

library(MSstatsTMT)
library(tidyverse)
library(ggrepel)
library(plotly)
library(scales)
library(ggthemes)
library(data.table)
library(utils)

# ---- data loading ----
 
proteinGroups <- data.table::fread(
  "data/default/proteinGroups.txt"
)

evidence <- data.table::fread(
  "data/default/evidence.txt"
)

annotation <- data.table::fread(
  "data/default/annotation.tsv"
)

# ---- default data input and preparation ----
# ---- maxquant input ----

if (
  menu(
    c("Yes", "No"), 
    title = "Do you want to start from scratch?"
  ) == 1
) {
  
  input_mq <- MaxQtoMSstatsTMTFormat(
    evidence, 
    proteinGroups, 
    annotation,
    use_log_file=T,
    append=T,
    verbose=T,
    log_file_path="logs/MSstats.log"
  )
  
  fwrite(input_mq, file = "data/default/input_mq")
  
} else {
  
  input_mq <- fread(
    "data/default/input_mq"
  )
  
}

input_mq_df <- as.data.frame(input_mq)

# ---- protein summarization ----

quant_msstats <- proteinSummarization(
  data = input_mq,
  method = "msstats",
  global_norm=T,
  reference_norm=T,
  remove_norm_channel=T,
  remove_empty_channel=T,
  use_log_file=T,
  append=T,
  verbose=T,
  log_file_path="logs/MSstats.log",
  msstats_log_path="logs/MSstats.log"
)

quant_msstats_protein <- quant_msstats$ProteinLevelData

quant_msstats_feature <- quant_msstats$FeatureLevelData

fwrite(quant_msstats$ProteinLevelData, file = "data/default/quant_msstats_protein")

fwrite(quant_msstats$FeatureLevelData, file = "data/default/quant_msstats_feature")

names(quant_msstats_merged) <- c("FeatureLevelData", "ProteinLevelData")

quant_msstats_merged <- list(
  fread("data/default/quant_msstats_feature"),
  fread("data/default/quant_msstats_protein")
)

quant_msstats$ProteinLevelData

cells_comparison <- groupComparisonTMT(
  data = quant_msstats_merged,
  contrast.matrix = cells_comparison_matrix,
  adj.method = "BH",
  moderated=T,
  use_log_file=T,
  append=T,
  verbose=T,
  log_file_path = "logs/MSstats.log"
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

dataProcessPlotsTMT(data=quant_msstats,
                    type='ProfilePlot', # choice of visualization
                    width = 21,
                    height = 7,
                    which.Protein = 'Q13268') 

df <- cells_comparison_results %>% 
  dplyr::filter(
    pvalue < 0.05
  )

googlesheets4::write_sheet(
  df,
  ss = "15uju3OURyXCH4ChNIhtnFjx5ysGk08keU1xID-m8N8M"
)
