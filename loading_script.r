# ---- libraries loading ----

library(MSstatsTMT)
library(tidyverse)
library(ggrepel)
library(plotly)
library(scales)
library(ggthemes)
library(data.table)
library(utils)
library(extrafont)

#---- utils script ----

# Custom transformation
custom_reverse_log_trans <- function(base = exp(1)) {
  transform <- function(x) -log(x, base)
  inverse <- function(x) base^(-x)
  scales::trans_new(
    paste("reverselog-", format(base), sep = ""),
    transform,
    inverse,
    log_breaks(base = base),
    domain = c(1e-1000, Inf)
  )
}

font_import(paths = c("utils/fonts/", prompt = F))

custom <- list(
  ggthemes::theme_base(),
  theme(
    legend.position = "bottom",
    text = element_text(family = "Google Sans", size = 12),
    plot.background = element_rect(color = NA)
  )
)


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
# maxquant input

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
    which.proteinid = "Leading.razor.proteins",
    use_log_file = TRUE,
    append = TRUE,
    verbose = TRUE,
    log_file_path = "logs/MSstats.log"
  )
  
  fwrite(input_mq, file = "data/default/std_input_mq")
  
} else {
  
  input_mq <- fread(
    "data/default/input_mq"
  )
  
}

input_mq_df <- as.data.frame(input_mq)

head(input_mq_df)


# ---- protein summarization ----

quant_msstats <- proteinSummarization(
  data = input_mq,
  method = "msstats",
  global_norm = T,
  reference_norm = T,
  remove_norm_channel = T,
  remove_empty_channel = T,
  use_log_file = T,
  append = T,
  verbose = T,
  log_file_path = "logs/MSstats.log",
  msstats_log_path = "logs/MSstats.log"
)

quant_msstats_protein <- quant_msstats$ProteinLevelData

quant_msstats_feature <- quant_msstats$FeatureLevelData

fwrite(quant_msstats$ProteinLevelData, file = "data/default/std_quant_msstats_protein")

fwrite(quant_msstats$FeatureLevelData, file = "data/default/std_quant_msstats_feature")

names(quant_msstats_merged) <- c("FeatureLevelData", "ProteinLevelData")

quant_msstats_merged <- list(
  fread("data/default/quant_msstats_feature"),
  fread("data/default/quant_msstats_protein")
)

# ---- samples comparison ----

sample_types <- levels(quant_msstats$ProteinLevelData$Condition)

cells_comparison_matrix <- matrix(c(-1, 1, 0, 0), nrow = 1)

row.names(cells_comparison_matrix) <- "cells_DMSO_vs_Pld-B"

colnames(cells_comparison_matrix) <- sample_types

cells_comparison_matrix

cells_comparison <- groupComparisonTMT(
  data = quant_msstats,
  contrast.matrix = cells_comparison_matrix,
  adj.method = "BH",
  moderated = T,
  use_log_file = T,
  append= T,
  verbose = T,
  log_file_path = "logs/MSstats.log"
)

cells_comparison_results <- cells_comparison$ComparisonResult

view(cells_comparison_results)

fwrite(
  cells_comparison_results, 
  file = "data/default/std_cells_comparison_results.csv"
)

dataProcessPlotsTMT(data = quant_msstats,
                    type = 'ProfilePlot', # choice of visualization
                    width = 21,
                    height = 7,
                    which.Protein = 'Q13268')
