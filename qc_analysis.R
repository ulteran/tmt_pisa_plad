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
library(romics)

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

# ---- 

evidence %>% 
  ggplot(aes(x = PEP)) +
  geom_histogram()

df_evidence <- evidence %>% 
  filter(if_any(58:73, ~. == 0 | is.na(.)))

df_evidence <- evidence %>% 
  dplyr::select(Sequence, Modifications, "m/z", Fraction, 58:73)
  
# ----
romics::hello()
