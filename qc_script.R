library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)

FeatureLevelData_df <- as.data.frame(quant_msstats$FeatureLevelData)

FeatureLevelData_df_wide <- FeatureLevelData_df %>% 
  dplyr::select(PSM, log2Intensity, Channel, Condition) %>% 
  tidyr::pivot_wider(names_from = PSM, values_from = log2Intensity)

FeatureLevelData_df_wide_pca <- FeatureLevelData_df_wide %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Channel") %>% 
  select(!Condition) %>% 
  select_if(~ !any(is.na(.)))

FeatureLevelData_df_wide_pca <- FeatureLevelData_df_wide_pca[ , colSums(is.na(FeatureLevelData_df_wide_pca))==0]

res.pca <- prcomp(FeatureLevelData_df_wide_pca, scale = TRUE)

fviz_eig(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

FeatureLevelData_df_wide_pca <- FeatureLevelData_df_wide[, ]

FeatureLevelData_df_wide <- FeatureLevelData_df %>% 
  dplyr::select(PSM, log2Intensity, Channel) %>% 
  dplyr::mutate(is_NA = is.na(log2Intensity)) %>% 
  dplyr::filter(is_NA == FALSE) %>% 
  tidyr::pivot_wider(names_from = Channel, values_from = log2Intensity)

FeatureLevelData_df %>% ggplot(aes(x = Channel, y = log2Intensity)) +
  geom_violin()


pca_result <- FeatureLevelData_df_wide[, 2:17] %>% 
  prcomp(scale. = TRUE)

FeatureLevelData_df_qc <- FeatureLevelData_df %>% 
  dplyr::mutate(
    is_NA = is.na(log2Intensity),
    is_non_finite = !is.finite(log2Intensity)
  ) %>% 
  dplyr::filter(
    is_non_finite == TRUE,
    is_NA == TRUE
  )

colSums(is.na(FeatureLevelData_df_wide))

FeatureLevelData_df <- FeatureLevelData_df %>% 
  dplyr::select(PSM, log2Intensity, Channel) %>% 
  dplyr::mutate(is_NA = is.na(log2Intensity)) %>% 
  dplyr::filter(is_NA != TRUE)

FeatureLevelData_df_wide <- FeatureLevelData_df %>%
  tidyr::pivot_wider(names_from = Channel, values_from = log2Intensity)

FeatureLevelData_df_wide_missing <- FeatureLevelData_df_wide %>%
  filter_at(vars(3:18), is.na(.))

FeatureLevelData_df_wide_missing <- FeatureLevelData_df_wide %>%
  rowwise() %>%
  mutate(has_non_finite = any(!is.finite(c_across(3:18))))

FeatureLevelData_df_wide_wo_missing <- FeatureLevelData_df_wide_missing %>% 
  dplyr::filter(has_non_finite != TRUE) %>% 
  dplyr::select(1, 3:18)


pca_result <- FeatureLevelData_df_wide_wo_missing[, 2:17] %>% 
  prcomp(scale. = TRUE)

rmarkdown::render_site()

dataProcessPlotsTMT(data=quant_msstats,
                    type='ProfilePlot', # choice of visualization
                    width = 21,
                    height = 7,
                    which.Protein = 'O75533') 