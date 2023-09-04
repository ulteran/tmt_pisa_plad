ProteinLevelData_df <- as.data.frame(quant_msstats$ProteinLevelData)

ProteinLevelData_df_wide <- ProteinLevelData_df %>% 
  dplyr::select(Channel, Protein, Abundance) %>% 
  tidyr::pivot_wider(names_from = Protein, values_from = Abundance)

ProteinLevelData_df_wide_pca <- ProteinLevelData_df_wide %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Channel")

ProteinLevelData_df_wide_pca <- ProteinLevelData_df_wide_pca[ , colSums(is.na(ProteinLevelData_df_wide_pca))==0]

# Finding distance matrix
ProteinLevelData_df_wide_pca_mat <- dist(ProteinLevelData_df_wide_pca, method = 'euclidean')
distance_mat

# Fitting Hierarchical clustering Model
# to training dataset
set.seed(240)  # Setting seed
Hierar_cl <- hclust(ProteinLevelData_df_wide_pca_mat, method = "average")
Hierar_cl

# Plotting dendrogram
plot(Hierar_cl)

# Choosing no. of clusters
# Cutting tree by height
abline(h = 110, col = "green")

# Cutting tree by no. of clusters
fit <- cutree(Hierar_cl, k = 3 )
fit

table(fit)
rect.hclust(Hierar_cl, k = 3, border = "green")

ProteinLevelData_df_wide <- ProteinLevelData_df %>% 
  dplyr::select(Channel, Protein, Abundance) %>% 
  tidyr::pivot_wider(names_from = Channel, values_from = Abundance)

write.csv(ProteinLevelData_df_wide, "ProteinLevelData.csv")
