```{r load}
# Install/load package
if(!require(geomorph)) install.packages("geomorph")
library(geomorph)
library(ggplot2)
if(!require(reshape2))install.packages("reshape2")
library(reshape2)
if(!require(abind))install.packages("abind")
library(abind)
```

```{r input}
# ------------------------------
# 1. Input your TPS file here
# ------------------------------
tps_files <- c("AFem_rand.TPS", "AMasc_rand.TPS", "5fem_rand.tps", "5masc_rand.TPS")  # add all your TPS files
```

```{r list}
all_gpa <- list()
all_centroid_sizes <- list()
all_pca <- list()
all_mean_shapes <- list()
```

```{r loop}
for(file in tps_files){
  
  # Read landmarks
  coords <- readland.tps(file, specID = NULL, readcurves = TRUE)
  
  # Perform GPA
  gpa <- gpagen(coords, ProcD = TRUE, print.progress = FALSE)
  
  # Store GPA
  all_gpa[[file]] <- gpa
  
  # Store centroid sizes
  all_centroid_sizes[[file]] <- gpa$Csize
  
  # PCA
  pca <- gm.prcomp(gpa$coords)
  all_pca[[file]] <- pca
  
  # Mean shape
  all_mean_shapes[[file]] <- mshape(gpa$coords)
  
  # Pairwise distances
  num_specimens <- dim(gpa$coords)[3]
  dist_matrix <- matrix(0, nrow=num_specimens, ncol=num_specimens)
  for(i in 1:num_specimens){
    for(j in i:num_specimens){
      diff <- gpa$coords[,,i] - gpa$coords[,,j]
      dist_matrix[i,j] <- sqrt(sum(diff^2))
      dist_matrix[j,i] <- dist_matrix[i,j]
    }
  }
  specimen_names <- paste0("Specimen", 1:num_specimens)
  rownames(dist_matrix) <- specimen_names
  colnames(dist_matrix) <- specimen_names
  write.csv(dist_matrix, paste0("Distances_", tools::file_path_sans_ext(file), ".csv"), row.names = TRUE)
  
}
```

```{r centroid}
library(ggplot2)
library(ggpubr)
library(dplyr)
library(emmeans)
library(stringr)
centroid_df <- centroid_df %>%
  mutate(StageSex = factor(paste(Stage, Sex), 
                           levels = c("5th Instar Female", "5th Instar Male",
                                      "Adult Female", "Adult Male")))
# Load necessary libraries
library(emmeans)

# ------------------------
# 1. Two-way ANOVA
# ------------------------
anova_res <- aov(Size ~ Stage * Sex, data = centroid_df)

# Show ANOVA table
summary(anova_res)

# ------------------------
# 2. Tukey HSD post-hoc pairwise comparisons
# ------------------------
# Get pairwise comparisons for all Stage * Sex combinations
tukey_res <- emmeans(anova_res, pairwise ~ Stage * Sex, adjust = "tukey")

# Extract results
pairwise_df <- as.data.frame(tukey_res$contrasts)

# Optional: show significance symbols
pairwise_df <- pairwise_df %>%
  dplyr::mutate(
    p.signif = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# View the pairwise comparisons
pairwise_df


# Pairwise comparisons with Tukey correction
emmeans_res <- emmeans(anova_res, pairwise ~ Stage * Sex, adjust = "tukey")
pairwise_df <- as.data.frame(emmeans_res$contrasts)

# Split contrast into group1 and group2
pairwise_df <- pairwise_df %>%
  mutate(
    group1 = str_extract(contrast, "^[^\\-]+") %>% str_trim(),
    group2 = str_extract(contrast, "[^\\-]+$") %>% str_trim(),
    p.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    y.position = max(centroid_df$Size) * 1.05 + row_number() * 0.2
  )

# Plot with asterisks
ggplot(centroid_df, aes(x = StageSex, y = Size, fill = Sex)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Centroid Size") +
  ggtitle("Centroid Size by Stage and Sex") +
  stat_pvalue_manual(
    pairwise_df,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    tip.length = 0.02,
    hide.ns = TRUE
  )
```
```{r PCA}
library(ggplot2)

# Custom names for datasets
custom_names <- c(
  "5fem_rand.tps" = "5th Instar Female",
  "5masc_rand.TPS" = "5th Instar Male",
  "AFem_rand.TPS" = "Adult Female",
  "AMasc_rand.TPS" = "Adult Male"
)

# Combine PCA scores
pca_scores <- do.call(rbind, lapply(names(all_pca), function(name){
  pca <- all_pca[[name]]
  data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    Dataset = name,
    Specimen = 1:nrow(pca$x)
  )
}))

# Apply custom labels
pca_scores$Dataset <- factor(pca_scores$Dataset, levels = names(custom_names), labels = custom_names)

# Get % variance explained (using first PCA object as example)
var_PC1 <- round(100 * all_pca[[1]]$d[1]^2 / sum(all_pca[[1]]$d^2), 1)
var_PC2 <- round(100 * all_pca[[1]]$d[2]^2 / sum(all_pca[[1]]$d^2), 1)

# Plot
ggplot(pca_scores, aes(x=PC1, y=PC2, color=Dataset, label=Specimen)) +
  geom_point(size=3) +
  geom_text(vjust=-0.5, hjust=-0.2, check_overlap=TRUE) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", var_PC1, "%)"),
    y = paste0("PC2 (", var_PC2, "%)"),
  )
```

```{r GPA}
# Custom names for datasets
custom_names <- c(
  "5fem_rand.tps" = "5th Instar Female",
  "5masc_rand.TPS" = "5th Instar Male",
  "AFem_rand.TPS" = "Adult Female",
  "AMasc_rand.TPS" = "Adult Male"
)

# Set up plotting area
par(mfrow=c(1, length(all_gpa)))  # one row, multiple columns

for(file in names(all_gpa)){
  plotAllSpecimens(all_gpa[[file]]$coords, mean=TRUE, label=TRUE)
  title(paste(custom_names[file]))  # use custom names instead of file
}

# Reset plotting layout
par(mfrow=c(1,1))

```
```{r mean shape}
library(ggplot2)
library(dplyr)
library(tidyr)

# Suppose all_mean_shapes is a list of mean shapes
mean_shapes_df <- bind_rows(lapply(names(all_mean_shapes), function(name){
  coords <- all_mean_shapes[[name]]  # [landmark, dimension]
  
  # For connected paths, create an order column (landmark order)
  data.frame(
    X = coords[,1],
    Y = coords[,2],
    Landmark = 1:nrow(coords),
    Dataset = name
  )
}))

ggplot(mean_shapes_df, aes(x = X, y = Y, group = Dataset, color = Dataset)) +
  geom_point(size = 3) +                  # plot landmarks
  geom_path(size = 1, alpha = 0.7) +      # connect landmarks in correct order
  coord_equal() +                          # preserve aspect ratio
  theme_minimal() +
  labs(title = "Comparison of Mean Shapes Across Datasets") +
  theme(legend.position = "bottom")


geom_text(aes(label = Landmark), hjust = -0.2, vjust = -0.2, size = 3)

# Original column is Dataset
mean_shapes_df$Dataset <- factor(mean_shapes_df$Dataset,
                            levels = c("5fem_rand.tps","5masc_rand.TPS", "AFem_rand.TPS", "AMasc_rand.TPS"),
                            labels = c("5th instar Female","5th instar Male", "Adult Female", "Adult Male"))

ggplot(mean_shapes_df, aes(x = X, y = Y , group = Dataset)) +
  geom_point(size = 3) +
  facet_wrap(~Dataset) +
  coord_equal() +
  theme_minimal() +
```
