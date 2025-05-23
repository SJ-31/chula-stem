library(tidyverse)
library(here)
library(vegan)
library(paletteer)
library(reticulate)
use_condaenv("too-predict")

data_file <- here("analyses", "h1_hmg_intensity", "h1_hmgb_intensity.csv")
outdir <- here("analyses", "h1_hmg_intensity")
skc <- import("sklearn.cluster")

add_kmeans <- function(tb, n_clusters) {
  kmeans <- skc$KMeans(n_clusters = as.integer(n_clusters))
  x <- select(tb, HMGB2_intensity_Mean, H1_intensity_Mean)
  kmeans$fit(x)
  tb %>% mutate(kmeans = as.factor(kmeans$predict(x)))
}

add_agg <- function(tb, distance_threshold = NULL, n_clusters = 3) {
  if (!is.null(distance_threshold)) {
    aggc <- skc$AgglomerativeClustering(distance_threshold = distance_threshold)
  } else {
    aggc <- skc$AgglomerativeClustering(n_clusters = as.integer(n_clusters))
  }
  x <- select(tb, HMGB2_intensity_Mean, H1_intensity_Mean)
  tb %>% mutate(agg_clust = as.factor(aggc$fit_predict(x)))
}

add_optics <- function(tb, min_samples = 100) {
  optics <- skc$OPTICS(min_samples = as.integer(min_samples))
  x <- select(tb, HMGB2_intensity_Mean, H1_intensity_Mean)
  tb %>% mutate(optics = as.factor(optics$fit_predict(x)))
}

data <- read_csv(data_file) |>
  rename_with(\(x) str_replace_all(x, " ", "_"))

## --- CODE BLOCK ---

data_longer <- data |>
  pivot_longer(-Object_No) |>
  group_by(name) |>
  mutate(bin = as.factor(ntile(value, 10))) |>
  ungroup()

against_plot <- ggplot(data, aes(
  x = HMGB2_intensity_Mean, y = H1_intensity_Mean
)) +
  geom_point()
against_plot


paired_plot <- data_longer |> ggplot(aes(x = Object_No, y = log(value), color = bin)) +
  geom_point() +
  facet_wrap(~name) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Object No.") +
  ylab("Mean Intensity")
paired_plot

twod_hist <- ggplot(data, aes(x = HMGB2_intensity_Mean, y = H1_intensity_Mean)) +
  geom_bin2d(bins = 70) +
  scale_fill_paletteer_c("ggthemes::Classic Green")
twod_hist

ggsave(here(outdir, "h1_hmgb2_two_histogram.png"), twod_hist, width = 8, height = 8)

# Try cluster objects


## data |> column_to_rownames(var = "Object_No") |> vegdist()
## distances <- gc
