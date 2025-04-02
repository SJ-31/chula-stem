library(tidyverse)
library(here)
library(paletteer)

data_file <- here("analyses", "h1_hmg_intensity", "h1_hmgb_intensity.csv")

data <- read_csv(data_file) |>
  rename_with(\(x) str_replace_all(x, " ", "_"))
## --- CODE BLOCK ---

data_longer <- data |>
  pivot_longer(-Object_No) |>
  group_by(name) |>
  mutate(bin = as.factor(ntile(value, 10))) |>
  ungroup()

against_plot <- ggplot(data, aes(x = HMGB2_intensity_Mean, y = H1_intensity_Mean)) +
  geom_point()

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
