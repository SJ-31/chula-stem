library(reticulate)
library(tidyverse)
library(here)
pre <- import("sklearn.preprocessing")

# Check that R scale is equivalent to StandardScaler
cup_ai_test <- read_csv(here("extras", "CUP-AI-Dx", "data", "ExternalDataMeta.csv"))
labels <- cup_ai_test$`tumor.type...20533`
cup_ai_test <- cup_ai_test |> select(where(is.numeric))

old_cols <- colnames(cup_ai_test)
colnames(cup_ai_test) <- seq(1, ncol(cup_ai_test))

transposed <- t(cup_ai_test)
scaler <- pre$StandardScaler()
scaler$fit(transposed)
sk_result <- scaler$transform(transposed)
r_result <- scale(transposed)

back <- t(r_result) |>
  `colnames<-`(old_cols) |>
  as_tibble() |>
  mutate(tumor_type = labels)
write_csv(back, here("extras", "CUP-AI-Dx", "data", "ExternalDataMetaCustom.csv"))

## <2025-02-04 Tue> You tested to confirm that "scale" should be the way to go, and
## doing this on their data works just fine
