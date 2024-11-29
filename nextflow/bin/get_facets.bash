#!/usr/bin/env bash

Rscript -e "
data <- readRDS(\"${1}\")
write.table(data\$segs, file = \"$2\", sep = \"\t\", row.names = FALSE, quote = FALSE)
"
