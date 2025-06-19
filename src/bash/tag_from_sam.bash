#!/usr/bin/env bash

file="${1}"
record_type="${2}" # The two-letter header record type, beginning each header line
# and starting with @
tag="${3}" # The two-character string whose value we want
samtools head "${1}" | grep "@${record_type}" | \
    sed -E "s/.*$tag:([-A-Za-z0-9_,]+).*/\1/" | \
    uniq
