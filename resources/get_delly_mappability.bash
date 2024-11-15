#!/usr/bin/env bash

source util.bash
out="$(conf ".dir.tool_ref")/delly_mappability_GRCh38.fna"
link=$(conf ".link.delly_mappability_GRCh38")
wget "$link" -O "$out"
bgzip "$out"
samtools faidx "$out"
