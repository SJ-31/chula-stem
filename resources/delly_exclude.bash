#!/usr/bin/env bash

link=$(yq ".links.delly_exclude" meta.yaml)
o=$(yq ".dir.tool_ref" meta.yaml)
tmp="${o}/delly_temp.tsv"
wget "$link" -O "$tmp"
sed 's/chr//' "$tmp" > "$o/hg38.excl.tsv"
rm "$tmp"
