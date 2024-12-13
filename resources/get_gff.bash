link=$(yq ".links.grch38_gff" meta.yaml)
o=$(yq ".dir.tool.ref")
name="GCF_000001405.40_GRCh38.p14_genomic.gff"
output="${o}/genomes/${name}"
rename=$(yq ".rename.grch38_nums" meta.yaml)
filtered="${o}/genomes/GRCh38.p14_filtered.gff"

wget "${link}" -O "${output}.gz"
gunzip "${output}.gz"
gffread "${output}" -m "${rename}" | awk '{
    if ($1 ~ /^#/) { print }
    else if (length($1) <= 2) { print };
}' > "${filtered}"
