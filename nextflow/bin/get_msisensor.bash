while getopts "i:o:g:s:" f; do
    case "$f" in
        i) input=${OPTARG} ;; # Output microsatellite information ("*_unstable") from msisensorpro
        o) output=${OPTARG} ;; # Output file
        s) summary=${OPTARG} ;;
        g) gff=${OPTARG} ;;
        *) echo "Flag not recognized"
           exit 1 ;;
    esac
done

# Convert results into a bed file
awk 'BEGIN { OFS = "\t"}; FNR != 1 {print $1,$2,(length($5) * $4) + $2}' "${input}" | \
    sort-bed - > msi_sites.bed
# Can merge this site with the *_unstable file by chr and location

# Extract overlapping transcripts from the given gff file using bedtools
# - This results in finding the repetitive sequences contained within transcript ranges
echo -e "chromosome\tlocation\tstop\tgene_start\tgene_stop\tgene_name" > overlapping.tsv

awk '{ if ($3 == "transcript") { print }}' "${gff}" | gff2bed | sort - | \
    bedtools intersect -a msi_sites.bed -b stdin -loj | \
    awk 'BEGIN { OFS = "\t" }
    {
        split($13, splits, ";")
        gsub("gene_name=", "", splits[3])
        print $1,$2,$3,$5,$6,splits[3]
    }' | awk -F "\t" 'BEGIN { OFS = "\t" }
    {
        combined = $1"-"$2"-"$3
        if (!_[combined]++ && $4 != -1)
            { print }
    }' >> overlapping.tsv
# Remove 'gene_name' string with awk and retrieve unique microsatellite locations, as some
# duplicates (alternative transcript names) were reported in the gff

Rscript -e "
library('tidyverse')
u <- read_tsv(\"${input}\")
summary <- read_tsv(\"${summary}\")
if (summary\$Number_of_Unstable_Sites > 0) {
    u <- u |> mutate(stop = nchar(repeat_unit_bases) * repeat_times + location)
    u\$chromosome <- as.character(u\$chromosome)
    o <- read_tsv('overlapping.tsv') |> select(-stop)
    o\$chromosome <- as.character(o\$chromosome)
    o\$location <- as.numeric(o\$location)
    result <- left_join(u, o) |> rename(start = location) |>
           relocate(all_of(c('stop', 'gene_start', 'gene_stop', 'gene_name')), .after = start)
} else {
  cn <- c('chromosome', 'start', 'stop', 'gene_start', 'gene_stop', 'gene_name', 'left_flank_bases', 'repeat_times', 'repeat_unit_bases', 'right_flank_bases', 'pro_p', 'pro_q', 'CovReads', 'threshold')
  empty <- as.list(rep(NA, length(cn)))
  names(empty) <- cn
  result <- as_tibble(empty) |> slice(2)
}
write_tsv(result, \"${output}\")
"
