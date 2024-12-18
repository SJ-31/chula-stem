#!/usr/bin/env bash

dir=$(yq ".dir.ref" meta.yaml)
rdir="${dir}/variants/gatk_resources"

gcloud storage cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz "${rdir}"
gcloud storage cp gs://genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf "${rdir}"
