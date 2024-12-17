#!/usr/bin/env bash

dir=$(yq ".dir.ref" meta.yaml)
rdir="${dir}/variants/gatk_resources"

# Download resources from gcloud for BQSR, following recommendations in https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-faqs/What_should_I_use_as_known_variants_sites_for_running_tool_X%3F.md

gcloud storage cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz "${rdir}"
gcloud storage cp gs://genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf "${rdir}"
