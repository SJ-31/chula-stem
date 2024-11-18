#!/usr/bin/env python

import rpy2.robjects as ro
from rpy2.robjects.packages import InstalledSTPackage, importr

ascat: InstalledSTPackage = importr("ASCAT")


# First need to use the hts data to get logR and B-allele frequency (BAF) files
# C
#
# For HTS data (WGS, WES and TS), gamma must be set to 1 in ascat.runASCAT (author instructions)
def _run_ascat(
    tumor_bam: str, normal_bam: str, sample_name: str, BED_file: str | None = None
):
    """
    For exome data, BED_file should be a target or baits file covering the sequenced regions
    """
    ascat.prepareHTS(
        tumourseqfile=tumor_bam,
        tumourname=sample_name,
        normalseqfile=normal_bam,
        BED_file=BED_file,
    )
    loaded = ascat.loadData(
        Tumor_LogR_file=f"{sample_name}_tumourLogR",
        Tumor_BAF_file=f"{sample_name}_tumourBAF",
        Germline_LogR_file=f"{sample_name}_normalLogR",
        Germline_BAF_file=f"{sample_name}_normalBAF",
    )
    segmented = ascat.aspcf(
        loaded
    )  # Alelle-specific piecewise constant fitting (ASPCF) algorithm segments the data
    output = ascat.runAscat(segmented)


# Will output four files
# - <sample_name>_tumourLogR.txt
# - <sample_name>_normalLogR.txt
# - <sample_name>_tumourBAF.txt
# - <sample_name>_normalBAF.txt
