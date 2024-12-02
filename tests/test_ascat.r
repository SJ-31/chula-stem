library(ASCAT)

ascat.prepareHTS(
  tumourseqfile = "/data/project/stemcell/shannc/tests/wes_test/ascat/c_recal.bam",
  normalseqfile = "/data/project/stemcell/shannc/tests/wes_test/ascat/b_recal.bam",
  tumourname = "patient10",
  allelecounter_exe = "~/miniforge3/bin/alleleCounter",
  loci.prefix = "/data/project/stemcell/shannc/tests/wes_test/ascat/G1000_loci_hg38/G1000_loci_hg38_chr",
  alleles.prefix = "/data/project/stemcell/shannc/tests/wes_test/ascat/G1000_alleles_hg38/G1000_alleles_hg38_chr"
)
# Wasn't working because alleleCounter can't seem to detect the files...
