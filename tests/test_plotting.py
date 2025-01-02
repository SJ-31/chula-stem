#!/usr/bin/env ipython
import chula_stem.plotting as plo


def test_plot_cnvkit():
    args = {
        "cns": "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns",
        "cnr": "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.cnr",
        "name": "cnvkit",
        "loci": "chr1,chr2",
        "output": "test_plot_cnvkit.png",
    }
    sizes = {"width": 30, "dpi": 500, "height": 10}
    plo.main(args, sizes)
