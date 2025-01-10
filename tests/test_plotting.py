#!/usr/bin/env ipython
import chula_stem.plotting as plo

cns = "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns"
cnr = "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.cnr"


def test_plot_cnvkit():
    plotdir = "/home/shannc/Bio_SDD/chula-stem/tests/plots"
    args = {
        "cns": "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns",
        "cnr": "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.cnr",
        "name": "cnvkit",
        "loci": "chr1,chr2",
        "output": f"{plotdir}/cnvkit.png",
    }
    sizes = {"width": 10, "dpi": 300, "height": 5}
    plo.plot_main(args, sizes)
    args["loci"] = None
    sizes = {"width": 20, "dpi": 300, "height": 30}
    args["output"] = f"{plotdir}/cnvkit_full.png"
    plo.plot_main(args, sizes)
