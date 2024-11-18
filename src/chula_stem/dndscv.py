import click
import rpy2.robjects as ro
from rpy2.robjects.packages import InstalledSTPackage, importr

dndscv: InstalledSTPackage = importr("dndscv")


@click.command()
@click.option("-c", "--cdsfile")
@click.option("-g", "--genomefile")
@click.option("-o", "--outfile")
def buildref(cdsfile: str, genomefile: str, outfile: str = "RefCDS.rda") -> None:
    dndscv.buildref(cdsfile=cdsfile, genomefile=genomefile, outfile=outfile)


def _dndscv(tsvfile: str, refdb: str, sample_id: str = "", **kwargs):
    import pandas as pd
    from rpy2.robjects import pandas2ri

    tsv = pd.read_csv(tsvfile, sep="\t")
    if "sampleID" not in tsv.columns and sample_id:
        tsv["sampleID"] = sample_id
    tsv: pd.DataFrame = tsv.loc[:, ["sampleID", "chr", "pos", "ref", "alt"]]
    tsv = tsv.drop_duplicates(subset=["chr", "pos", "ref", "alt"])
    with (ro.default_converter + pandas2ri.converter).context():
        mutations = ro.conversion.get_conversion().py2rpy(tsv)
        if kwargs:
            dndscv.dndscv(mutations, refdb=refdb, **kwargs)
        else:
            dndscv.dndscv(mutations, refdb=refdb)
