import click
from rpy2.robjects.packages import InstalledSTPackage, importr

dndscv: InstalledSTPackage = importr("dndscv")


@click.command()
@click.option("-c", "--cdsfile")
@click.option("-g", "--genomefile")
@click.option("-o", "--outfile")
def buildref(cdsfile: str, genomefile: str, outfile: str = "RefCDS.rda") -> None:
    dndscv.buildref(cdsfile, genomefile, outfile)
