from . import align
from . import coverage

def calcAbun(config: str, outdir: str):
    results = align.AlignBySamp(config, outdir)
    results += coverage.TPMBySamp(config, outdir)
    return results

