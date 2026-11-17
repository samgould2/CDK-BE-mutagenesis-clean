"""Download every structure used by the CDK figures into ./structures/.

Kept local (rather than letting ChimeraX fetch) so the figures are reproducible
offline and byte-identical between runs.
"""
import os, sys, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cdk_config import STRUCTURES

DEST = os.path.join(HERE, "structures")
os.makedirs(DEST, exist_ok=True)

for gene, cfg in STRUCTURES.items():
    fn = cfg["file"]
    path = os.path.join(DEST, fn)
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        print(f"{gene:6s} have {fn}")
        continue
    # CDK19 has no crystal structure in these panels: they use the AlphaFold
    # v4 model, which EBI no longer serves. That file is committed under
    # structures/ for exactly this reason, so reaching this line for CDK19
    # means it has been deleted and there is nowhere left to fetch it from.
    url = cfg["url"]
    print(f"{gene:6s} fetching {url}")
    urllib.request.urlretrieve(url, path)
    print(f"       -> {fn} ({os.path.getsize(path)/1e6:.1f} MB)")
