"""
Taxon-stratified deep homolog search for the screened CDKs.

A single BLAST against all of UniProtKB saturates its 1000-hit cap inside the
vertebrates, so the deep branches of the tree are never reached. We therefore
run one search per taxonomic division, which forces the search to report the
best homologs that exist in fungi, plants, protists, archaea and bacteria as
well as animals.
"""
import requests, time, os, json, sys

from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)

BASE  = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
EMAIL = "sgould@wi.mit.edu"          # required by the EBI job dispatcher
OUT   = BLAST

CDKS = ["CDK7", "CDK8", "CDK19", "CDK9", "CDK12", "CDK13", "CDK2", "CDK4", "CDK6"]
# uniprotkb (metazoa-dense) already run in round 1
DBS = ["uniprotkb_swissprot", "uniprotkb_fungi", "uniprotkb_viridiplantae",
       "uniprotkb_archaea", "uniprotkb_bacteria"]

MAX_CONCURRENT = 12


def submit(gene, db):
    seq = open(f"{OUT}/{gene}_query.fasta").read()
    r = requests.post(f"{BASE}/run", timeout=120, data={
        "email": EMAIL, "program": "blastp", "matrix": "BLOSUM62",
        "alignments": 1000, "scores": 1000, "exp": "1e-5", "filter": "F",
        "gapalign": "true", "compstats": "F", "align": 0,
        "stype": "protein", "sequence": seq, "database": db})
    r.raise_for_status()
    return r.text.strip()


def main():
    todo = [(g, db) for g in CDKS for db in DBS
            if not os.path.exists(f"{OUT}/{g}__{db}_hits.tsv")]
    print(f"{len(todo)} searches to run", flush=True)
    running, done, failed = {}, [], []

    while todo or running:
        while todo and len(running) < MAX_CONCURRENT:
            g, db = todo.pop(0)
            try:
                running[(g, db)] = submit(g, db)
                print(f"  submit {g} / {db}", flush=True)
            except Exception as e:
                print(f"  SUBMIT FAIL {g}/{db}: {e}", flush=True)
                failed.append((g, db))
            time.sleep(2)

        time.sleep(25)
        for key, jid in list(running.items()):
            g, db = key
            try:
                st = requests.get(f"{BASE}/status/{jid}", timeout=60).text.strip()
            except Exception:
                continue
            if st == "FINISHED":
                tsv = requests.get(f"{BASE}/result/{jid}/tsv", timeout=300).text
                open(f"{OUT}/{g}__{db}_hits.tsv", "w").write(tsv)
                n = max(len(tsv.splitlines()) - 1, 0)
                print(f"DONE {g} / {db}: {n} hits  "
                      f"({len(done)+1}/{len(done)+len(running)+len(todo)})", flush=True)
                done.append(key); del running[key]
            elif st in ("ERROR", "FAILURE", "NOT_FOUND"):
                print(f"FAILED {g} / {db}: {st}", flush=True)
                failed.append(key); del running[key]

    print(f"ALL DONE  ok={len(done)} failed={len(failed)} {failed}", flush=True)


if __name__ == "__main__":
    main()
