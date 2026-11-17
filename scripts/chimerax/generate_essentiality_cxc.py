#!/usr/bin/env python3
"""
ChimeraX panels of the early (T0 vs plasmid) residue essentiality score.

Paints the score onto the cartoon of each CDK, in the same orientation, with the
same colour scale and the same limits as the linear tracks in
`figures/T0_plasmid_essentiality/essentiality_tracks.pdf`, so a structural panel
and its track can be placed side by side.

Reads  screening_data/03_LFC_FDR_tables/T0_plasmid_essentiality/residue_essentiality.csv
Writes ./cxc_essentiality/
         attr/<GENE>.defattr             per-residue attribute files
         <GENE>.cxc                      one panel per CDK (cyclin surface shown)
         <GENE>_nosurface.cxc            same, cyclin hidden
         00_all_cdks_essentiality.cxc    all nine superposed on CDK2
         render_all.cxc                  renders every panel to figures/chimerax/essentiality

Usage:  python3 generate_essentiality_cxc.py [--score ess_early] [--out cxc_essentiality]
        then, from the ChimeraX GUI (it needs a real GPU context):
            open cxc_essentiality/render_all.cxc
        or   ChimeraX --exit --cmd "open <abs path>/cxc_essentiality/render_all.cxc"
"""
import argparse, csv, gzip, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
from cdk_config import (STRUCTURES, ALIGN_REFERENCE, CDK_COLORS, CYCLIN_COLOR, CYCLIN_TRANSPARENCY,
                        OTHER_CHAIN_COLOR, RENDER_PREAMBLE, CANONICAL_VIEW, ZOOM_OUT,
                        WINDOW_SIZE)

sys.path.insert(0, os.path.join(ROOT, 'src'))
from cdkbe import paths as _paths  # noqa: E402
from cdkbe import io as _io        # noqa: E402

# Was screening_data/03_LFC_FDR_tables/... in the original tree; `io.resolve`
# accepts the table whether or not it was gzipped on the way into data/.
SCORES = _io.resolve(_paths.LFC / 'T0_plasmid_essentiality' / 'residue_essentiality.csv')
DEPLETED_CSV = os.path.join(HERE, 'essentiality_depleted_residues.csv')
FIGDIR = os.path.join(ROOT, 'figures/chimerax/essentiality')

# identical to essentiality_figures.py
DEPLETED, MID, ENRICHED = '#2166ac', '#f2f2f2', '#b2182b'
MID_LO, MID_HI = '#d1e5f0', '#fddbc7'
NO_DATA = '#9e9e9e'
VLIM = 1.0        # ~1.6 x the permutation-null sd of a residue score:
                  # fully saturated colour is about where FDR < 0.05 starts
PALETTE = (f'-{VLIM},{DEPLETED}:-{VLIM/2},{MID_LO}:0,{MID}:'
           f'{VLIM/2},{MID_HI}:{VLIM},{ENRICHED}')

SIG_FDR = 0.05
ORDER = [ALIGN_REFERENCE] + [g for g in STRUCTURES if g != ALIGN_REFERENCE]

_CHAINS = {"CDK7": list("HIJ"), "CDK8": list("AB"), "CDK19": ["A"], "CDK9": list("AB"),
           "CDK12": list("ABCDEFGHI"), "CDK13": list("ABCD"), "CDK2": list("ABCDEF"),
           "CDK4": list("AB"), "CDK6": list("AB")}


def load_scores(score_col):
    # SCORES may be the gzipped spelling (see io.resolve above), so pick the opener
    # from the suffix rather than assuming plain text.
    _open = gzip.open if str(SCORES).endswith('.gz') else open
    with _open(SCORES, 'rt') as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for r in rows:
        v = r.get(score_col, '')
        if v in ('', 'nan', 'NaN'):
            continue
        fdr = r.get('FDR_early', '')
        out.setdefault(r['gene' if 'gene' in r else 'Gene'], []).append(
            (int(r['Codon']), float(v),
             float(fdr) if fdr not in ('', 'nan', 'NaN') else None))
    return out


ATTR_REF = None       # path to the attr/ directory as the .cxc files must see it


def write_defattr(out, gene, mid, vals, score_col):
    cfg = STRUCTURES[gene]
    kc = cfg['kinase_chain']
    dr = cfg.get('display_range')
    path = os.path.join(out, 'attr', f'{gene}.defattr')
    with open(path, 'w') as fh:
        fh.write(f'# early essentiality score ({score_col}) for {gene}\n')
        fh.write('# written by generate_essentiality_cxc.py\n')
        fh.write('attribute: essentiality\n')
        fh.write('match mode: any\n')
        fh.write('recipient: residues\n')
        n = 0
        for pos, v, _ in sorted(vals):
            if dr and not (dr[0] <= pos <= dr[1]):
                continue
            fh.write(f'\t#{mid}/{kc}:{pos}\t{v:.4f}\n')
            n += 1
    return path, n


def model_block(mid, gene, vals, lines, out, score_col, show_surface=True):
    cfg = STRUCTURES[gene]
    kc, cyc, other = cfg['kinase_chain'], cfg['cyclin_chains'], cfg['other_chains']
    dr = cfg.get('display_range')

    _, n = write_defattr(out, gene, mid, vals, score_col)

    lines += [f'\n# ---------------------------------------------------------------- {gene}',
              f"# {cfg['pdb']}: {' '.join(cfg['note'].split())}",
              f'open structures/{cfg["file"]} name {gene}',
              f'hide #{mid} target acsb',
              f'hide #{mid} models']

    if dr:
        lines.append(f'# only residues {dr[0]}-{dr[1]} are confidently modelled')
        lines.append(f'show #{mid}/{kc}:{dr[0]}-{dr[1]} cartoon')
    else:
        lines.append(f'show #{mid}/{kc} cartoon')

    lines += [f'# {n} residues carry a score; everything else stays {NO_DATA}',
              f'open {ATTR_REF}/{gene}.defattr',
              f'color byattribute r:essentiality #{mid}/{kc} target c '
              f'palette {PALETTE} noValueColor {NO_DATA}']

    if cyc and show_surface:
        cs = f'#{mid}/' + ','.join(cyc)
        lines += [f"# {cfg['cyclin_name']} as a wheat surface, the convention used by",
                  '# generate_cxc.py for every other structural figure in this repo',
                  f'show {cs} atoms', f'surface {cs}', f'hide {cs} atoms',
                  f'color {cs} {CYCLIN_COLOR} target s']
        if CYCLIN_TRANSPARENCY:
            lines.append(f'transparency {cs} {CYCLIN_TRANSPARENCY} target s')
    elif cyc:
        lines.append(f'# cyclin chain(s) {",".join(cyc)} left hidden')

    if cfg['ligand']:
        ls = f'#{mid}/{kc} & :{cfg["ligand"]}'
        lines += [f'show {ls} target ab', f'style {ls} stick',
                  f'color {ls} #4d4d4d target ab', f'color {ls} byhetero']
        if gene == 'CDK2':
            lines += [f'show #{mid}/{kc} & :MG atoms', f'style #{mid}/{kc} & :MG sphere',
                      f'size #{mid}/{kc} & :MG atomRadius 0.8']

    # significance is one-sided: FDR < 0.05 *and* depleted. The two-sided list is not
    # usable at the default level, where most of its calls are GC-driven enrichments.
    sig = sorted({p for p, v, q in vals if q is not None and q < SIG_FDR and v < 0
                  and (not dr or dr[0] <= p <= dr[1])})
    if sig:
        spec = f'#{mid}/{kc}:' + ','.join(str(p) for p in sig)
        lines += [f'# {len(sig)} significantly DEPLETED residues (FDR < {SIG_FDR}, '
                  'score < 0): side chains as sticks,',
                  '# still coloured by the score so the figure has one colour meaning',
                  '# target `ab` matters: the blanket `hide ... target acsb` above turns bond',
                  '# display off, and `show ... atoms` alone would draw isolated dots with',
                  '# no sticks between them (sticks are bonds)',
                  f'show {spec} target ab', f'cartoon {spec} suppressBackboneDisplay false',
                  f'style {spec} stick', f'size {spec} stickRadius 0.28',
                  f'color byattribute r:essentiality {spec} target a '
                  f'palette {PALETTE} noValueColor {NO_DATA}',
                  f'name {gene}_significant {spec}']

    if other:
        lines.append(f'# hidden partners in this file: {", ".join(other)}')
        lines.append(f'color #{mid}/' + ','.join(other) + f' {OTHER_CHAIN_COLOR}')
    unused = [c for c in _CHAINS[gene] if c not in [kc] + cyc + other]
    if unused:
        lines.append(f'# extra copies in the asymmetric unit, left hidden: {", ".join(unused)}')
    lines.append(f'name {gene}_kinase #{mid}/{kc}')


def borrowed_block(mid, gene, host_mid, lines):
    """Superpose the donor CDK onto `gene` and show the donor's cyclin.

    Ported from generate_cxc.py: CDK19 is the CDK8 paralogue and shares cyclin C, but no
    CDK19-cyclin C structure exists, so cyclin C is borrowed from 5CEI and drawn more
    transparent to make clear it is inferred rather than observed.
    """
    cfg = STRUCTURES[gene]
    b = cfg['borrowed_cyclin']
    donor = STRUCTURES[b['donor']]
    lines += [f'\n# ------------------------------- borrowed cyclin for {gene}',
              f"# {gene} has no structure with a cyclin bound; {b['donor']} "
              f"({donor['pdb']}) is superposed",
              f"# onto it and its {donor['cyclin_name']} shown at "
              f"{b['transparency']}% transparency.",
              f"open structures/{donor['file']} name "
              f"{gene}_borrowed_{donor['cyclin_name'].replace(' ', '')}",
              f'hide #{mid} target acsb',
              f"matchmaker #{mid}/{donor['kinase_chain']} to "
              f"#{host_mid}/{cfg['kinase_chain']}"]
    cs = f'#{mid}/' + ','.join(donor['cyclin_chains'])
    lines += [f'show {cs} atoms', f'surface {cs}', f'hide {cs} atoms',
              f'color {cs} {CYCLIN_COLOR} target s',
              f"transparency {cs} {b['transparency']} target s",
              f'name {gene}_borrowed_cyclin {cs}']


def footer(lines):
    lines.append('\n# ---------------------------------------------------- look and camera')
    lines += RENDER_PREAMBLE.strip().splitlines()
    lines += ['show #* models', 'label delete', 'hide #* target p',
              'hide @@element=H target ab',
              f'view matrix camera {CANONICAL_VIEW}',
              f'zoom {ZOOM_OUT}', '']


def header(title, extra=''):
    return ['# ' + '=' * 74, f'# {title}',
            '# Generated by generate_essentiality_cxc.py -- edit that script, not this file.',
            f'# colour: {PALETTE}   (no value: {NO_DATA})',
            f'# {extra}' if extra else '# ', '# ' + '=' * 74,
            'close session', f'cd {HERE}', '']


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--score', default='ess_early',
                    help='column of residue_essentiality.csv to paint (default ess_early)')
    ap.add_argument('--out', default=os.path.join(HERE, 'cxc_essentiality'))
    ap.add_argument('--figdir', default=FIGDIR)
    a = ap.parse_args()
    global ATTR_REF
    out = os.path.abspath(a.out)
    os.makedirs(os.path.join(out, 'attr'), exist_ok=True)
    # the .cxc files `cd` to this script's directory first, so the attribute files
    # must be referenced relative to that, not to the .cxc file itself
    ATTR_REF = os.path.join(os.path.relpath(out, HERE), 'attr')
    os.makedirs(a.figdir, exist_ok=True)

    scores = load_scores(a.score)
    print(f'scores: {SCORES}  ({a.score})')
    for g in ORDER:
        print(f'  {g:<6} {len(scores.get(g, [])):>5} residues scored')

    # ------------------------------------------------------------- combined
    lines = header('CDK-family residue essentiality (early dropout), superposed on '
                   f'{ALIGN_REFERENCE}',
                   'Model numbering: ' + ', '.join(f'#{i+1}={g}' for i, g in enumerate(ORDER)))
    for i, gene in enumerate(ORDER):
        model_block(i + 1, gene, scores.get(gene, []), lines, out, a.score,
                    show_surface=False)
    lines.append('\n# ------------------------------------------------ structural alignment')
    ref = f'#1/{STRUCTURES[ALIGN_REFERENCE]["kinase_chain"]}'
    for i, gene in enumerate(ORDER):
        if gene != ALIGN_REFERENCE:
            lines.append(f'matchmaker #{i+1}/{STRUCTURES[gene]["kinase_chain"]} to {ref}')
    footer(lines)
    open(os.path.join(out, '00_all_cdks_essentiality.cxc'), 'w').write('\n'.join(lines) + '\n')

    # --------------------------------------------------------- single panels
    panels = ['00_all_cdks_essentiality']
    for gene in ORDER:
        for suffix, surf in [('', True), ('_nosurface', False)]:
            cfg = STRUCTURES[gene]
            lines = header(f'{gene} residue essentiality ({cfg["pdb"]})',
                           ' '.join(cfg['note'].split()))
            if gene == ALIGN_REFERENCE:
                model_block(1, gene, scores.get(gene, []), lines, out, a.score, surf)
            else:
                # the reference is loaded and hidden, so every panel shares one frame
                lines += [f'# {ALIGN_REFERENCE} is loaded only to put this panel in the '
                          'same frame as the others',
                          f'open structures/{STRUCTURES[ALIGN_REFERENCE]["file"]} '
                          f'name {ALIGN_REFERENCE}_reference',
                          'hide #1 target acsb', 'hide #1 models']
                model_block(2, gene, scores.get(gene, []), lines, out, a.score, surf)
                lines.append(f'matchmaker #2/{cfg["kinase_chain"]} to '
                             f'#1/{STRUCTURES[ALIGN_REFERENCE]["kinase_chain"]}')
                lines.append('hide #1 models')
                if surf and 'borrowed_cyclin' in cfg:
                    borrowed_block(3, gene, 2, lines)
            footer(lines)
            if gene != ALIGN_REFERENCE:
                lines.append('# keep the alignment reference out of the render')
                lines.append('hide #1 models')
            name = f'{gene}{suffix}'
            open(os.path.join(out, f'{name}.cxc'), 'w').write('\n'.join(lines) + '\n')
            panels.append(name)

    # ------------------------------------- significantly depleted residues as atoms
    PALE = '#D9D9D9'

    def depleted_block(mid, gene, vals, lines, show_surface, per_cdk_colour):
        """Pale cartoon, with every significantly depleted residue drawn as spheres."""
        cfg = STRUCTURES[gene]
        kc, cyc = cfg['kinase_chain'], cfg['cyclin_chains']
        dr = cfg.get('display_range')
        lines += [f'\n# ---------------------------------------------------------- {gene}',
                  f"# {cfg['pdb']}: {' '.join(cfg['note'].split())}",
                  f'open structures/{cfg["file"]} name {gene}',
                  f'hide #{mid} target acsb', f'hide #{mid} models']
        if dr:
            lines.append(f'show #{mid}/{kc}:{dr[0]}-{dr[1]} cartoon')
        else:
            lines.append(f'show #{mid}/{kc} cartoon')
        lines.append(f'color #{mid}/{kc} {PALE} target c')
        if cyc and show_surface:
            cs = f'#{mid}/' + ','.join(cyc)
            lines += [f'show {cs} atoms', f'surface {cs}', f'hide {cs} atoms',
                      f'color {cs} {CYCLIN_COLOR} target s',
                      f'transparency {cs} 55 target s']
        if cfg['ligand']:
            ls = f'#{mid}/{kc} & :{cfg["ligand"]}'
            lines += [f'show {ls} target ab', f'style {ls} stick',
                      f'color {ls} #707070 target ab', f'color {ls} byhetero']
        dep = sorted({p for p, v, q in vals
                      if q is not None and q < SIG_FDR and v < 0
                      and (not dr or dr[0] <= p <= dr[1])})
        if dep:
            spec = f'#{mid}/{kc}:' + ','.join(str(p) for p in dep)
            lines += [f'# {len(dep)} significantly depleted residues, drawn as atoms',
                      f'show {spec} target ab',
                      f'cartoon {spec} suppressBackboneDisplay false',
                      f'style {spec} sphere', f'size {spec} atomRadius 0.75']
            if per_cdk_colour:
                lines.append(f'color {spec} {CDK_COLORS[gene]} target ac')
            else:
                lines.append(f'color byattribute r:essentiality {spec} target ac '
                             f'palette {PALETTE} noValueColor {NO_DATA}')
                lines.insert(len(lines) - 1, f'open {ATTR_REF}/{gene}.defattr')
            lines.append(f'name {gene}_depleted {spec}')
        return len(dep)

    for gene in ORDER:
        cfg = STRUCTURES[gene]
        lines = header(f'{gene} — significantly depleted residues ({cfg["pdb"]})',
                       'pale cartoon; depleted residues (FDR < 0.05, score < 0) as spheres '
                       'coloured by score')
        if gene == ALIGN_REFERENCE:
            n = depleted_block(1, gene, scores.get(gene, []), lines, True, False)
        else:
            lines += [f'open structures/{STRUCTURES[ALIGN_REFERENCE]["file"]} '
                      f'name {ALIGN_REFERENCE}_reference',
                      'hide #1 target acsb', 'hide #1 models']
            n = depleted_block(2, gene, scores.get(gene, []), lines, True, False)
            lines.append(f'matchmaker #2/{cfg["kinase_chain"]} to '
                         f'#1/{STRUCTURES[ALIGN_REFERENCE]["kinase_chain"]}')
            if 'borrowed_cyclin' in cfg:
                borrowed_block(3, gene, 2, lines)
        footer(lines)
        if gene != ALIGN_REFERENCE:
            lines += ['# keep the alignment reference out of the render', 'hide #1 models']
        open(os.path.join(out, f'{gene}_depleted.cxc'), 'w').write('\n'.join(lines) + '\n')
        panels.append(f'{gene}_depleted')
        print(f'  {gene}_depleted: {n} residues')

    # all nine at once, superposed, each CDK's depleted residues in its own colour
    lines = header('All nine CDKs superposed — significantly depleted residues',
                   'pale cartoons; spheres coloured by CDK')
    for i, gene in enumerate(ORDER):
        depleted_block(i + 1, gene, scores.get(gene, []), lines, False, True)
    lines.append('\n# ----------------------------------------------- structural alignment')
    ref = f'#1/{STRUCTURES[ALIGN_REFERENCE]["kinase_chain"]}'
    for i, gene in enumerate(ORDER):
        if gene != ALIGN_REFERENCE:
            lines.append(f'matchmaker #{i+1}/{STRUCTURES[gene]["kinase_chain"]} to {ref}')
    footer(lines)
    open(os.path.join(out, '02_all_cdks_depleted.cxc'), 'w').write('\n'.join(lines) + '\n')
    panels.append('02_all_cdks_depleted')

    # ------------------------------------------------- family-consensus peaks
    peak_csv = os.path.join(HERE, 'essentiality_peak_residues.csv')
    if os.path.exists(peak_csv):
        prs = list(csv.DictReader(open(peak_csv)))
        if prs:
            ref = ALIGN_REFERENCE
            cfg = STRUCTURES[ref]
            kc = cfg['kinase_chain']
            lines = header(f'Family-recurrent essentiality peaks on {ref} ({cfg["pdb"]})',
                           'from essentiality_peaks.py; a peak is a window of the CDK '
                           'alignment where depleted residues recur across the family')
            model_block(1, ref, scores.get(ref, []), lines, out, a.score,
                        show_surface=False)
            by_rank = {}
            for r in prs:
                by_rank.setdefault(r['peak_rank'], []).append(int(r['position']))
            for rank, pos in sorted(by_rank.items()):
                spec = f'#1/{kc}:' + ','.join(str(p) for p in sorted(pos))
                klifs = sorted({r['klifs'] for r in prs
                                if r['peak_rank'] == rank and r['klifs']})
                lines += [f'\n# peak {rank}: {len(pos)} residues'
                          + (f'  ({", ".join(klifs)})' if klifs else ''),
                          f'show {spec} target ab',
                          f'cartoon {spec} suppressBackboneDisplay false',
                          f'style {spec} sphere', f'size {spec} atomRadius 0.9',
                          f'color {spec} {DEPLETED} target ac',
                          f'name peak{rank} {spec}']
            footer(lines)
            open(os.path.join(out, 'peaks_consensus.cxc'), 'w').write('\n'.join(lines) + '\n')
            panels.append('peaks_consensus')
            print(f'  peaks_consensus: {sum(len(v) for v in by_rank.values())} '
                  f'{ref} residues in {len(by_rank)} peak(s)')

    # -------------------------------------------------------------- renderer
    r = ['# ' + '=' * 74,
         '# Render every essentiality panel at publication size.',
         '# ChimeraX needs a real GPU context for `save`, so run this from the GUI',
         f'#   ChimeraX --exit --cmd "open {out}/render_all.cxc"',
         '# ' + '=' * 74,
         'close session', f'cd {HERE}',
         f'windowsize {WINDOW_SIZE[0]} {WINDOW_SIZE[1]}', '']
    rel = os.path.relpath(out, HERE)
    for name in panels:
        r += [f'open {rel}/{name}.cxc',
              f'save {a.figdir}/{name}.png width 3000 supersample 3 '
              f'transparentBackground true', '']
    r.append('close session')
    open(os.path.join(out, 'render_all.cxc'), 'w').write('\n'.join(r) + '\n')
    print(f'\nwrote {len(panels)} panels to {out}')
    print(f'render with:  ChimeraX --exit --cmd "open {out}/render_all.cxc"')


if __name__ == '__main__':
    main()
