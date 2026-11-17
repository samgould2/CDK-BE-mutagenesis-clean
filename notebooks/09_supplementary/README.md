# Supplementary figures

One notebook per supplementary figure, generated from
`docs/source_data_mapping_supp.xlsx`:

```bash
python3 scripts/make_supplementary_notebooks.py
```

Each notebook carries its panel table — source file, whether that file is in
the repository, and your notes — and nothing else yet. **No panel is ported.**

Current state across the 16 figures (103 panels):

| | panels |
|---|---:|
| source file named and present | 54 |
| awaited (marked `Missing`) | 29 |
| no underlying data (blots, Sanger traces) | 8 |
| unassigned | 12 |

Regenerate the notebooks after editing the mapping; they are derived files, so
edit the spreadsheet rather than the notebook.

To port a panel, find the notebook or script that wrote its source file, move
the reduction into `cdkbe` and the drawing into `cdkbe.source_data`, then add a
cell that calls it. Most layouts already exist — `cdkbe.lollipop`,
`cdkbe.heatmap` and `cdkbe.logo_panel` cover the lollipops, the clustered
heatmaps and the bar-plus-logo panels respectively.
