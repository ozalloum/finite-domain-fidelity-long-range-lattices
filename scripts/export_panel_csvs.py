"""Split each combined figure CSV into one CSV per plotted panel.

Most figures already carry a ``panel`` column. Figure 6 stores four heat-map
fields in one tidy table, so this script additionally exports one source CSV
for each of its four rendered panels.
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FD = ROOT / 'figure_data'

for p in sorted(FD.glob('figure[0-9][0-9].csv')):
    df = pd.read_csv(p)
    if 'panel' not in df.columns:
        continue
    for panel, g in df.groupby('panel', dropna=False):
        name = str(panel).strip().replace(' ', '_').replace('/', '_')
        out = FD / f'{p.stem}_panel_{name}.csv'
        g.to_csv(out, index=False)
        print(out.relative_to(ROOT))

# Figure 6: one CSV per heat-map panel, in addition to the combined tidy map.
f6 = FD / 'figure06.csv'
if f6.exists():
    m = pd.read_csv(f6)
    if {'Z', 'n', 'reference_intensity', 'open_intensity', 'corrected_intensity'}.issubset(m.columns):
        panel_specs = {
            'a_reference': ('reference_intensity', 'intensity'),
            'b_open': ('open_intensity', 'intensity'),
            'c_corrected': ('corrected_intensity', 'intensity'),
        }
        for name, (src, dst) in panel_specs.items():
            out = FD / f'figure06_panel_{name}.csv'
            m[['Z', 'n', src]].rename(columns={src: dst}).to_csv(out, index=False)
            print(out.relative_to(ROOT))
        d = m[['Z', 'n', 'open_intensity', 'reference_intensity']].copy()
        d['absolute_intensity_error'] = (d['open_intensity'] - d['reference_intensity']).abs()
        out = FD / 'figure06_panel_d_open_minus_reference.csv'
        d[['Z', 'n', 'absolute_intensity_error']].to_csv(out, index=False)
        print(out.relative_to(ROOT))
