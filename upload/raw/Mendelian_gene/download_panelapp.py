#!/usr/bin/env python3
# Data source: Genomics England PanelApp public API, https://panelapp.genomicsengland.co.uk/api/v1 (accessed 2026-09-08).
# The API serves the current panel versions; the snapshot used in the paper is included as
# raw/Mendelian_gene/panelapp_green_gene_disease_pairs.csv (panel versions recorded per row).
# Outputs: raw/Mendelian_gene/panelapp_green_gene_disease_pairs.csv, raw/Mendelian_gene/panelapp_all_panels_metadata.json
import json, csv, subprocess, time

OUT = 'raw/Mendelian_gene'
API = 'https://panelapp.genomicsengland.co.uk/api/v1'
CATCHALL = {484, 478, 285}
RECLASS_IDS = {1350}
RECLASS_NAMES = {'Rare genetic inflammatory skin disorders',
                 'Intestinal failure or congenital diarrhoea'}

def get(url):
    for attempt in range(6):
        out = subprocess.run(['curl', '-s', '-A', 'Mozilla/5.0', '--max-time', '120', url],
                             capture_output=True, text=True).stdout
        try:
            return json.loads(out)
        except Exception:
            time.sleep(30 * (attempt + 1))
    raise RuntimeError(f'failed after retries: {url}')

panels, url = [], f'{API}/panels/?page=1'
while url:
    d = get(url)
    panels += d['results']
    url = d['next']
json.dump(panels, open(f'{OUT}/panelapp_all_panels_metadata.json', 'w'))
print('panels total:', len(panels))

gms = [p for p in panels if any(t['name'] == 'GMS signed-off' for t in p['types'])]
supers = {p['id'] for p in gms if any(t['name'] == 'Super Panel' for t in p['types'])}
kept = [p['id'] for p in gms if p['id'] not in supers and p['id'] not in CATCHALL]
print(f'GMS signed-off: {len(gms)} | super panels dropped: {len(supers)} | kept: {len(kept)}')

rows = []
for pid in sorted(kept):
    d = get(f'{API}/panels/{pid}/')
    grp = d.get('disease_group') or 'NA'
    immune = int(grp == 'Immunology' or pid in RECLASS_IDS or d['name'] in RECLASS_NAMES)
    for g in d['genes']:
        if g['confidence_level'] != '3':
            continue
        gd = g['gene_data']
        rows.append([pid, d['name'], d['version'], grp, immune, gd['gene_symbol'],
                     gd['ensembl_genes']['GRch38']['90']['ensembl_id'],
                     g.get('mode_of_inheritance') or 'NA'])

with open(f'{OUT}/panelapp_green_gene_disease_pairs.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['panel_id', 'panel_name', 'panel_version', 'disease_group', 'immune',
                'gene_symbol', 'ensembl_id', 'mode_of_inheritance'])
    w.writerows(rows)

imm = {r[6] for r in rows if r[4] == 1}
non = {r[6] for r in rows if r[4] == 0}
print(f'green pairs {len(rows)} | immune {len(imm)} | non-immune {len(non)} | overlap {len(imm & non)}')
