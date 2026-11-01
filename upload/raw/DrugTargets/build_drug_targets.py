#!/usr/bin/env python3
# Inputs not included in this repository: Open Targets Platform 25.03 known_drug parquet dataset
# (253,459 drug-target-indication records) from https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.03/output/known_drug/ ;
# place the dataset folder in raw/DrugTargets/known_drug/.
# Also requires raw/GWASgene/study_index.tsv, built by raw/GWASgene/build_ot_inputs.py.
from pathlib import Path
import pyarrow.parquet as pq
import pandas as pd

U = set(open('data/onek1k/GECO_universe.txt').read().split())

panel = pd.read_csv('data/5b-ldsc/sensitivity/FinnGen_phenotypes_final.csv')
si = pd.read_csv('raw/GWASgene/study_index.tsv', sep='\t', usecols=['studyId', 'efo_id'])
panel = panel.merge(si, left_on='FinnGen ID', right_on='studyId')

kd = pq.read_table('raw/DrugTargets/known_drug', columns=['targetId', 'diseaseId', 'phase', 'prefName']).to_pandas()
kd = kd[kd.phase >= 4].merge(panel[['efo_id', 'Category']].drop_duplicates('efo_id'),
                             left_on='diseaseId', right_on='efo_id')
kd = (kd.rename(columns={'prefName': 'drug'})[['drug', 'targetId', 'diseaseId', 'Category']]
        .drop_duplicates().sort_values(['Category', 'drug', 'targetId', 'diseaseId']))
kd.to_csv('data/others/known_drug_approved_panel.tsv', sep='\t', index=False)

for cat in ['Immune', 'Non-Immune']:
    genes = sorted(set(kd.loc[kd.Category == cat, 'targetId']) & U)
    Path(f"data/5e-drug_target/targets_{cat.replace('-', '')}_approved.txt").write_text('\n'.join(genes) + '\n')

# CLUES: metformin excluded (METFORMIN / METFORMIN HYDROCHLORIDE records dropped), targets in the CLUES universe
Uc = set(open('data/clues/allgenes_CLUES.txt').read().split())
raw = pq.read_table('raw/DrugTargets/known_drug', columns=['targetId', 'diseaseId', 'phase', 'prefName']).to_pandas()
raw = raw[~raw.prefName.str.upper().str.contains('METFORMIN', na=False)]
kdc = raw.groupby(['targetId', 'diseaseId'], as_index=False).phase.max()
kdc = kdc[kdc.phase >= 4].merge(panel[['efo_id', 'Category']].drop_duplicates('efo_id'), left_on='diseaseId', right_on='efo_id')
for cat in ['Immune', 'Non-Immune']:
    genes = sorted(set(kdc.loc[kdc.Category == cat, 'targetId']) & Uc)
    Path(f"data/7d_drug_target/targets_{cat.replace('-', '')}_approved.txt").write_text('\n'.join(genes) + '\n')
