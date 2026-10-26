#!/usr/bin/env python3
# Inputs not included in this repository: Open Targets Platform 26.03 parquet datasets study/, disease/, credible_set/
# and l2g_prediction/ from https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/ ; place the four
# dataset folders in raw/GWASgene/otg_l2g_2603/.
# Outputs (not included, ~100 MB): raw/GWASgene/study_index.tsv and raw/GWASgene/gene_locus_features.tsv.gz
from pathlib import Path
import pyarrow.parquet as pq
import pandas as pd, numpy as np

SRC  = Path('raw/GWASgene/otg_l2g_2603')
OUT  = Path('raw/GWASgene')
geco = set(open('data/onek1k/GECO_universe.txt').read().splitlines())

sf = pq.read_table(SRC / 'study', filters=[('studyType', '=', 'gwas')]).to_pandas()
sf['efo_id'] = sf['traitFromSourceMappedIds'].apply(lambda a: a[0] if a is not None and len(a) else None)
study_idx = sf.loc[sf['efo_id'].notna(),
                   ['studyId', 'projectId', 'traitFromSource', 'efo_id', 'nCases', 'nControls', 'nSamples']].copy()

dz = pq.read_table(SRC / 'disease').to_pandas()
ta = dz[dz['ontology'].apply(lambda x: bool(x is not None and x.get('isTherapeuticArea', False)))]
ta_label = dict(zip(ta['id'],
                    ta['name'].str.replace(r'[^A-Za-z0-9]+', '_', regex=True).str.strip('_').str.lower()))
def ta_of(ids):
    labels = [ta_label[t] for t in ids if t in ta_label] if ids is not None else []
    return labels
efo_tas = dict(zip(dz['id'], dz['therapeuticAreas'].apply(ta_of)))
study_idx['ta'] = study_idx['efo_id'].map(efo_tas).apply(lambda x: x if isinstance(x, list) else [])
study_idx.to_csv(OUT / 'study_index.tsv', sep='\t', index=False)

cs = pq.read_table(SRC / 'credible_set', columns=['studyLocusId', 'studyId']).to_pandas().drop_duplicates()
cs = cs.merge(study_idx[['studyId', 'efo_id']].drop_duplicates(), on='studyId')

l2g = pq.read_table(SRC / 'l2g_prediction',
                    columns=['studyLocusId', 'geneId', 'score', 'features']).to_pandas()
l2g = l2g.merge(cs, on='studyLocusId')
l2g = l2g[l2g['geneId'].isin(geco)].copy()

CHANNELS = {'L2G_score': None, 'eQtl_coloc': 'eQtlColocH4Maximum', 'sQtl_coloc': 'sQtlColocH4Maximum',
            'pQtl_coloc': 'pQtlColocH4Maximum', 'enhancer_to_gene': 'e2gMean', 'nearest_tss': 'distanceTssMean'}
feat_names = [v for v in CHANNELS.values() if v is not None]
def pull(feats):
    out = {n: np.nan for n in feat_names}
    if feats is None: return out
    for f in feats:
        if f['name'] in out: out[f['name']] = f['value']
    return out
ext = pd.DataFrame(list(l2g['features'].apply(pull)))
for ch, fn in CHANNELS.items():
    l2g[ch] = l2g['score'] if fn is None else ext[fn].values
l2g = l2g.drop(columns='features').sort_values(['studyId', 'studyLocusId', 'geneId'])
l2g.to_csv(OUT / 'gene_locus_features.tsv.gz', sep='\t', index=False, compression='gzip')
