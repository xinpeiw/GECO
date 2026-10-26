#!/usr/bin/env python3
# Inputs not included in this repository: Open Targets Platform 26.03 parquet datasets study/, credible_set/ and
# l2g_prediction/ in raw/GWASgene/otg_l2g_2603/ (see raw/GWASgene/build_ot_inputs.py).
# Output (not included, ~70 MB): raw/CLUES_GWASgene/gene_locus_L2G_allgenes.tsv.gz
from pathlib import Path
import pyarrow.parquet as pq
import pandas as pd

SRC = Path('raw/GWASgene/otg_l2g_2603')
OUT = Path('raw/CLUES_GWASgene')

sf = pq.read_table(SRC / 'study', filters=[('studyType', '=', 'gwas')]).to_pandas()
sf['efo_id'] = sf['traitFromSourceMappedIds'].apply(
    lambda a: a[0] if a is not None and len(a) else None)
study_idx = sf.loc[sf['efo_id'].notna(), ['studyId', 'efo_id']]

cs = pq.read_table(SRC / 'credible_set',
                   columns=['studyLocusId', 'studyId']).to_pandas().drop_duplicates()
cs = cs.merge(study_idx, on='studyId')

l2g = pq.read_table(SRC / 'l2g_prediction',
                    columns=['studyLocusId', 'geneId', 'score']).to_pandas()
l2g = l2g.rename(columns={'score': 'L2G_score'}).merge(cs, on='studyLocusId')
l2g = l2g.sort_values(['studyId', 'studyLocusId', 'geneId'])
l2g.to_csv(OUT / 'gene_locus_L2G_allgenes.tsv.gz', sep='\t',
           index=False, compression='gzip')
print(len(l2g), 'rows,', l2g.geneId.nunique(), 'genes,',
      l2g.efo_id.nunique(), 'disease terms')
