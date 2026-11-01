from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

OUT = Path('results/others'); OUT.mkdir(parents=True, exist_ok=True)

DATASETS = {
    'GECO': dict(universe=Path('data/onek1k/GECO_universe.txt'),
                 hubs={'Shared': Path('data/onek1k/GECO_Shared.txt'),
                       'Specific': Path('data/onek1k/GECO_Specific.txt')}),
    'CLUES mega': dict(universe=Path('data/clues/allgenes_CLUES.txt'),
                       hubs={'Shared': Path('data/clues/hubs_alpha.txt'),
                             'Specific': Path('data/clues/hubs_specific.txt')}),
}

kd = pd.read_csv('data/others/known_drug_approved_panel.tsv', sep='\t')


def fit(universe, hub, genes):
    y = np.isin(universe, list(genes)).astype(int)
    x = np.isin(universe, hub).astype(int)
    r = sm.GLM(y, sm.add_constant(x), family=sm.families.Binomial()).fit()
    return r.params[1] / np.log(2), r.pvalues[1]


rows = []
for ds, cfg in DATASETS.items():
    U = np.array(cfg['universe'].read_text().split())
    Uset = set(U)
    hubs = {name: path.read_text().split() for name, path in cfg['hubs'].items()}
    for cat in ['Immune', 'Non-Immune']:
        sub = kd[kd.Category == cat]
        full = set(sub.targetId) & Uset
        base = {h: fit(U, g, full) for h, g in hubs.items()}
        for drug in sub.drug.unique():
            keep = set(sub[sub.drug != drug].targetId) & Uset
            if keep == full:
                continue
            for h, g in hubs.items():
                lodo_or, lodo_p = fit(U, g, keep)
                rows.append(dict(dataset=ds, group=cat, hub=h, drug=drug,
                                 n_dropped=len(full) - len(keep),
                                 base_or=base[h][0], lodo_or=lodo_or,
                                 base_p=base[h][1], lodo_p=lodo_p))

res = pd.DataFrame(rows)
res.to_csv(OUT / 'DrugTargets_LODO.tsv', sep='\t', index=False)
print(f'{len(res)} refits -> {OUT / "DrugTargets_LODO.tsv"}')
for (ds, cat, h), g in res.groupby(['dataset', 'group', 'hub']):
    worst = g.loc[(g.lodo_or - g.base_or).abs().idxmax()]
    print(f'{ds:11s} {cat:11s} {h:9s} base={g.base_or.iloc[0]:+.2f} '
          f'LODO range [{g.lodo_or.min():+.2f}, {g.lodo_or.max():+.2f}] '
          f'| most influential: {worst.drug} ({int(worst.n_dropped)} targets) -> {worst.lodo_or:+.2f}')
