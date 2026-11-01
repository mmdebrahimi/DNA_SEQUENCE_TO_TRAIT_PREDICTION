"""De-confound the in-distribution-length question by holding MASK DENSITY fixed.

Why this exists: the first probe held mask COUNT fixed, so density = count/length varied 70% -> 20%
across the lengths tested. At 20 masks the accuracy curve peaked at 150 bp; at the true pretraining
count of 60 it rose monotonically to 300 bp. A ranking that inverts when a nuisance parameter
changes is measuring the nuisance parameter. Density is the confound, so pin it.

If accuracy is FLAT across lengths at fixed density, the model is genuinely length-robust and the
arm's 150->200 geometry is safe. If it PEAKS, that peak is the in-distribution window.
"""
import os
#: Where the public MIT PromoterAtlas clone lives. Override with PA_REPO.
_PA_DEFAULT = r'C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas'
PA_REPO = os.environ.get('PA_REPO', _PA_DEFAULT)
import importlib.util, json, pathlib, random
import torch

REPO = pathlib.Path(PA_REPO)
NUC = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

spec = importlib.util.spec_from_file_location(
    'pa_dt', REPO / 'src/promoter_atlas/models/dna_transformer.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
model = mod.DNATransformer(dim=128, n_blocks=8, n_heads=8)

pathlib.PosixPath = pathlib.WindowsPath
ck = torch.load(REPO / 'trained_weights/base_model/promoteratlas-base.pt',
                map_location='cpu', weights_only=False)   # justified by pickle_audit.py
model.load_state_dict(ck['model_state_dict'], strict=True)
model.eval()

txt = (REPO / 'data/genomes/NC_000913.1.gb').read_text(encoding='utf-8', errors='replace')
genome = ''.join(''.join(c for c in ln if c.isalpha()).upper()
                 for ln in txt.split('ORIGIN', 1)[-1].splitlines())
genome = ''.join(c if c in 'ACGT' else 'A' for c in genome)
print('K-12 parsed: %d bp' % len(genome))


def onehot(s):
    t = torch.zeros((4, len(s)))
    for i, c in enumerate(s):
        t[NUC[c], i] = 1.0
    return t


@torch.no_grad()
def acc_at(L, density, n=32, seed=11, shuffle=False):
    rng = random.Random(seed)
    mrng = random.Random(seed + 7919)
    seqs = []
    for _ in range(n):
        p = rng.randrange(0, len(genome) - L)
        s = genome[p:p + L]
        if shuffle:
            ls = list(s); mrng.shuffle(ls); s = ''.join(ls)
        seqs.append(s)
    X = torch.stack([onehot(s) for s in seqs])
    truth = X.argmax(dim=1)
    k = max(1, int(round(density * L)))
    pos = [mrng.sample(range(L), k) for _ in range(n)]
    Xm = X.clone()
    for b, ps in enumerate(pos):
        for p in ps:
            Xm[b, :, p] = 0.0
    logits, _ = model(Xm)
    c = t = 0
    for b, ps in enumerate(pos):
        for p in ps:
            c += int(logits[b, :, p].argmax().item() == truth[b, p].item())
            t += 1
    return c / t, k


LENGTHS = (86, 120, 150, 200, 250, 300, 400)
out = {}
for density in (0.10, 0.20, 0.30):
    print('\n--- mask density %.0f%% (held CONSTANT across lengths) ---' % (100 * density))
    print('   %-7s %-7s %-11s %-11s %s' % ('len', 'n_mask', 'real', 'shuffled', 'real-shuffled'))
    rows = {}
    for L in LENGTHS:
        a, k = acc_at(L, density)
        s, _ = acc_at(L, density, shuffle=True)
        rows[L] = {'n_mask': k, 'real': round(a, 4), 'shuffled': round(s, 4),
                   'delta': round(a - s, 4)}
        print('   %-7d %-7d %-11.4f %-11.4f %+.4f' % (L, k, a, s, a - s))
    best = max(rows, key=lambda L: rows[L]['real'])
    spread = max(r['real'] for r in rows.values()) - min(r['real'] for r in rows.values())
    rows['_best'] = best
    rows['_spread'] = round(spread, 4)
    out['density_%.2f' % density] = rows
    print('   best length: %d    spread across lengths: %.4f' % (best, spread))

p = pathlib.Path(os.environ.get('PA_OUT', 'wiki/promoteratlas_density_probe_2026-11-01.json'))
p.write_text(json.dumps(out, indent=1), encoding='utf-8')

bests = [out[k]['_best'] for k in out]
spreads = [out[k]['_spread'] for k in out]
print('\n--- verdict ---')
print('best length per density setting:', bests)
print('max spread across lengths:      %.4f' % max(spreads))
if len(set(bests)) == 1:
    print('STABLE: the same length wins at every density -> %d is the in-distribution window.' % bests[0])
elif max(spreads) < 0.05:
    print('FLAT: accuracy varies < 0.05 across 86-400 bp at every density. The model is genuinely')
    print('length-robust in this range, so the arm\'s 150->200 geometry is SAFE, and reconstruction')
    print('cannot (and need not) identify a single in-distribution window.')
else:
    print('UNSTABLE: the winner moves with density and the spread is large. Reconstruction accuracy')
    print('does NOT identify an in-distribution length. The bar must treat length as UNRESOLVED.')
print('\nwrote', p)
