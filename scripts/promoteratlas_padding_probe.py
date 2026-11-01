"""Does zero-padding a 150 bp tile to 200 bp change its frozen embedding?

This is a design point the arm's own NAME pre-judges. `...-frozen-150bp-padded-to-200` was declared
before anyone had established that the base model is length-agnostic. It is, so padding is now
OPTIONAL -- and zero-padding is not neutral here: an all-zero one-hot column is exactly how
pretraining represented a MASKED base (mask multiplies the one-hot). Padding to 200 therefore asks
the model to fill 50 masked positions, which can leak into the pooled latent.

Measured three ways:
  1. masked mean-pool at native 150 vs masked mean-pool of the 150->200 padded input
     (masked pooling excludes padded positions, so a difference is CONTAMINATION, not geometry)
  2. unmasked mean-pool of the padded input (what a careless implementation would do)
  3. whether the 128-d embedding is non-degenerate at all (per-dim variance across tiles)
"""
import os
#: Where the public MIT PromoterAtlas clone lives. Override with PA_REPO.
_PA_DEFAULT = r'C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas'
PA_REPO = os.environ.get('PA_REPO', _PA_DEFAULT)
import importlib.util, json, pathlib, random
import torch

REPO = pathlib.Path(PA_REPO)
NUC = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
L_NAT, L_PAD, N = 150, 200, 256

spec = importlib.util.spec_from_file_location(
    'pa_dt', REPO / 'src/promoter_atlas/models/dna_transformer.py')
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
model = mod.DNATransformer(dim=128, n_blocks=8, n_heads=8)
pathlib.PosixPath = pathlib.WindowsPath
ck = torch.load(REPO / 'trained_weights/base_model/promoteratlas-base.pt',
                map_location='cpu', weights_only=False)
model.load_state_dict(ck['model_state_dict'], strict=True); model.eval()

txt = (REPO / 'data/genomes/NC_000913.1.gb').read_text(encoding='utf-8', errors='replace')
genome = ''.join(''.join(c for c in ln if c.isalpha()).upper()
                 for ln in txt.split('ORIGIN', 1)[-1].splitlines())
genome = ''.join(c if c in 'ACGT' else 'A' for c in genome)

rng = random.Random(4242)
tiles = []
for _ in range(N):
    p = rng.randrange(0, len(genome) - L_NAT)
    tiles.append(genome[p:p + L_NAT])


def onehot_batch(seqs, width):
    X = torch.zeros((len(seqs), 4, width))
    for b, s in enumerate(seqs):
        for i, c in enumerate(s):
            X[b, NUC[c], i] = 1.0
    return X


@torch.no_grad()
def latents(X):
    _, _, z = model(X, return_latent=True)
    return z                                        # (B, 128, L)


Xn = onehot_batch(tiles, L_NAT)
Xp = onehot_batch(tiles, L_PAD)                     # positions 150..199 are all-zero = MASKED
zn, zp = latents(Xn), latents(Xp)

emb_native = zn.mean(dim=2)                                  # (B,128) masked-pool == plain here
emb_pad_masked = zp[:, :, :L_NAT].mean(dim=2)                # padded positions EXCLUDED
emb_pad_plain = zp.mean(dim=2)                               # padded positions INCLUDED

out = {}


def compare(a, b, label):
    d = (a - b).abs()
    cos = torch.nn.functional.cosine_similarity(a, b, dim=1)
    # per-tile correlation across the 128 dims
    am = a - a.mean(dim=1, keepdim=True); bm = b - b.mean(dim=1, keepdim=True)
    r = (am * bm).sum(1) / (am.norm(dim=1) * bm.norm(dim=1) + 1e-12)
    rec = {'max_abs_delta': round(d.max().item(), 6),
           'mean_abs_delta': round(d.mean().item(), 6),
           'median_cosine': round(cos.median().item(), 6),
           'min_cosine': round(cos.min().item(), 6),
           'median_pearson_over_dims': round(r.median().item(), 6)}
    out[label] = rec
    print('   %-34s max|d| %-10.5f mean|d| %-10.5f cos med %-9.6f min %-9.6f'
          % (label, rec['max_abs_delta'], rec['mean_abs_delta'],
             rec['median_cosine'], rec['min_cosine']))
    return rec


print('1-2. does padding to 200 change the pooled embedding of a 150 bp tile? (n=%d)' % N)
a = compare(emb_native, emb_pad_masked, 'native150 vs padded200-MASKED-pool')
b = compare(emb_native, emb_pad_plain, 'native150 vs padded200-plain-pool')

# 3. non-degeneracy
var = emb_native.var(dim=0)
out['3_nondegeneracy'] = {
    'n_dims': int(emb_native.shape[1]),
    'dims_with_var_below_1e-8': int((var < 1e-8).sum().item()),
    'var_min': round(var.min().item(), 8), 'var_median': round(var.median().item(), 6),
    'var_max': round(var.max().item(), 6),
    'mean_pairwise_cosine_between_tiles': round(
        torch.nn.functional.cosine_similarity(
            emb_native[:128].unsqueeze(1), emb_native[:128].unsqueeze(0), dim=2
        ).fill_diagonal_(float('nan')).nanmean().item(), 6),
}
print('\n3. non-degeneracy of the 128-d masked mean-pooled embedding')
for k, v in out['3_nondegeneracy'].items():
    print('   %-38s %s' % (k, v))

print('\n--- verdict ---')
PAD_TOL = 1e-3
if a['max_abs_delta'] < PAD_TOL:
    print('PADDING IS INERT under masked pooling (max|delta| %.2e < %.0e).' % (a['max_abs_delta'], PAD_TOL))
    print('-> native 150 and padded-200-masked are interchangeable; prefer NATIVE 150 (simpler,')
    print('   injects no mask tokens). The padded cell becomes a formality, not a design risk.')
else:
    print('PADDING CONTAMINATES the embedding (max|delta| %.2e >= %.0e) EVEN WITH MASKED POOLING,'
          % (a['max_abs_delta'], PAD_TOL))
    print('-> attention lets the 50 masked columns influence the real positions. Run the arm at')
    print('   NATIVE 150 as primary; the "padded to 200" in the arm name must be reported as a')
    print('   declared-but-superseded geometry, with this measurement as the reason.')
if out['3_nondegeneracy']['dims_with_var_below_1e-8'] > 0:
    print('WARNING: %d of 128 dims are constant across tiles -- effective width is lower than 128.'
          % out['3_nondegeneracy']['dims_with_var_below_1e-8'])

p = pathlib.Path(os.environ.get('PA_OUT', 'wiki/promoteratlas_padding_probe_2026-11-01.json'))
p.write_text(json.dumps(out, indent=1), encoding='utf-8')
print('\nwrote', p)
