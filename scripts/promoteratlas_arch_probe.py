"""PromoterAtlas base-model FEASIBILITY GATE for the frozen-embeddings ceiling arm.

Run BEFORE freezing any acceptance bar. This project has repeated one lesson four times:
a threshold set before asking what the arithmetic and the substrate permit is a threshold that
has to be amended mid-run. So this answers, by measurement:

  A. Does the architecture I inferred from source load the published weights EXACTLY?
     (strict state_dict load -- no missing, no unexpected keys)
  B. Is the base model really length-agnostic? Latent shape at 86 / 150 / 200 / 300 bp.
  C. Is 150-200 bp IN-DISTRIBUTION? Measured via masked-nucleotide reconstruction accuracy
     vs length -- the model's own pretraining objective.
  D. Is it deterministic? (same input twice, bit-identical)
  E. NON-VACUITY: does real sequence reconstruct better than shuffled? If not, the model is
     doing nothing here and C says nothing.

Nothing is installed. dna_transformer.py imports only torch, so it is loaded as a standalone
module file -- the package __init__ of a downloaded repo never runs.
"""
import os
#: Where the public MIT PromoterAtlas clone lives. Override with PA_REPO.
_PA_DEFAULT = r'C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas'
PA_REPO = os.environ.get('PA_REPO', _PA_DEFAULT)
import importlib.util, json, pathlib, random, sys
import torch

REPO = pathlib.Path(PA_REPO)
MODEL_SRC = REPO / 'src/promoter_atlas/models/dna_transformer.py'
WEIGHTS = REPO / 'trained_weights/base_model/promoteratlas-base.pt'
GENOME = REPO / 'data/genomes/NC_000913.1.gb'
NUC = {'A': 0, 'C': 1, 'G': 2, 'T': 3}          # from src/promoter_atlas/utils/genomics.py
N_POINT_MASKS = int(sys.argv[1]) if len(sys.argv) > 1 else 20   # pretraining used 60 (wm0_pm60)
out = {}

# ---------- load the architecture as a standalone module ----------
spec = importlib.util.spec_from_file_location('pa_dna_transformer', MODEL_SRC)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
model = mod.DNATransformer(dim=128, n_blocks=8, n_heads=8)

# ---------- A. strict load ----------
# The checkpoint was saved on Linux and embeds a pathlib.PosixPath, which (a) is not on
# weights_only's allowlist and (b) cannot instantiate on Windows at all. Aliasing PosixPath does
# not help, because safe_globals keys on the NAME "pathlib.PosixPath", not the resolved object.
#
# So weights_only=False is required here -- and it is justified by MEASUREMENT, not by trusting
# the source. pickle_audit.py disassembled the pickle without executing it and found exactly four
# referenced globals: collections.OrderedDict, pathlib.PosixPath, torch.FloatStorage,
# torch._utils._rebuild_tensor_v2. No callable that could execute arbitrary logic is named.
# Aliasing PosixPath -> WindowsPath is still needed so the (unused) metadata path can instantiate.
pathlib.PosixPath = pathlib.WindowsPath
ck = torch.load(WEIGHTS, map_location='cpu', weights_only=False)
sd = ck['model_state_dict']
res = model.load_state_dict(sd, strict=True)      # raises if my inference is wrong
model.eval()
n_par = sum(p.numel() for p in model.parameters())
out['A_strict_load'] = {
    'ok': True, 'missing': list(res.missing_keys), 'unexpected': list(res.unexpected_keys),
    'n_params': n_par, 'ckpt_extra_keys': [k for k in ck if k != 'model_state_dict'],
}
print('A. strict state_dict load: OK  (%d params, 0 missing, 0 unexpected)' % n_par)
print('   checkpoint also carries:', out['A_strict_load']['ckpt_extra_keys'])

# ---------- real E. coli K-12 sequence ----------
txt = GENOME.read_text(encoding='utf-8', errors='replace')
seq = []
for line in txt.split('ORIGIN', 1)[-1].splitlines():
    s = ''.join(c for c in line if c.isalpha()).upper()
    if s:
        seq.append(s)
genome = ''.join(seq)
print('\nK-12 sequence parsed from the repo GenBank: %d bp' % len(genome))
assert len(genome) > 500_000, 'genome parse looks wrong'
out['genome_bp'] = len(genome)


def onehot(s):
    t = torch.zeros((4, len(s)))
    for i, c in enumerate(s):
        j = NUC.get(c)
        if j is not None:
            t[j, i] = 1.0
    return t


def windows(L, n, seed, shuffle=False):
    rng = random.Random(seed)
    xs = []
    for _ in range(n):
        p = rng.randrange(0, len(genome) - L)
        s = genome[p:p + L]
        if not set(s) <= set('ACGT'):
            s = ''.join(c if c in 'ACGT' else 'A' for c in s)
        if shuffle:
            ls = list(s); rng.shuffle(ls); s = ''.join(ls)
        xs.append(s)
    return xs


@torch.no_grad()
def recon_accuracy(L, n=24, seed=11, shuffle=False):
    """Mask N_POINT_MASKS positions per sequence (held CONSTANT across lengths so the task is
    comparable), then score argmax at exactly those positions."""
    rng = random.Random(seed + 7919)
    seqs = windows(L, n, seed, shuffle=shuffle)
    X = torch.stack([onehot(s) for s in seqs])
    truth = X.argmax(dim=1)
    pos = [rng.sample(range(L), N_POINT_MASKS) for _ in range(n)]
    Xm = X.clone()
    for b, ps in enumerate(pos):
        for p in ps:
            Xm[b, :, p] = 0.0
    logits, _, latent = model(Xm, return_latent=True)
    correct = tot = 0
    for b, ps in enumerate(pos):
        for p in ps:
            correct += int(logits[b, :, p].argmax().item() == truth[b, p].item())
            tot += 1
    return correct / tot, tuple(latent.shape), tuple(logits.shape)


# ---------- B + C + E ----------
print('\nB/C/E. length-agnosticism, in-distribution length, and the non-vacuity control')
print('   %-7s %-22s %-12s %-12s %s' % ('len', 'latent shape', 'real acc', 'shuffled', 'real-shuffled'))
rows = {}
for L in (86, 150, 200, 300):
    a_real, lat, lg = recon_accuracy(L)
    a_shuf, _, _ = recon_accuracy(L, shuffle=True)
    rows[L] = {'latent_shape': lat, 'logits_shape': lg,
               'recon_acc_real': round(a_real, 4), 'recon_acc_shuffled': round(a_shuf, 4),
               'delta': round(a_real - a_shuf, 4)}
    print('   %-7d %-22s %-12.4f %-12.4f %+.4f' % (L, str(lat), a_real, a_shuf, a_real - a_shuf))
out['BCE_by_length'] = rows

# ---------- D. determinism ----------
X = torch.stack([onehot(s) for s in windows(200, 4, 5)])
with torch.no_grad():
    l1, _, z1 = model(X, return_latent=True)
    l2, _, z2 = model(X, return_latent=True)
dz = (z1 - z2).abs().max().item()
out['D_determinism_max_abs_latent_delta'] = dz
print('\nD. determinism: max |latent delta| over a repeat forward = %.3e  -> %s'
      % (dz, 'DETERMINISTIC' if dz == 0.0 else 'NOT bit-identical'))

# ---------- verdicts, stated as conditions not conclusions ----------
v = {}
v['architecture_loads_exactly'] = not (res.missing_keys or res.unexpected_keys)
v['length_agnostic'] = all(rows[L]['latent_shape'] == (24, 128, L) for L in rows)
v['model_is_doing_something'] = all(rows[L]['delta'] > 0.02 for L in rows)
best = max(rows, key=lambda L: rows[L]['recon_acc_real'])
v['best_length_by_reconstruction'] = best
v['acc_at_150_minus_acc_at_86'] = round(rows[150]['recon_acc_real'] - rows[86]['recon_acc_real'], 4)
v['acc_at_200_minus_acc_at_86'] = round(rows[200]['recon_acc_real'] - rows[86]['recon_acc_real'], 4)
out['verdicts'] = v

print('\n--- feasibility verdicts ---')
for k, val in v.items():
    print('   %-34s %s' % (k, val))

p = pathlib.Path(os.environ.get('PA_OUT', 'wiki/promoteratlas_arch_probe_2026-11-01.json'))
p.write_text(json.dumps(out, indent=1), encoding='utf-8')
print('\nwrote', p)
if not v['model_is_doing_something']:
    print('\nREFUSED: real sequence does not reconstruct better than shuffled at every length.')
    print('The reconstruction metric is then uninformative and says nothing about in-distribution')
    print('length. Do NOT freeze a bar on it.')
    sys.exit(3)
