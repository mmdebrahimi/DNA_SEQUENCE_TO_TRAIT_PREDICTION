"""Is 'PromoterAtlas benchmark genome' the same thing as 'held out of pretraining'?

Decides one sentence in the arm's acceptance bar: whether a result on E. coli tiles must be
labelled 'sequence possibly seen during pretraining'.

Method: the 13 benchmark GenBank files name a BioSample. Resolve each BioSample to its RefSeq
assembly accession (grounded via NCBI, never from memory), then test membership in the
PromoterAtlas base-model corpus file list (pure string membership, no further API calls).

REFUSES a verdict if too few BioSamples resolve -- a failed lookup must not read as 'absent'.
That exact failure already fired once tonight on the 3,371-accession sweep.
"""
import json, re, sys, time, urllib.request, pathlib

BENCH = pathlib.Path(r'C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas/data/genomes')
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'


def get(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode('utf-8', 'replace')
        except Exception as e:
            if attempt == 2:
                return ''
            time.sleep(2)
    return ''


# 1. benchmark genome -> BioSample, read locally
pairs = []
for f in sorted(BENCH.glob('*.gb')):
    if '_PA_' in f.name:
        continue
    head = f.read_text(encoding='utf-8', errors='replace')[:4000]
    m = re.search(r'BioSample:\s*(SAM\w+)', head)
    org = re.search(r'^SOURCE\s+(.+)$', head, re.M)
    pairs.append((f.stem, m.group(1) if m else None,
                  org.group(1).strip() if org else '?'))

print('benchmark genomes: %d' % len(pairs))

# 2. BioSample -> assembly accession, grounded
resolved = {}
for stem, bs, org in pairs:
    if not bs:
        continue
    xml = get('%s/esearch.fcgi?db=assembly&term=%s[BioSample]&retmode=json' % (EUTILS, bs))
    ids = []
    try:
        ids = json.loads(xml)['esearchresult'].get('idlist', [])
    except Exception:
        pass
    if not ids:
        continue
    summ = get('%s/esummary.fcgi?db=assembly&id=%s&retmode=json' % (EUTILS, ids[0]))
    try:
        d = json.loads(summ)['result'][ids[0]]
        acc = d.get('refseq_category') and d.get('assemblyaccession') or d.get('assemblyaccession')
        syn = (d.get('synonym') or {})
        resolved[stem] = (syn.get('refseq') or acc or '').split('.')[0]
    except Exception:
        pass
    time.sleep(0.4)

rate = len(resolved) / max(1, len([p for p in pairs if p[1]]))
print('resolved %d of %d BioSamples (%.0f%%)' % (len(resolved), len([p for p in pairs if p[1]]), 100 * rate))

if rate < 0.70:
    print('\nREFUSED: resolution rate %.0f%% < 70%%. No held-out claim is made.' % (100 * rate))
    sys.exit(3)

# 3. corpus membership, local string test
corpus = json.load(urllib.request.urlopen(
    'https://huggingface.co/api/datasets/LCoppens/PromoterAtlas-data', timeout=60))
accs = {s['rfilename'].split('/')[-1][:-6].split('.')[0]
        for s in corpus['siblings'] if s['rfilename'].endswith('.fasta')}
print('corpus assemblies: %d' % len(accs))
print()

inside = outside = 0
print('%-16s %-52s %-18s %s' % ('benchmark', 'organism', 'assembly', 'in pretraining corpus?'))
print('-' * 112)
for stem, bs, org in pairs:
    acc = resolved.get(stem)
    if not acc:
        print('%-16s %-52s %-18s %s' % (stem[:16], org[:52], '(unresolved)', 'CANNOT SAY'))
        continue
    hit = acc in accs
    inside += hit
    outside += not hit
    print('%-16s %-52s %-18s %s' % (stem[:16], org[:52], acc, 'PRESENT' if hit else 'absent'))

print()
print('benchmark genomes IN the pretraining corpus:     %d' % inside)
print('benchmark genomes ABSENT from the corpus:        %d' % outside)
print()
if inside == 0:
    print('VERDICT: consistent with a HELD-OUT benchmark design -- no benchmark assembly appears in')
    print('the pretraining corpus. Strengthens the separate measurement that K-12 GCF_000005845 is')
    print('absent. Still NOT proof the corpus holds no E. coli at all (different assembly of the')
    print('same species would not be caught by this test).')
elif inside == len(resolved):
    print('VERDICT: benchmark does NOT mean held out -- every resolved benchmark assembly is in the')
    print('corpus. The K-12 absence then carries much less weight and the arm MUST carry a')
    print('sequence-exposure caveat.')
else:
    print('VERDICT: MIXED (%d in / %d out) -- benchmark membership does not cleanly imply held-out.' % (inside, outside))
