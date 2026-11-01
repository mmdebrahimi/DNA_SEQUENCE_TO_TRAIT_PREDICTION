"""Resolve ALL PromoterAtlas base-corpus assemblies to organism names.

Settles the sequence-exposure residual: is any E. coli / Shigella assembly in the pretraining
corpus? If yes, the frozen encoder saw the tile substrate's SEQUENCES (never its labels) and the
arm must say so. If no, the arm evaluates on sequence absent from pretraining.

Earlier attempts and why they failed, so this is not re-litigated:
  - POST /genome/accession/dataset_report            -> 404, wrong endpoint shape
  - GET  /genome/accession/{150 accs}/dataset_report -> 400, URL too long
  - GET  /genome/accession/{2 accs}/dataset_report   -> works
So: small chunks. ~3,371 accessions / 20 = ~169 requests.

REFUSES a verdict below a 90% resolution rate. A failed lookup printing 'no E. coli found' is a
broken lookup wearing the costume of a finding -- that fired twice already tonight.
"""
import json, sys, time, urllib.request, urllib.error, collections

HF = 'https://huggingface.co/api/datasets/LCoppens/PromoterAtlas-data'
BASE = ('https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/%s'
        '/dataset_report?returned_content=ASSEMBLY_INFO&page_size=100')

d = json.load(urllib.request.urlopen(HF, timeout=90))
accs = sorted({s['rfilename'].split('/')[-1][:-6]
               for s in d['siblings'] if s['rfilename'].endswith('.fasta')})
print('corpus assemblies: %d' % len(accs), flush=True)

names, failed = {}, []
CH = 20
for i in range(0, len(accs), CH):
    chunk = accs[i:i + CH]
    ok = False
    for attempt in range(3):
        try:
            r = json.load(urllib.request.urlopen(BASE % ','.join(chunk), timeout=90))
            for rep in r.get('reports', []):
                names[rep.get('accession')] = (rep.get('organism') or {}).get('organism_name', '?')
            ok = True
            break
        except Exception as e:
            if attempt == 2:
                failed.extend(chunk)
                print('  chunk %d failed: %s %s' % (i, type(e).__name__, getattr(e, 'code', '')),
                      file=sys.stderr, flush=True)
            time.sleep(2)
    if i % 400 == 0:
        print('  ... %d/%d resolved' % (len(names), len(accs)), flush=True)
    time.sleep(0.34)

rate = len(names) / len(accs)
print(flush=True)
print('resolved %d / %d  (%.1f%%)   unresolved: %d' % (len(names), len(accs), 100 * rate,
                                                        len(accs) - len(names)), flush=True)

json.dump(names, open(os.environ.get('PA_OUT', 'wiki/promoteratlas_corpus_names.json'), 'w'), indent=0)
print('saved -> C:/Users/LAPTOP~1/AppData/Local/Temp/pa_corpus_names.json', flush=True)

if rate < 0.90:
    print('\nREFUSED: resolution rate %.1f%% < 90%%. No presence/absence claim is made.' % (100 * rate))
    sys.exit(3)

# --- the question ---
genus = collections.Counter(n.split()[0] for n in names.values() if n)
ec = sorted((a, n) for a, n in names.items()
            if n.lower().startswith(('escherichia', 'shigella')))

print()
print('Escherichia / Shigella assemblies in the PromoterAtlas base corpus: %d' % len(ec))
for a, n in ec:
    print('   %-18s %s' % (a, n))
print()
print('K-12 MG1655 reference GCF_000005845 present:',
      any(a.startswith('GCF_000005845') for a in names))
print()
print('top 12 genera in the corpus:')
for g, c in genus.most_common(12):
    print('   %-26s %d' % (g, c))
print()
print('distinct organism names: %d  (corpus claims ">3,000 species")' % len(set(names.values())))
print()
if ec:
    print('VERDICT: SEQUENCE EXPOSURE IS REAL. An E. coli/Shigella assembly is in the pretraining')
    print('corpus, so the frozen encoder saw genomic sequence from this species -- never the tile')
    print('LABELS. The arm must be reported as evaluating on sequence the encoder has seen.')
else:
    print('VERDICT: NO E. coli/Shigella assembly in the pretraining corpus. The arm evaluates on')
    print('sequence absent from pretraining. Note this is an ASSEMBLY-level test: a near-identical')
    print('relative (e.g. another enteric) could still carry homologous sequence.')
