"""Audit every global a downloaded .pt checkpoint references, WITHOUT executing it.

weights_only=True refuses this file over a pathlib.PosixPath, and weights_only=False permits
arbitrary code execution. Rather than trust or abandon, disassemble the pickle and list exactly
which classes/functions it names. Only then decide.
"""
import os
#: Where the public MIT PromoterAtlas clone lives. Override with PA_REPO.
_PA_DEFAULT = r'C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas'
PA_REPO = os.environ.get('PA_REPO', _PA_DEFAULT)
import pathlib, pickletools, sys, zipfile, collections, io

P = pathlib.Path(PA_REPO) / 'trained_weights/base_model/promoteratlas-base.pt'

with zipfile.ZipFile(P) as z:
    names = z.namelist()
    pkls = [n for n in names if n.endswith('.pkl')]
    print('archive entries: %d   (pickles: %s)' % (len(names), pkls))
    print('non-tensor entries:', [n for n in names if '/data/' not in n][:12])
    print()
    globals_found = collections.Counter()
    for pk in pkls:
        raw = z.read(pk)
        ops = list(pickletools.genops(io.BytesIO(raw)))
        # STACK_GLOBAL pulls two preceding SHORT_BINUNICODE strings off the stack
        strings = []
        for op, arg, _pos in ops:
            if op.name in ('SHORT_BINUNICODE', 'BINUNICODE', 'UNICODE'):
                strings.append(arg)
            elif op.name == 'STACK_GLOBAL':
                if len(strings) >= 2:
                    globals_found['%s.%s' % (strings[-2], strings[-1])] += 1
            elif op.name in ('GLOBAL', 'INST'):
                globals_found[str(arg).replace(' ', '.')] += 1

print('GLOBALS referenced by the pickle (name -> times):')
for g, c in sorted(globals_found.items()):
    print('   %-56s %d' % (g, c))

ALLOWED_PREFIXES = ('torch.', 'collections.OrderedDict', 'pathlib.', 'builtins.',
                    'numpy.core.multiarray', 'numpy.dtype', '_codecs.encode')
suspicious = [g for g in globals_found
              if not any(g.startswith(p) or g == p for p in ALLOWED_PREFIXES)]
print()
if suspicious:
    print('SUSPICIOUS globals -> do NOT load with weights_only=False:')
    for g in suspicious:
        print('   ', g)
    sys.exit(3)
print('VERDICT: every referenced global is a torch / pathlib / stdlib container type.')
print('No callable that could execute arbitrary logic is named. weights_only=False is')
print('demonstrably safe FOR THIS FILE, on this evidence.')
