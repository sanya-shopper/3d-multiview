#!/usr/bin/env python3
"""Bibliography consistency unit test (parent CLAUDE.md policy).

Checks, exiting nonzero on any failure:
  1. every refs.bib entry is cited somewhere in the body (via the .aux);
  2. every \\cite key in the .tex resolves to a refs.bib entry;
  3. every entry carries a local-copy note, and every claimed local copy
     exists in bibsrc/ and is a real PDF;
  4. the last bibtex run produced zero warnings (via the .blg).
Run after `make doc`; wired into `make check`.
"""
import os
import re
import sys

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(root)
fail = []

bib = open('doc/refs.bib').read()
keys = set(re.findall(r'@\w+\{(\w+),', bib))

aux_cited = set()
AUX = os.path.join(
    os.environ.get('BUILD_TARGET_PREFIX', os.path.expanduser('~/Claude/Projects')),
    '_buildoutput', '3d-multiview', 'doc', 'multiview.aux')
if os.path.exists(AUX):
    for m in re.findall(r'\\citation\{([^}]+)\}', open(AUX).read()):
        aux_cited.update(m.split(','))
else:
    fail.append('doc/multiview.aux missing (run make doc first)')

tex_cited = set()
for m in re.findall(r'\\cite\{([^}]+)\}', open('doc/multiview.tex').read()):
    tex_cited.update(k.strip() for k in m.split(','))

for k in sorted(keys - aux_cited):
    fail.append(f'bib entry never cited: {k}')
for k in sorted(tex_cited - keys):
    fail.append(f'cite without bib entry: {k}')

# per-entry note discipline
for m in re.finditer(r'@\w+\{(\w+),', bib):
    key = m.group(1)
    j = bib.index('{', m.start())
    depth, k2 = 0, j
    while True:
        if bib[k2] == '{':
            depth += 1
        elif bib[k2] == '}':
            depth -= 1
            if depth == 0:
                break
        k2 += 1
    entry = bib[m.start():k2 + 1]
    nm = re.search(r'note = \{([^}]*)\}', entry)
    if not nm:
        fail.append(f'entry missing local-copy note: {key}')
        continue
    note = nm.group(1)
    pm = re.search(r'Local copy: (\.\./\\\\?_refs/3d-multiview/[\w.-]+\.pdf)', note)
    if pm:
        # The note is typeset, so the underscore in _refs is escaped there;
        # the filesystem path is not.
        path = pm.group(1).replace('\\_', '_')
        if not os.path.exists(path):
            fail.append(f'{key}: claimed local copy missing: {path}')
        elif open(path, 'rb').read(4) != b'%PDF':
            fail.append(f'{key}: local copy is not a PDF: {path}')
    elif 'no local copy' not in note.lower():
        fail.append(f'{key}: note neither names a local copy nor '
                    f'declares its absence: {note!r}')

if os.path.exists('doc/multiview.blg'):
    blg = open('doc/multiview.blg').read()
    wm = re.search(r'warning\$? -- (\d+)', blg)
    nwarn = blg.lower().count('warning--')
    if nwarn:
        fail.append(f'bibtex reported {nwarn} warning(s); see doc/multiview.blg')

if fail:
    print('BIBLIOGRAPHY CHECK FAILED:')
    for f in fail:
        print('  -', f)
    sys.exit(1)
print(f'bibliography ok: {len(keys)} entries, all cited, notes consistent, '
      f'{sum(1 for k in keys if os.path.exists(f"../_refs/3d-multiview/{k}.pdf"))} local copies verified')
