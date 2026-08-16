# refs retrieval and screening record

Relocated 2026-08-12: the PDFs now live OUTSIDE the repository, in the
sibling folder `../_refs/multiview/` (one shared refs tree across
projects, one subdirectory per repo). This file stays in-repo as the
screening record; bibliography local-copy notes and tests/check_bib.py
point at the new location. The git history was rewritten on the same
day (git filter-repo, bibsrc/*.pdf purged from every commit), so the
repository no longer carries the PDFs at any point in its history; a
pre-rewrite backup bundle is kept locally outside the repo.

Retrieved 2026-08-04 from publisher/author/institutional open-access
sources (arXiv, NIST, university pages) via curl. Exact per-file URLs
were not recorded at retrieval time (an earlier note here claimed the
git history had them; it did not). On 2026-08-12 canonical sources
were re-located and verified BYTE-IDENTICAL by SHA-256 for 26 of the
27 files (the exception, Maze2018, exists publicly only as a same-work
preprint); the verified URLs are recorded in the bibliography notes. VirusTotal was not available in this environment; screening
was heuristic: PDF magic verified, and object-level scan for active
content (JavaScript/Launch actions, embedded files), including inside
decompressed streams. /OpenAction and /AA tokens alone are ordinary
viewer navigation actions and were not treated as findings.

| file | size | active-content findings | embedded files |
|---|---|---|---|
| CurlessLevoy1996.pdf | 677 K | none | no |
| Deng2019.pdf | 11186 K | none | no |
| Deng2020.pdf | 5910 K | none | no |
| FischlerBolles1981.pdf | 1204 K | none | no |
| Ge2021.pdf | 851 K | none | no |
| Grother2019.pdf | 28858 K | none | no |
| Guo2021.pdf | 4955 K | none | no |
| Hartley1997.pdf | 4408 K | none | no |
| HartleySturm1997.pdf | 431 K | none | no |
| Hoppe1992.pdf | 562 K | none | yes |
| Horn1987.pdf | 1436 K | none | no |
| Hornung2013.pdf | 3592 K | none | no |
| Huang2007.pdf | 218 K | none | no |
| Kazhdan2006.pdf | 4463 K | none | yes |
| Keller2013.pdf | 2390 K | none | no |
| KutulakosSeitz2000.pdf | 861 K | none | no |
| Lowe2004.pdf | 444 K | none | no |
| Maze2018.pdf | 4057 K | none | no |
| Newcombe2011.pdf | 8099 K | none | no |
| Redmon2018.pdf | 2398 K | none | no |
| ScharsteinSzeliski2002.pdf | 1470 K | none | no |
| Schroff2015.pdf | 4594 K | none | no |
| Seitz2006.pdf | 1786 K | none | no |
| Snavely2006.pdf | 1681 K | none | no |

Files with real active-content findings: 0.

## Additions 2026-08-11

Same screening method as above (PDF magic verified, object-level scan
for JavaScript/Launch/RichMedia/XFA actions and embedded files,
including inside decompressed streams; VirusTotal unavailable).
Amanatides1987 from the author's York University page; Cadena2016 from
arXiv (1606.05830); Elfes1989 via the Internet Archive mirror of a
Brooklyn College course page --- the live copy has a corrupt xref table
(unreadable by poppler) and was replaced by the archived one, which
renders correctly and matches the published paper.

| file | size | active-content findings | embedded files |
|---|---|---|---|
| Amanatides1987.pdf | 17 K | none | no |
| Cadena2016.pdf | 6535 K | none | no |
| Elfes1989.pdf | 1176 K | none | no |

Files with real active-content findings: 0.
