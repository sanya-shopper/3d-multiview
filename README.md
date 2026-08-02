# multiview

A portable, dependency-free **C99** library for multiview-geometry scene
analysis and 3D model extraction from the video of **two stationary,
calibrated cameras** observing an overlapping volume — optionally augmented
by a collection of still photos, and later by range-finder data.

Full documentation — capabilities, underlying theory, numerical notes,
measured results, and bibliography — is in [`doc/multiview.tex`](doc/multiview.tex)
(build with `make doc`, requires pdflatex + bibtex). A running log of
development sessions is in [`CONVERSATION_LOG.md`](CONVERSATION_LOG.md).

## Layout

- `include/mv/` — public headers (`mv/mv.h` is the umbrella); one header per module
- `src/` — implementations: linear algebra + Jacobi SVD, camera model,
  epipolar geometry, DLT triangulation, rectification, dense stereo,
  PGM image I/O, PLY point clouds
- `demo/synthetic.c` — deterministic synthetic two-camera experiment; its
  output is quoted in the documentation's results section
- `tests/test_mv.c` — exact-property unit tests on noiseless data
- `doc/` — self-contained LaTeX documentation with local bibliography

## Build & verify

```sh
make          # builds libmv.a, the demo, and the tests (strict C99, libm only)
make check    # runs the unit tests
make demo     # runs the synthetic experiment, writes out_cloud.ply
make doc      # builds doc/multiview.pdf
```

Everything is deterministic: identical inputs produce bit-identical outputs,
so all reported numbers are reproducible.
