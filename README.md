# multiview

A portable, dependency-free **C99** library for multiview-geometry scene
analysis and 3D model extraction from the video of **two stationary
cameras** observing an overlapping volume — optionally augmented by a
collection of still photos, and later by range-finder data.

Calibration is part of the system: the rig calibrates itself (Zhang's
method) from views of a **printed letter-page checkerboard** of known square
size, or from a **computer display** showing synchronized Gray-code
patterns (dense correspondences with no corner detector).

Full documentation — capabilities, the underlying theory as a self-study
overview with hand-computable worked examples, numerical notes, measured
results, and bibliography — is in [`doc/multiview.tex`](doc/multiview.tex)
(build with `make doc`, requires pdflatex + bibtex). A running log of
development sessions is in [`CONVERSATION_LOG.md`](CONVERSATION_LOG.md).

## Layout

- `include/mv/` — public headers (`mv/mv.h` is the umbrella); one header per module
- `src/` — implementations: linear algebra + Jacobi SVD, camera model,
  planar targets + Gray-code display patterns, Zhang calibration,
  the spec-v1 calibration screen (M-array pattern generator), a synthetic
  plane renderer, the blind pattern reader,
  epipolar geometry, DLT triangulation, rectification, dense stereo,
  PGM image I/O, PLY point clouds
- `demo/synthetic.c` — deterministic two-camera reconstruction experiment
- `demo/calibrate.c` — deterministic calibration experiment; also renders
  the printable target (`target_letter.pgm`)
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
