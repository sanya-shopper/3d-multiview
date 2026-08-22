# multiview

A portable, dependency-free **C99** library for multiview-geometry scene
analysis and 3D model extraction from the video of **fixed (stationary)
cameras** — two in the reference rig — observing an overlapping volume;
optionally augmented by a collection of still photos, later by
range-finder data, and extensible to tracked mobile cameras that borrow
their poses from the fixed backbone.

**Interactive companions** (GitHub Pages):
[two cameras watching a molecule](https://sanya-shopper.github.io/3d-multiview/web/index.html)
(explore a solved rig) and
[earn the rig](https://sanya-shopper.github.io/3d-multiview/web/rig.html)
(solve one yourself) — see
[the landing page](https://sanya-shopper.github.io/3d-multiview/) for both,
plus the built document.

Calibration is part of the system: the rig calibrates itself (Zhang's
method) from views of a **printed letter-page checkerboard** of known square
size, or from a **computer display** showing synchronized Gray-code
patterns (dense correspondences with no corner detector).

Full documentation — capabilities, the underlying theory as a self-study
overview with hand-computable worked examples, numerical notes, measured
results, and bibliography — is in [`doc/multiview.tex`](doc/multiview.tex)
(build with `make doc`, requires pdflatex + bibtex). A running log of
development sessions is in [`CONVERSATION_LOG.md`](CONVERSATION_LOG.md).

Code and paper cross-reference each other: every public header opens
with a `Paper:` comment naming the section it implements, and the
paper's final appendix ("Where to read: the code–paper map") maps each
section to its module and its validating test or experiment.

## Layout

- `include/mv/` — public headers (`mv/mv.h` is the umbrella); one header per module
- `src/` — implementations: linear algebra + Jacobi SVD, camera model,
  planar targets + Gray-code display patterns, Zhang calibration,
  the spec-v1 calibration screen (M-array pattern generator), a synthetic
  plane renderer, the blind pattern reader, TSDF fusion with mesh extraction,
  epipolar geometry, DLT triangulation, rectification, dense stereo,
  PGM image I/O, PLY point clouds
- `demo/synthetic.c` — deterministic two-camera reconstruction experiment
- `demo/calibrate.c` — deterministic calibration experiment; also renders
  the printable target (`target_letter.pgm`)
- further deterministic experiments: tracking (`track_robot`, `track_insects`,
  `track_people`), photometric analysis (`lightlog`), error diagnostics
  (`diagnose`), the calibration-pattern pipeline (`patternsim`), TSDF
  fusion (`tsdfsim`), and the dense room model with change detection
  (`roomsim`) — `make demo` runs them all
- `tests/` — 17 C suites (`test_mv.c` is the core exact-property suite
  on noiseless data), two headless JS suites for the web apps, and
  Python consistency checks (build targets, bibliography, web stamp)
- `doc/` — self-contained LaTeX documentation with local bibliography

## Build & verify

```sh
make          # builds libmv.a, the demo, and the tests (strict C99, libm only)
make check    # runs all test suites and consistency checks (bib check needs TeX)
make demo     # runs all the deterministic experiments (writes out_cloud.ply etc.)
make doc      # builds doc/multiview.pdf
```

Everything is deterministic: identical inputs produce bit-identical outputs,
so all reported numbers are reproducible.
