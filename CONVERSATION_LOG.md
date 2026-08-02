# Conversation Log

Running log of development sessions in this repo, updated as we go.

## 2026-08-02

- Started the project: initialized the git repository in an empty `multiview` directory.
- Created the GitHub repo and pushed the initial commit (README + this log).
- Decision: keep a running conversation log in `CONVERSATION_LOG.md`, updated and committed as work happens.
- Briefly added a product gap-research doc here by mistake (instructions were meant for the `commercial` repo); removed it and moved the work to `~/Claude/Projects/commercial/RESEARCH.md`.
- Project intent defined: a portable C99 multiview-geometry system for scene analysis and 3D model extraction from two stationary, calibrated cameras (video) with overlapping views; later augmented by photo collections and possibly range-finder data. Self-contained, for later learning/self-development in a controlled environment.
- Built the v0.1.0 core: module design (`include/mv/*.h`), implementations (`src/`) — linear algebra with one-sided Jacobi SVD, camera model with Brown distortion, analytic + 8-point fundamental matrix, N-view DLT triangulation, Fusiello rectification, SAD block-matching stereo, PGM and PLY I/O.
- Validation: 19 exact-property unit tests all pass (`make check`); deterministic synthetic experiment (`make demo`, seed 42) — 250 points, 0.5 m baseline, 0.3 px noise → RMS 3D error 2.65 cm (matches the σ_Z = Z²σ_d/(fB) prediction of ~2.3 cm), epipolar residual 0.34 px, 8-point F within 4.5e-3 of analytic, rectified rows exactly aligned.
- Wrote `doc/multiview.tex` — self-contained: capabilities, theory (camera model → epipolar → triangulation → stereo), planned extensions (photo SfM, range-finder fusion), numerical notes, measured results table, TikZ epipolar figure, local `refs.bib` (18 entries), hyperref cross-linking. Compiles clean (`make doc`).
- Decisions: strict C99 + libm only; row-major doubles; determinism as a feature; Jacobi SVD as the single decomposition kernel; calibration treated as offline input; stationary-rig quantities (F, rectifying homographies) computed once and reused per frame.
