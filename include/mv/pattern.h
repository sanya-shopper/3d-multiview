#ifndef MV_PATTERN_H
#define MV_PATTERN_H

/* Paper: doc/multiview.tex, section "The calibration subsystem" (pattern spec v1, M-array
 * construction). */

/* Pattern spec v1 (doc section "The calibration subsystem"): the composite
 * calibration screen for a 1920x1080 reference display. The M-array is a
 * fixed constant (generated once by seeded search, verified: every 4x4
 * window is unique under all four rotations with pairwise Hamming
 * distance >= 2), so any decoded window fixes position AND orientation.
 *
 * Setup & run: mv_pattern_render emits full frames for a given counter
 * value (the display app will blit these); mv_pattern_selftest verifies
 * the shipped M-array exhaustively. Exercised end to end by
 *     make demo_patternsim && ./demo_patternsim
 * and covered by make check. */

#define MV_PAT_SPEC_VERSION 1
#define MV_PAT_W 1920
#define MV_PAT_H 1080

#define MV_PAT_GRID_COLS 19   /* checkerboard squares */
#define MV_PAT_GRID_ROWS 10
#define MV_PAT_CELL 80        /* px per square */
#define MV_PAT_GRID_X0 192
#define MV_PAT_GRID_Y0 160
#define MV_PAT_CORNER_COLS 18 /* inner corner lattice */
#define MV_PAT_CORNER_ROWS 9
#define MV_PAT_DOT 24         /* M-array mark size */
#define MV_PAT_WIN 4          /* M-array window */

#define MV_PAT_CTR_X0 192
#define MV_PAT_CTR_Y0 984
#define MV_PAT_CTR_CELL 64
#define MV_PAT_CTR_BITS 20    /* Gray-coded counter bits */
#define MV_PAT_CTR_CELLS 24   /* [1,0] marker + 20 bits + parity + ~parity */

/* M-array bit of grid square (col 0..18, row 0..9). */
int mv_marray_bit(int col, int row);

/* Look up a 4x4 window code (row-major, MSB = top-left bit). On success
 * returns MV_OK and sets the window's pattern position (col,row of its
 * top-left square) and rot = number of clockwise 90-degree rotations
 * that must be APPLIED TO THE OBSERVED window to reproduce the pattern
 * window (the reader applies the same rot to lattice coordinates).
 * MV_ERR if absent. */
int mv_marray_lookup(unsigned code, int *col, int *row, int *rot);

/* Render the full spec-v1 frame for display-frame counter value.
 * img must hold MV_PAT_W * MV_PAT_H bytes; values are 0 or 255 only. */
void mv_pattern_render(unsigned char *img, unsigned counter);

/* Display-pixel coordinates of inner corner (i 0..17, j 0..8). */
void mv_pattern_corner_px(int i, int j, double xy[2]);

/* Display-pixel center of counter cell c (0..23). */
void mv_pattern_ctr_cell_px(int c, double xy[2]);

/* Verify the shipped M-array properties exhaustively; MV_OK if valid. */
int mv_pattern_selftest(void);

/* ---- Coarse tier (pattern spec v2) --------------------------------
 *
 * Room-range layout: the fine tier above needs ~17 px per 80 px cell
 * on the sensor (usable to ~1.1 m for a phone, ~0.6 m for a 720p
 * webcam -- doc eq. "patrange"); the coarse tier scales every feature
 * up ~3.4x for the same display, extending reach by the same factor.
 * 6x3 checkerboard of 270 px cells (10 inner corners), an asymmetric
 * L of three 90 px marks in the top-left cells for orientation (the
 * mark set maps to a disjoint set under every non-identity rotation),
 * and a 12-cell counter strip along the bottom: [1,0] sync + 8 Gray
 * bits + parity + ~parity, wrapping every 256 refreshes (~4.3 s at
 * 60 Hz; consumers disambiguate wraps by capture-time proximity or
 * rig-rigidity consistency). */

#define MV_PAT2_GRID_COLS 6
#define MV_PAT2_GRID_ROWS 3
#define MV_PAT2_CELL 270
#define MV_PAT2_GRID_X0 150
#define MV_PAT2_GRID_Y0 90
#define MV_PAT2_CORNER_COLS 5 /* inner corner lattice */
#define MV_PAT2_CORNER_ROWS 2
#define MV_PAT2_MARK 90       /* orientation mark size */
#define MV_PAT2_CTR_X0 150
#define MV_PAT2_CTR_Y0 1005  /* strip centered 1.5 lattice steps below
                              * the bottom corner row: the gutters
                              * between adjacent squares are true
                              * saddles, so they must sit where grid
                              * growth (integer steps) cannot reach */
#define MV_PAT2_CTR_CELL 135  /* counter cell width */
#define MV_PAT2_CTR_H 60      /* square height (width MV_PAT2_MARK) */
#define MV_PAT2_CTR_BITS 8
#define MV_PAT2_CTR_CELLS 12  /* [1,0] + 8 bits + parity + ~parity */

/* 1 if grid cell (col 0..5, row 0..2) carries an orientation mark. */
int mv_pattern2_mark(int col, int row);

/* Render the full coarse frame for a counter value (low 8 bits used). */
void mv_pattern2_render(unsigned char *img, unsigned counter);

/* Display-pixel coordinates of coarse inner corner (i 0..4, j 0..1). */
void mv_pattern2_corner_px(int i, int j, double xy[2]);

#endif /* MV_PATTERN_H */
