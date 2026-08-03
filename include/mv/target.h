#ifndef MV_TARGET_H
#define MV_TARGET_H

/* Paper: doc/multiview.tex, section "The printed target". */

/* Planar calibration targets.
 *
 * The standard printed target is a checkerboard on a letter/A4 page:
 * 8 x 10 squares of 24 mm (192 x 240 mm printed area), giving 7 x 9 inner
 * corners. Print without scaling and measure one printed square to confirm
 * the physical size (printers rescale); pass the measured size to
 * mv_target_checkerboard. */

#define MV_TARGET_LETTER_COLS   7     /* inner corners across */
#define MV_TARGET_LETTER_ROWS   9     /* inner corners down */
#define MV_TARGET_LETTER_SQUARE 0.024 /* metres, nominal */

/* Fill obj with the 2D plane coordinates (X,Y in metres, Z=0 implicit) of
 * the inner corners, row-major: obj[2*(r*inner_cols+c)] = c*square, ...
 * obj must hold 2*inner_cols*inner_rows doubles. Returns the corner count. */
int mv_target_checkerboard(double *obj, int inner_cols, int inner_rows,
                           double square);

/* Render a printable checkerboard as a PGM image. square_px sets the print
 * resolution: at 10 px/mm (254 dpi) a 24 mm square is square_px = 240.
 * Returns MV_ERR on I/O failure. */
int mv_target_render_pgm(const char *path, int inner_cols, int inner_rows,
                         int square_px, int margin_px);

#endif /* MV_TARGET_H */
