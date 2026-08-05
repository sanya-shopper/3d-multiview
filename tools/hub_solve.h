#ifndef HUB_SOLVE_H
#define HUB_SOLVE_H

/* OWNERSHIP (parallel build): this header and tools/hub_solve.c belong
 * to the JOINT-SOLVE work item ONLY.  livehub.c is FROZEN during the
 * parallel phase; it collects gated anchor-pose pairs and calls
 * hub_solve_pair() -- everything downstream of that call is this
 * module's to improve (chordal mean -> robust joint estimate). */

/* one matched observation: the display-plane pose seen by camera a and
 * by camera b at (nearly) the same display counter */
typedef struct {
    double Ra[9], ta[3]; /* display -> camera a */
    double Rb[9], tb[3]; /* display -> camera b */
} hub_obs;

/* Estimate the rigid transform  x_a = R x_b + t  from n >= 3 matched
 * pose pairs.  Writes the translation scatter (mean |t_i - t|, mm) to
 * *dev_mm as the quality figure.  Returns 0 on success, -1 on failure
 * (degenerate input). */
int hub_solve_pair(double R[9], double t[3], double *dev_mm,
                   const hub_obs *obs, int n);

/* Per-observation diagnostics from the robust solve (additive API;
 * hub_solve_pair keeps its frozen signature and behavior). */
typedef struct {
    double weight;  /* final trim weight: 1 = kept, 0 = rejected     */
    double rot_err; /* geodesic rotation residual vs estimate, rad   */
    double t_err;   /* translation residual |t_i - t|, m             */
} hub_solve_diag;

/* Same estimate as hub_solve_pair; additionally fills diag[i] for the
 * first min(n, 256) observations when diag is non-NULL (the caller
 * provides at least that many entries).  Intended for logging which
 * anchor pairs the robust estimator distrusted. */
int hub_solve_pair_diag(double R[9], double t[3], double *dev_mm,
                        hub_solve_diag *diag, const hub_obs *obs, int n);

#endif /* HUB_SOLVE_H */
