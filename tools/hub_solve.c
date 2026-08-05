/* Rig pair solver -- v0 is the chordal-mean + per-axis-median code
 * moved verbatim out of livehub.c so the joint-solve work item can
 * replace it (bundle-adjust / Huber-IRLS over all observations)
 * without touching the frozen hub.
 * OWNERSHIP (parallel build): joint-solve work item ONLY. */

#include <math.h>
#include <string.h>

#include "mv/mv.h"

#include "hub_solve.h"

#define SOLVE_MAXOBS 256

int hub_solve_pair(double R[9], double t[3], double *dev_mm,
                   const hub_obs *obs, int n)
{
    double Rsum[9] = { 0 }, tacc[3][SOLVE_MAXOBS];
    double U[9], S[3], V[9], Vt[9];
    double dev = 0.0;
    int i, k;

    if (n < 3)
        return -1;
    if (n > SOLVE_MAXOBS)
        n = SOLVE_MAXOBS;
    for (i = 0; i < n; i++) {
        double Rbt[9], Ri[9], ti[3];
        mv_mat_transpose(Rbt, obs[i].Rb, 3, 3);
        mv_mat_mul(Ri, obs[i].Ra, Rbt, 3, 3, 3);
        mv_mat_mul(ti, Ri, obs[i].tb, 3, 3, 1);
        for (k = 0; k < 3; k++)
            ti[k] = obs[i].ta[k] - ti[k];
        for (k = 0; k < 9; k++)
            Rsum[k] += Ri[k];
        for (k = 0; k < 3; k++)
            tacc[k][i] = ti[k];
    }
    memcpy(U, Rsum, sizeof(Rsum));
    if (mv_svd(U, S, V, 3, 3) != MV_OK)
        return -1;
    mv_mat_transpose(Vt, V, 3, 3);
    mv_mat_mul(R, U, Vt, 3, 3, 3);
    for (k = 0; k < 3; k++) {
        double col[SOLVE_MAXOBS];
        int q, r;
        memcpy(col, tacc[k], (size_t)n * sizeof(double));
        for (q = 1; q < n; q++)
            for (r = q; r > 0 && col[r] < col[r - 1]; r--) {
                double tmp = col[r];
                col[r] = col[r - 1];
                col[r - 1] = tmp;
            }
        t[k] = col[n / 2];
    }
    for (i = 0; i < n; i++) {
        double d = 0.0;
        for (k = 0; k < 3; k++)
            d += (tacc[k][i] - t[k]) * (tacc[k][i] - t[k]);
        dev += sqrt(d);
    }
    if (dev_mm)
        *dev_mm = 1000.0 * dev / n;
    return 0;
}
