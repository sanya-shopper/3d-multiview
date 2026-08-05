#include <math.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/rot.h"

void mv_rot_exp(double R[9], const double r[3])
{
    double th = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    double kx, ky, kz, c, s, v;
    if (th < 1e-12) {
        /* first-order: I + [r]_x (exact at r = 0) */
        R[0] = 1.0;   R[1] = -r[2]; R[2] = r[1];
        R[3] = r[2];  R[4] = 1.0;   R[5] = -r[0];
        R[6] = -r[1]; R[7] = r[0];  R[8] = 1.0;
        return;
    }
    kx = r[0] / th; ky = r[1] / th; kz = r[2] / th;
    c = cos(th); s = sin(th); v = 1.0 - c;
    R[0] = c + kx * kx * v;       R[1] = kx * ky * v - kz * s;
    R[2] = kx * kz * v + ky * s;
    R[3] = ky * kx * v + kz * s;  R[4] = c + ky * ky * v;
    R[5] = ky * kz * v - kx * s;
    R[6] = kz * kx * v - ky * s;  R[7] = kz * ky * v + kx * s;
    R[8] = c + kz * kz * v;
}

void mv_rot_log(double r[3], const double R[9])
{
    double c = (R[0] + R[4] + R[8] - 1.0) / 2.0, th, s;
    if (c > 1.0)
        c = 1.0;
    if (c < -1.0)
        c = -1.0;
    th = acos(c);
    s = sin(th);
    if (fabs(s) >= 1e-9) {
        r[0] = th * (R[7] - R[5]) / (2.0 * s);
        r[1] = th * (R[2] - R[6]) / (2.0 * s);
        r[2] = th * (R[3] - R[1]) / (2.0 * s);
        return;
    }
    if (c > 0.0) { /* theta ~ 0: no rotation */
        r[0] = r[1] = r[2] = 0.0;
        return;
    }
    /* theta ~ pi: R = 2 a a^T - I, so a_i = sqrt((R_ii + 1)/2); recover
     * the largest component robustly and fix the others' signs from the
     * symmetric off-diagonals (R_ij = 2 a_i a_j) */
    {
        double d0 = (R[0] + 1.0) * 0.5, d1 = (R[4] + 1.0) * 0.5,
               d2 = (R[8] + 1.0) * 0.5, ax, ay, az;
        const double PI = 3.14159265358979323846;
        if (d0 < 0.0) d0 = 0.0;
        if (d1 < 0.0) d1 = 0.0;
        if (d2 < 0.0) d2 = 0.0;
        if (d0 >= d1 && d0 >= d2) {
            ax = sqrt(d0);
            ay = (R[1] + R[3]) / (4.0 * ax);
            az = (R[2] + R[6]) / (4.0 * ax);
        } else if (d1 >= d2) {
            ay = sqrt(d1);
            ax = (R[1] + R[3]) / (4.0 * ay);
            az = (R[5] + R[7]) / (4.0 * ay);
        } else {
            az = sqrt(d2);
            ax = (R[2] + R[6]) / (4.0 * az);
            ay = (R[5] + R[7]) / (4.0 * az);
        }
        r[0] = PI * ax;
        r[1] = PI * ay;
        r[2] = PI * az;
    }
}

int mv_rot_project(double R[9], const double M[9])
{
    double U[9], S[3], V[9], Vt[9], det;
    memcpy(U, M, 9 * sizeof(double));
    if (mv_svd(U, S, V, 3, 3) != MV_OK)
        return MV_ERR;
    mv_mat_transpose(Vt, V, 3, 3);
    mv_mat_mul(R, U, Vt, 3, 3, 3);
    det = R[0] * (R[4] * R[8] - R[5] * R[7])
        - R[1] * (R[3] * R[8] - R[5] * R[6])
        + R[2] * (R[3] * R[7] - R[4] * R[6]);
    if (det < 0.0) { /* reflection: flip the last column of U and redo */
        int i;
        for (i = 0; i < 3; i++)
            U[i * 3 + 2] = -U[i * 3 + 2];
        mv_mat_mul(R, U, Vt, 3, 3, 3);
    }
    return MV_OK;
}
