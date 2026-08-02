#include <math.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/rectify.h"

int mv_rectify_pair(const mv_camera *c1, const mv_camera *c2,
                    mv_camera *r1, mv_camera *r2,
                    double H1[9], double H2[9])
{
    double C1[3], C2[3], x[3], y[3], z[3], k[3];
    double Rnew[9], Knew[9], KR[9];
    const mv_camera *orig[2];
    mv_camera *rect[2];
    double *H[2];
    int i, v;

    mv_cam_center(C1, c1);
    mv_cam_center(C2, c2);

    /* new x axis: the baseline */
    for (i = 0; i < 3; i++)
        x[i] = C2[i] - C1[i];
    if (mv_normalize(x, 3) != MV_OK)
        return MV_ERR;

    /* new y: perpendicular to old optical axis of camera 1 and to x */
    k[0] = c1->R[6]; k[1] = c1->R[7]; k[2] = c1->R[8];
    mv_cross3(y, k, x);
    if (mv_normalize(y, 3) != MV_OK)
        return MV_ERR;
    mv_cross3(z, x, y);

    for (i = 0; i < 3; i++) {
        Rnew[0 + i] = x[i];
        Rnew[3 + i] = y[i];
        Rnew[6 + i] = z[i];
    }

    /* shared intrinsics: average of the two, zero skew */
    for (i = 0; i < 9; i++)
        Knew[i] = 0.5 * (c1->K[i] + c2->K[i]);
    Knew[1] = 0.0;

    orig[0] = c1; orig[1] = c2;
    rect[0] = r1; rect[1] = r2;
    H[0] = H1;   H[1] = H2;
    for (v = 0; v < 2; v++) {
        double C[3], Rt[9], Kinv[9], tmp[9];
        memcpy(rect[v]->K, Knew, sizeof(Knew));
        memcpy(rect[v]->R, Rnew, sizeof(Rnew));
        memset(rect[v]->k, 0, sizeof(rect[v]->k));
        mv_cam_center(C, orig[v]);
        for (i = 0; i < 3; i++)
            rect[v]->t[i] = -(Rnew[i * 3 + 0] * C[0] + Rnew[i * 3 + 1] * C[1]
                              + Rnew[i * 3 + 2] * C[2]);
        /* H = Knew Rnew R^T K^-1 (rotation about the shared center) */
        mv_mat_transpose(Rt, orig[v]->R, 3, 3);
        if (mv_mat_inv3(Kinv, orig[v]->K) != MV_OK)
            return MV_ERR;
        mv_mat_mul(KR, Knew, Rnew, 3, 3, 3);
        mv_mat_mul(tmp, KR, Rt, 3, 3, 3);
        mv_mat_mul(H[v], tmp, Kinv, 3, 3, 3);
    }
    return MV_OK;
}

double mv_baseline(const mv_camera *c1, const mv_camera *c2)
{
    double C1[3], C2[3], d[3];
    int i;
    mv_cam_center(C1, c1);
    mv_cam_center(C2, c2);
    for (i = 0; i < 3; i++)
        d[i] = C2[i] - C1[i];
    return mv_norm(d, 3);
}

double mv_disp_to_depth(double focal_px, double baseline, double disp)
{
    if (disp <= 0.0)
        return -1.0;
    return focal_px * baseline / disp;
}
