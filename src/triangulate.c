#include <math.h>
#include <stdlib.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/triangulate.h"

int mv_triangulate(double X[3], const mv_camera *const *cams,
                   const double *uv, int nviews)
{
    double *A;
    double x4[4];
    int v, j, ret;

    if (nviews < 2)
        return MV_ERR;
    A = (double *)malloc((size_t)(2 * nviews) * 4 * sizeof(double));
    if (!A)
        return MV_ERR;

    for (v = 0; v < nviews; v++) {
        double P[12];
        double u = uv[2 * v], w = uv[2 * v + 1];
        double *r0 = A + (2 * v) * 4;
        double *r1 = A + (2 * v + 1) * 4;
        mv_cam_projection(P, cams[v]);
        for (j = 0; j < 4; j++) {
            r0[j] = u * P[8 + j] - P[0 + j];
            r1[j] = w * P[8 + j] - P[4 + j];
        }
    }
    ret = mv_nullvec(x4, A, 2 * nviews, 4);
    free(A);
    if (ret != MV_OK || !(fabs(x4[3]) >= 1e-14))
        return MV_ERR; /* the >= form also rejects a NaN homogeneous w */
    X[0] = x4[0] / x4[3];
    X[1] = x4[1] / x4[3];
    X[2] = x4[2] / x4[3];
    if (!(isfinite(X[0]) && isfinite(X[1]) && isfinite(X[2])))
        return MV_ERR;
    return MV_OK;
}

double mv_reproj_rms(const mv_camera *const *cams, const double *uv,
                     int nviews, const double X[3])
{
    double se = 0.0;
    int v, valid = 0;
    for (v = 0; v < nviews; v++) {
        double p[2], du, dv;
        if (mv_cam_project(p, cams[v], X) != MV_OK)
            continue;
        du = p[0] - uv[2 * v];
        dv = p[1] - uv[2 * v + 1];
        se += du * du + dv * dv;
        valid++;
    }
    if (!valid)
        return HUGE_VAL;
    return sqrt(se / valid);
}
