#include <math.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/cam.h"

void mv_cam_set_K(mv_camera *c, double fx, double fy, double cx, double cy)
{
    memset(c->K, 0, sizeof(c->K));
    c->K[0] = fx;
    c->K[2] = cx;
    c->K[4] = fy;
    c->K[5] = cy;
    c->K[8] = 1.0;
}

void mv_cam_set_identity_pose(mv_camera *c)
{
    memset(c->R, 0, sizeof(c->R));
    c->R[0] = c->R[4] = c->R[8] = 1.0;
    c->t[0] = c->t[1] = c->t[2] = 0.0;
}

void mv_cam_set_pose_yaw(mv_camera *c, double angle, const double C[3])
{
    double ca = cos(angle), sa = sin(angle);
    int i;
    c->R[0] = ca;  c->R[1] = 0.0; c->R[2] = -sa;
    c->R[3] = 0.0; c->R[4] = 1.0; c->R[5] = 0.0;
    c->R[6] = sa;  c->R[7] = 0.0; c->R[8] = ca;
    /* t = -R C */
    for (i = 0; i < 3; i++)
        c->t[i] = -(c->R[i * 3 + 0] * C[0] + c->R[i * 3 + 1] * C[1]
                    + c->R[i * 3 + 2] * C[2]);
}

void mv_cam_projection(double P[12], const mv_camera *c)
{
    double KR[9], Kt[3];
    int i;
    mv_mat_mul(KR, c->K, c->R, 3, 3, 3);
    mv_mat_mul(Kt, c->K, c->t, 3, 3, 1);
    for (i = 0; i < 3; i++) {
        P[i * 4 + 0] = KR[i * 3 + 0];
        P[i * 4 + 1] = KR[i * 3 + 1];
        P[i * 4 + 2] = KR[i * 3 + 2];
        P[i * 4 + 3] = Kt[i];
    }
}

void mv_cam_center(double C[3], const mv_camera *c)
{
    /* C = -R^T t */
    int i;
    for (i = 0; i < 3; i++)
        C[i] = -(c->R[0 * 3 + i] * c->t[0] + c->R[1 * 3 + i] * c->t[1]
                 + c->R[2 * 3 + i] * c->t[2]);
}

static void distort(double *xd, double *yd, const double k[5],
                    double x, double y)
{
    double r2 = x * x + y * y;
    double radial = 1.0 + r2 * (k[0] + r2 * (k[1] + r2 * k[4]));
    *xd = x * radial + 2.0 * k[2] * x * y + k[3] * (r2 + 2.0 * x * x);
    *yd = y * radial + k[2] * (r2 + 2.0 * y * y) + 2.0 * k[3] * x * y;
}

int mv_cam_project(double uv[2], const mv_camera *c, const double X[3])
{
    double Xc[3], xn, yn, xd, yd;
    int i;
    for (i = 0; i < 3; i++)
        Xc[i] = c->R[i * 3 + 0] * X[0] + c->R[i * 3 + 1] * X[1]
              + c->R[i * 3 + 2] * X[2] + c->t[i];
    /* reject non-finite depth explicitly: NaN <= 1e-12 is false, so a
     * NaN pose/point would otherwise pass and emit NaN pixels */
    if (!(Xc[2] > 1e-12))
        return MV_ERR;
    xn = Xc[0] / Xc[2];
    yn = Xc[1] / Xc[2];
    distort(&xd, &yd, c->k, xn, yn);
    uv[0] = c->K[0] * xd + c->K[1] * yd + c->K[2];
    uv[1] = c->K[4] * yd + c->K[5];
    return MV_OK;
}

int mv_cam_ray(double orig[3], double dir[3], const mv_camera *c,
               const double uv[2])
{
    double Kinv[9], xn, yn, x, y, d[3];
    int i, it;
    if (mv_mat_inv3(Kinv, c->K) != MV_OK)
        return MV_ERR;
    xn = Kinv[0] * uv[0] + Kinv[1] * uv[1] + Kinv[2];
    yn = Kinv[3] * uv[0] + Kinv[4] * uv[1] + Kinv[5];
    /* fixed-point undistortion: find (x,y) with distort(x,y) = (xn,yn) */
    x = xn;
    y = yn;
    for (it = 0; it < 20; it++) {
        double r2 = x * x + y * y;
        double radial = 1.0 + r2 * (c->k[0] + r2 * (c->k[1] + r2 * c->k[4]));
        double dx = 2.0 * c->k[2] * x * y + c->k[3] * (r2 + 2.0 * x * x);
        double dy = c->k[2] * (r2 + 2.0 * y * y) + 2.0 * c->k[3] * x * y;
        x = (xn - dx) / radial;
        y = (yn - dy) / radial;
    }
    d[0] = x;
    d[1] = y;
    d[2] = 1.0;
    /* world direction = R^T d */
    for (i = 0; i < 3; i++)
        dir[i] = c->R[0 * 3 + i] * d[0] + c->R[1 * 3 + i] * d[1]
               + c->R[2 * 3 + i] * d[2];
    if (mv_normalize(dir, 3) != MV_OK)
        return MV_ERR;
    mv_cam_center(orig, c);
    return MV_OK;
}
