#ifndef MV_CAM_H
#define MV_CAM_H

/* Calibrated pinhole camera with Brown radial-tangential distortion.
 * World-to-camera: Xc = R * X + t; camera center C = -R^T t. */

typedef struct {
    double K[9]; /* intrinsics [fx s cx; 0 fy cy; 0 0 1], row-major */
    double R[9]; /* world-to-camera rotation, row-major */
    double t[3]; /* translation, Xc = R X + t */
    double k[5]; /* distortion k1 k2 p1 p2 k3 (all zero = ideal pinhole) */
} mv_camera;

/* Convenience initializers. */
void mv_cam_set_K(mv_camera *c, double fx, double fy, double cx, double cy);
void mv_cam_set_identity_pose(mv_camera *c);
/* R = rotation of `angle` radians about the world y axis, centered at C. */
void mv_cam_set_pose_yaw(mv_camera *c, double angle, const double C[3]);

/* P (3x4) = K [R | t]. Ignores distortion. */
void mv_cam_projection(double P[12], const mv_camera *c);

void mv_cam_center(double C[3], const mv_camera *c);

/* Project world point X to pixel uv (applies distortion).
 * Returns MV_ERR if X is at or behind the camera plane. */
int mv_cam_project(double uv[2], const mv_camera *c, const double X[3]);

/* Back-project pixel uv to a world-space ray (undistorts iteratively).
 * orig = camera center, dir = unit direction. */
int mv_cam_ray(double orig[3], double dir[3], const mv_camera *c,
               const double uv[2]);

#endif /* MV_CAM_H */
