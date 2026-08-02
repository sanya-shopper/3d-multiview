#ifndef MV_CLOUD_H
#define MV_CLOUD_H

/* Growable 3D point cloud with optional per-point color, and an ASCII PLY
 * writer so results open directly in MeshLab / CloudCompare / Blender. */

typedef struct {
    double *xyz;        /* 3*n doubles */
    unsigned char *rgb; /* 3*n bytes, or NULL if uncolored */
    int n;
    int cap;
    int has_rgb;
} mv_cloud;

void mv_cloud_init(mv_cloud *c, int with_rgb);
/* rgb may be NULL when the cloud is uncolored. */
int mv_cloud_push(mv_cloud *c, const double xyz[3], const unsigned char *rgb);
int mv_cloud_write_ply(const char *path, const mv_cloud *c);
void mv_cloud_free(mv_cloud *c);

#endif /* MV_CLOUD_H */
