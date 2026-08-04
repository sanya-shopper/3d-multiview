#ifndef MV_TSDF_H
#define MV_TSDF_H

/* Paper: doc/multiview.tex, section "The TSDF, in detail" (fusion
 * interface, update equation, budget example). */

/* Truncated signed distance function over a fixed axis-aligned volume.
 * Fusion is a pure, commutative fold over enriched samples
 * (point, origin, weight): each sample updates the voxels in the
 * truncation band along its ray by the weighted running average of the
 * signed along-ray offset (positive in front of the surface). Band-only
 * in v1: free-space carving beyond +tau is not yet applied, which is
 * sufficient for static-scene reconstruction. Mesh extraction is the
 * table-free marching-tetrahedra variant of marching cubes (six
 * tetrahedra per cube around the main diagonal, linear interpolation on
 * edges, triangles oriented outward). */

typedef struct {
    int nx, ny, nz;
    double x0, y0, z0;  /* corner of voxel (0,0,0) */
    double voxel;       /* edge length, metres */
    double tau;         /* truncation distance, metres */
    float *d;           /* nx*ny*nz signed distances */
    float *w;           /* accumulated weights, 0 = unobserved */
} mv_tsdf;

int mv_tsdf_init(mv_tsdf *t, double x0, double y0, double z0,
                 double sx, double sy, double sz, double voxel, double tau);
void mv_tsdf_free(mv_tsdf *t);

/* Fuse one enriched sample: surface point p observed from origin o with
 * weight wgt (use 1/sigma_Z^2, or 1.0 for unweighted). */
void mv_tsdf_fuse(mv_tsdf *t, const double p[3], const double o[3],
                  double wgt);

/* Trilinearly interpolated signed distance at q; HUGE_VAL if any of the
 * 8 surrounding voxels is unobserved or q is outside the volume. */
double mv_tsdf_query(const mv_tsdf *t, const double q[3]);

/* Extract the zero level set as a triangle soup: 9 doubles per triangle
 * (three xyz vertices), malloc'd into *tris (caller frees). Only cubes
 * with all 8 corners observed contribute. Returns MV_OK. */
int mv_tsdf_mesh(const mv_tsdf *t, double **tris, int *ntri);

/* Extract and write an ASCII PLY mesh (vertices + faces). */
int mv_tsdf_write_ply(const mv_tsdf *t, const char *path);

#endif /* MV_TSDF_H */
