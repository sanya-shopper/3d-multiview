#ifndef MV_TSDF_H
#define MV_TSDF_H

/* Paper: doc/multiview.tex, section "The TSDF, in detail" (fusion
 * interface, update equation, voxel-major fast path, budget example).
 * OWNERSHIP (parallel build): this header, src/tsdf.c, and
 * tests/test_tsdffast.c belong to the TSDF-completion work item ONLY. */

#include "mv/cam.h"

/* Truncated signed distance function over a fixed axis-aligned volume.
 * Fusion is a pure, commutative fold over enriched samples
 * (point, origin, weight): each sample updates the voxels in the
 * truncation band along its ray by the weighted running average of the
 * signed along-ray offset (positive in front of the surface). The
 * sample fold is band-only: free-space carving beyond +tau is applied
 * by the depth-map fast path below, not by mv_tsdf_fuse -- band-only is
 * sufficient for static-scene reconstruction. Mesh extraction is the
 * table-free marching-tetrahedra variant of marching cubes (six
 * tetrahedra per cube around the main diagonal, linear interpolation on
 * edges, triangles oriented outward).
 * Note: the band walk steps at 0.7 voxel, so one ray may update a voxel
 * twice with correlated values -- weights are mildly optimistic; treat W
 * as relative confidence, not an exact observation count.
 *
 * Setup & run: mv_tsdf_init sizes the volume (origin corner, extents,
 * voxel, tau ~ 3 sigma_Z), mv_tsdf_fuse folds in samples in any order,
 * mv_tsdf_write_ply extracts the mesh. Exercise end to end with
 *     make demo_tsdf && ./demo_tsdf
 * (fuses depth-law-noisy views of a known sphere at two frame counts,
 * scores surface RMS and the sqrt(n) scaling; writes out_tsdf.ply --
 * open in MeshLab/CloudCompare/Blender). Unit coverage: make check
 * (sign, unobserved, and noiseless-mesh assertions). */

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

/* ---- voxel-major fast path (stationary camera) --------------------- */

/* Because the rig never moves, each voxel's projection into a camera --
 * the depth-map pixel it reads and its distance along the viewing ray
 * -- is constant.  mv_tsdf_table_build computes it once; each frame is
 * then one depth lookup + compare + multiply-add per voxel, including
 * the free-space carving (clamped +tau update) that the band-only
 * sample fold does not apply.  Semantics per voxel and frame:
 *     eta = Z*scale - lv        (along-ray signed offset, + in front)
 *     eta < -tau                -> untouched (occluded, no information)
 *     eta > +tau                -> carved with +tau
 *     else                      -> band, true signed offset
 * folded by the same running average as mv_tsdf_fuse.  The along-ray
 * metric matches the sample fold exactly on-axis and everywhere on the
 * ray; it exceeds the Euclidean distance to a slanted surface by the
 * incidence secant, same as the fold. */

typedef struct {
    long nvox;    /* nx*ny*nz of the volume the table was built for */
    int w, h;     /* depth-map size the pixel indices assume */
    int *pix;     /* per-voxel pixel index v*w + u; -1 = no projection */
    float *lv;    /* camera-center-to-voxel distance along the ray */
    float *scale; /* along-ray length per unit depth (>= 1) */
} mv_tsdf_table;

/* Build the projection table of volume t into camera cam for w x h
 * depth maps (applies cam's distortion; nearest-pixel rounding).
 * Voxels behind the camera or off the image get pix = -1. */
int mv_tsdf_table_build(const mv_tsdf *t, const mv_camera *cam,
                        int w, int h, mv_tsdf_table *tab);
void mv_tsdf_table_free(mv_tsdf_table *tab);

/* Fuse one w x h row-major depth map (metres; depth <= 0, NaN or inf =
 * no measurement at that pixel).  Per-sample weight is 1/sigma_Z^2 with
 * the stereo depth law sigma_Z = sigma0 * Z^2; pass sigma0 <= 0 for
 * unit weights.  Returns MV_ERR on volume/table/image size mismatch. */
int mv_tsdf_fuse_depthmap(mv_tsdf *t, const mv_tsdf_table *tab,
                          const float *depth, int w, int h, double sigma0);

#endif /* MV_TSDF_H */
