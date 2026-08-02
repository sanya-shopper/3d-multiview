#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/cloud.h"

void mv_cloud_init(mv_cloud *c, int with_rgb)
{
    memset(c, 0, sizeof(*c));
    c->has_rgb = with_rgb ? 1 : 0;
}

int mv_cloud_push(mv_cloud *c, const double xyz[3], const unsigned char *rgb)
{
    if (c->n == c->cap) {
        int ncap = c->cap ? c->cap * 2 : 256;
        double *nxyz = (double *)realloc(c->xyz,
                                         (size_t)ncap * 3 * sizeof(double));
        if (!nxyz)
            return MV_ERR;
        c->xyz = nxyz;
        if (c->has_rgb) {
            unsigned char *nrgb = (unsigned char *)realloc(c->rgb,
                                                           (size_t)ncap * 3);
            if (!nrgb)
                return MV_ERR;
            c->rgb = nrgb;
        }
        c->cap = ncap;
    }
    memcpy(c->xyz + 3 * c->n, xyz, 3 * sizeof(double));
    if (c->has_rgb) {
        if (rgb)
            memcpy(c->rgb + 3 * c->n, rgb, 3);
        else
            memset(c->rgb + 3 * c->n, 200, 3);
    }
    c->n++;
    return MV_OK;
}

int mv_cloud_write_ply(const char *path, const mv_cloud *c)
{
    FILE *f = fopen(path, "w");
    int i;
    if (!f)
        return MV_ERR;
    fprintf(f, "ply\nformat ascii 1.0\nelement vertex %d\n", c->n);
    fprintf(f, "property float x\nproperty float y\nproperty float z\n");
    if (c->has_rgb)
        fprintf(f, "property uchar red\nproperty uchar green\n"
                   "property uchar blue\n");
    fprintf(f, "end_header\n");
    for (i = 0; i < c->n; i++) {
        fprintf(f, "%g %g %g", c->xyz[3 * i], c->xyz[3 * i + 1],
                c->xyz[3 * i + 2]);
        if (c->has_rgb)
            fprintf(f, " %d %d %d", c->rgb[3 * i], c->rgb[3 * i + 1],
                    c->rgb[3 * i + 2]);
        fputc('\n', f);
    }
    fclose(f);
    return MV_OK;
}

void mv_cloud_free(mv_cloud *c)
{
    free(c->xyz);
    free(c->rgb);
    memset(c, 0, sizeof(*c));
}
