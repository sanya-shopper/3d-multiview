#include <stdio.h>
#include <stdlib.h>

#include "mv/core.h"
#include "mv/img.h"

/* Read the next integer token, skipping whitespace and '#' comments. */
static int next_int(FILE *f, int *out)
{
    int c, v = 0, got = 0;
    for (;;) {
        c = fgetc(f);
        if (c == '#') {
            while (c != '\n' && c != EOF)
                c = fgetc(f);
            continue;
        }
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            if (got)
                break;
            continue;
        }
        if (c < '0' || c > '9') {
            if (got)
                break;
            return MV_ERR;
        }
        if (v > 100000000) /* cap before signed overflow (UB) on a
                            * pathological digit run; real dims are small */
            return MV_ERR;
        v = v * 10 + (c - '0');
        got = 1;
    }
    *out = v;
    return MV_OK;
}

int mv_pgm_read(const char *path, unsigned char **data, int *w, int *h)
{
    FILE *f = fopen(path, "rb");
    int maxval;
    size_t sz;
    if (!f)
        return MV_ERR;
    if (fgetc(f) != 'P' || fgetc(f) != '5')
        goto fail;
    if (next_int(f, w) != MV_OK || next_int(f, h) != MV_OK
        || next_int(f, &maxval) != MV_OK)
        goto fail;
    if (*w <= 0 || *h <= 0 || maxval <= 0 || maxval > 255)
        goto fail;
    /* bound total pixels before allocating: a crafted/corrupt header
     * (e.g. 66350x66350) would otherwise malloc gigabytes -- a DoS on
     * any tool reading files. 100 Mpx covers any real camera frame. */
    if ((double)(*w) * (double)(*h) > 100.0e6)
        goto fail;
    sz = (size_t)(*w) * (size_t)(*h);
    *data = (unsigned char *)malloc(sz);
    if (!*data)
        goto fail;
    if (fread(*data, 1, sz, f) != sz) {
        free(*data);
        *data = NULL;
        goto fail;
    }
    fclose(f);
    return MV_OK;
fail:
    fclose(f);
    return MV_ERR;
}

int mv_pgm_write(const char *path, const unsigned char *data, int w, int h)
{
    FILE *f = fopen(path, "wb");
    size_t sz = (size_t)w * (size_t)h;
    if (!f)
        return MV_ERR;
    fprintf(f, "P5\n%d %d\n255\n", w, h);
    if (fwrite(data, 1, sz, f) != sz) {
        fclose(f);
        return MV_ERR;
    }
    fclose(f);
    return MV_OK;
}
