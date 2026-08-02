#include <stdlib.h>

#include "mv/core.h"
#include "mv/img.h"
#include "mv/target.h"

int mv_target_checkerboard(double *obj, int inner_cols, int inner_rows,
                           double square)
{
    int r, c;
    for (r = 0; r < inner_rows; r++) {
        for (c = 0; c < inner_cols; c++) {
            obj[2 * (r * inner_cols + c)] = c * square;
            obj[2 * (r * inner_cols + c) + 1] = r * square;
        }
    }
    return inner_cols * inner_rows;
}

int mv_target_render_pgm(const char *path, int inner_cols, int inner_rows,
                         int square_px, int margin_px)
{
    int sq_cols = inner_cols + 1, sq_rows = inner_rows + 1;
    int w = sq_cols * square_px + 2 * margin_px;
    int h = sq_rows * square_px + 2 * margin_px;
    unsigned char *img;
    int x, y, ret;

    if (inner_cols < 2 || inner_rows < 2 || square_px < 1 || margin_px < 0)
        return MV_ERR;
    img = (unsigned char *)malloc((size_t)w * h);
    if (!img)
        return MV_ERR;
    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            int bx = x - margin_px, by = y - margin_px;
            unsigned char v = 255;
            if (bx >= 0 && by >= 0 && bx < sq_cols * square_px
                && by < sq_rows * square_px)
                v = ((bx / square_px + by / square_px) % 2) ? 255 : 0;
            img[y * w + x] = v;
        }
    }
    ret = mv_pgm_write(path, img, w, h);
    free(img);
    return ret;
}
