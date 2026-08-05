/* libFuzzer harness: the pattern readers eat untrusted camera bytes on
 * every live frame -- fuzz them against arbitrary input so a malformed
 * or adversarial frame can never crash the hub.  The fuzzer synthesizes
 * an image of bounded size from the input and runs both tier readers.
 *
 * Build (clang):
 *   clang -std=c99 -g -O1 -fsanitize=fuzzer,address,undefined \
 *     -Iinclude fuzz/fuzz_reader.c src/all.c -lm -o fuzz_reader
 *   ./fuzz_reader -max_total_time=60 */
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    /* first 4 bytes pick a plausible image geometry; the rest is pixels
     * (tiled to fill), so the fuzzer explores decode paths not just the
     * size guard */
    int w, h;
    unsigned char *img;
    mv_read_result rr;
    size_t i;
    if (size < 8)
        return 0;
    w = 64 + (int)((data[0] | (data[1] << 8)) % 1217);   /* 64..1280 */
    h = 64 + (int)((data[2] | (data[3] << 8)) % 656);    /* 64..720  */
    img = (unsigned char *)malloc((size_t)w * h);
    if (!img)
        return 0;
    for (i = 0; i < (size_t)w * h; i++)
        img[i] = data[4 + (i % (size - 4))];
    mv_read_pattern(&rr, img, w, h);
    mv_read_coarse(&rr, img, w, h);
    /* if a read claimed corners, exercise the downstream pose path the
     * hub would run -- a decode that "succeeds" on garbage must not
     * blow up the geometry */
    if (rr.n >= 4 && rr.n <= MV_READ_MAXC) {
        double obj[2 * MV_READ_MAXC], K[9], kr[2] = { 0, 0 };
        mv_camera cam;
        int j, ok = 1;
        memset(K, 0, sizeof(K));
        K[0] = K[4] = 800.0; K[2] = w / 2.0; K[5] = h / 2.0; K[8] = 1.0;
        for (j = 0; j < rr.n; j++) {
            int id = rr.id[j];
            if (id < 0 || id >= MV_PAT_CORNER_COLS * MV_PAT_CORNER_ROWS) {
                ok = 0;
                break;
            }
            {
                double xy[2];
                mv_pattern_corner_px(id % MV_PAT_CORNER_COLS,
                                     id / MV_PAT_CORNER_COLS, xy);
                obj[2 * j] = xy[0] * 2.7e-4;
                obj[2 * j + 1] = xy[1] * 2.7e-4;
            }
        }
        if (ok)
            mv_calib_plane_pose(&cam, K, obj, rr.uv, rr.n);
    }
    free(img);
    return 0;
}
