/* libFuzzer harness: mv_pgm_read parses files from disk (recorded
 * frames, replay input). Fuzz the header/parser against arbitrary bytes
 * via an in-memory temp file so a corrupt PGM cannot crash the tools.
 *
 * Build (clang):
 *   clang -std=c99 -g -O1 -fsanitize=fuzzer,address,undefined \
 *     -Iinclude fuzz/fuzz_pgm.c src/all.c -lm -o fuzz_pgm */
#define _POSIX_C_SOURCE 200809L
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "mv/mv.h"

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    char path[] = "/tmp/fuzzpgmXXXXXX";
    int fd = mkstemp(path);
    unsigned char *img = NULL;
    int w = 0, h = 0;
    if (fd < 0)
        return 0;
    if (write(fd, data, size) != (ssize_t)size) {
        close(fd);
        unlink(path);
        return 0;
    }
    close(fd);
    if (mv_pgm_read(path, &img, &w, &h) == MV_OK)
        free(img);
    unlink(path);
    return 0;
}
