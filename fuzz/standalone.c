/* Portable fuzz driver: runs LLVMFuzzerTestOneInput on deterministic
 * pseudo-random inputs so the harnesses run anywhere (no libFuzzer
 * runtime needed -- macOS Xcode clang lacks it). CI runs the real
 * coverage-guided libFuzzer; this is the everywhere-reproducible net.
 *   cc -O1 -g -fsanitize=address,undefined -Iinclude fuzz/standalone.c \
 *      fuzz/<harness>.c src/*.c -lm -o drv && ./drv <iters> <seed>  */
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size);

int main(int argc, char **argv)
{
    unsigned long long s = argc > 2 ? strtoull(argv[2], NULL, 10)
                                    : 20260807ULL;
    long iters = argc > 1 ? atol(argv[1]) : 200000;
    long i;
    unsigned char buf[8192];
    for (i = 0; i < iters; i++) {
        size_t n, j;
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        n = (size_t)((s >> 33) % sizeof(buf));
        for (j = 0; j < n; j++) {
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            buf[j] = (unsigned char)(s >> 40);
        }
        LLVMFuzzerTestOneInput(buf, n);
        if ((i & 0x3fff) == 0) {
            printf("\r%ld/%ld", i, iters);
            fflush(stdout);
        }
    }
    printf("\rstandalone fuzz done: %ld iterations, no crash\n", iters);
    return 0;
}
