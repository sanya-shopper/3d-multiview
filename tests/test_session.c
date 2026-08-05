#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

static int failures = 0;

#define CHECK(cond, name) \
    do { \
        if (cond) \
            printf("ok:   %s\n", name); \
        else { \
            printf("FAIL: %s\n", name); \
            failures++; \
        } \
    } while (0)

/* frontal synthetic camera over the spec-v1 display (as test_mv's
 * test_pattern): 640x480, f=800, display 0.85 m away, pitch 0.2745 mm */
static void setup_cam(mv_camera *cam)
{
    double center[3] = { MV_PAT_W / 2.0 * 0.0002745,
                         MV_PAT_H / 2.0 * 0.0002745, 0.0 };
    int i;
    mv_cam_set_K(cam, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(cam);
    memset(cam->k, 0, sizeof(cam->k));
    for (i = 0; i < 3; i++)
        cam->t[i] = -center[i];
    cam->t[2] += 0.85;
}

static unsigned char pat[MV_PAT_W * MV_PAT_H];
static unsigned char img[640 * 480];

/* render counter, blind-read, return read result */
static int render_and_read(mv_read_result *rr, unsigned counter)
{
    mv_camera cam;
    unsigned long long seed = 1;
    setup_cam(&cam);
    mv_pattern_render(pat, counter);
    if (mv_render_plane(img, 640, 480, &cam, pat, MV_PAT_W, MV_PAT_H,
                        0.0002745, 128, 0.0, &seed) != MV_OK)
        return MV_ERR;
    return mv_read_pattern(rr, img, 640, 480);
}

static void test_samplers(void)
{
    mv_read_result rr;
    double phs[5];
    unsigned counter;
    int par;

    /* both parities of counter bit 0 (and bit 1 = 0 for both) */
    for (par = 0; par < 2; par++) {
        counter = 123456u + (unsigned)par;
        CHECK(render_and_read(&rr, counter) == MV_OK,
              par ? "samplers: blind read (odd counter)"
                  : "samplers: blind read (even counter)");
        CHECK(rr.counter_valid && rr.counter == counter,
              "samplers: counter decodes");

        CHECK(mv_read_version(&rr, img, 640, 480) == MV_PAT_SPEC_VERSION,
              "samplers: version block decodes to spec version");

        CHECK(mv_read_phases(phs, &rr, img, 640, 480) == MV_OK,
              "samplers: phases sampled");
        {
            /* set bit renders black (g=0); clear bit white (g=1) */
            int i, ok = 1;
            for (i = 0; i < 5; i++) {
                unsigned on = (i < 4) ? (counter & 1u)
                                      : ((counter >> 1) & 1u);
                double ideal = on ? 0.0 : 1.0;
                if (fabs(phs[i] - ideal) > 0.25)
                    ok = 0;
            }
            CHECK(ok, "samplers: all 5 phases within 0.25 of parity");
        }
    }
}

static void test_roundtrip(void)
{
    mv_read_result rr;
    double phs[5];
    FILE *f;
    char line[256];
    mv_session_record rec;
    const double t_mono = 5182.304125;
    const int frame_idx = 1234;
    int nver = 0, nfrm = 0, nctr = 0, ncrn = 0, nphs = 0, nbad = 0;
    int crn_ok = 1, phs_ok = 1, ctr_ok = 1, frm_ok = 1, ver_ok = 1;
    signed char seen[MV_READ_MAXC];

    CHECK(render_and_read(&rr, 123456u) == MV_OK, "roundtrip: read");
    CHECK(mv_read_phases(phs, &rr, img, 640, 480) == MV_OK,
          "roundtrip: phases");

    f = tmpfile();
    CHECK(f != NULL, "roundtrip: tmpfile");
    if (!f)
        return;
    /* deliberately NOT the reference pitch: VER must record the pitch
     * of the display actually used (the live hub once stamped the
     * 0.2745 reference constant while calibrating at 0.1133 -- a 2.4x
     * metric-scale error for any replay) */
    mv_session_ver(f, MV_PAT_SPEC_VERSION, 0.1133, "0");
    mv_session_frm(f, "0", t_mono, frame_idx);
    mv_session_read(f, "0", t_mono, &rr, phs);
    rewind(f);

    memset(seen, 0, sizeof(seen));
    while (fgets(line, (int)sizeof(line), f)) {
        double x, y;
        if (mv_session_parse(&rec, line) != MV_OK) {
            nbad++;
            continue;
        }
        if (strcmp(rec.type, "VER") == 0) {
            nver++;
            if (mv_session_num(&rec, "spec", &x) != MV_OK
                || x != (double)MV_PAT_SPEC_VERSION
                || mv_session_num(&rec, "w", &y) != MV_OK
                || y != (double)MV_PAT_W
                || mv_session_num(&rec, "pitch_mm", &y) != MV_OK
                || fabs(y - 0.1133) > 1e-9
                || !mv_session_field(&rec, "Td_ms"))
                ver_ok = 0;
        } else if (strcmp(rec.type, "FRM") == 0) {
            const char *cam = mv_session_field(&rec, "cam");
            nfrm++;
            if (mv_session_num(&rec, "k", &x) != MV_OK
                || x != (double)frame_idx
                || mv_session_num(&rec, "t", &y) != MV_OK
                || fabs(y - t_mono) > 1e-6
                || !cam || strcmp(cam, "0") != 0)
                frm_ok = 0;
        } else if (strcmp(rec.type, "CTR") == 0) {
            nctr++;
            if (mv_session_num(&rec, "count", &x) != MV_OK
                || x != (double)rr.counter
                || mv_session_num(&rec, "conf", &y) != MV_OK
                || fabs(y - rr.counter_conf) > 1e-4)
                ctr_ok = 0;
        } else if (strcmp(rec.type, "CRN") == 0) {
            double id;
            ncrn++;
            if (mv_session_num(&rec, "id", &id) != MV_OK
                || mv_session_num(&rec, "u", &x) != MV_OK
                || mv_session_num(&rec, "v", &y) != MV_OK) {
                crn_ok = 0;
                continue;
            }
            {
                int i, found = 0;
                for (i = 0; i < rr.n; i++)
                    if (rr.id[i] == (int)id) {
                        found = 1;
                        seen[i] = 1;
                        if (fabs(x - rr.uv[2 * i]) > 1e-3
                            || fabs(y - rr.uv[2 * i + 1]) > 1e-3)
                            crn_ok = 0;
                    }
                if (!found)
                    crn_ok = 0;
            }
        } else if (strcmp(rec.type, "PHS") == 0) {
            double p;
            nphs++;
            if (mv_session_num(&rec, "patch", &p) != MV_OK
                || mv_session_num(&rec, "g", &x) != MV_OK
                || (int)p < 0 || (int)p > 4
                || fabs(x - phs[(int)p]) > 1e-6)
                phs_ok = 0;
        } else {
            nbad++;
        }
    }
    fclose(f);

    CHECK(nbad == 0, "roundtrip: every line parses to a known type");
    CHECK(nver == 1 && ver_ok, "roundtrip: VER fields");
    CHECK(nfrm == 1 && frm_ok, "roundtrip: FRM cam/k/t");
    CHECK((nctr == 1) == (rr.counter_valid != 0) && ctr_ok,
          "roundtrip: CTR presence = validity; count and conf match");
    CHECK(ncrn == rr.n && crn_ok,
          "roundtrip: every corner id and u/v to 1e-3 px");
    {
        int i, all = 1;
        for (i = 0; i < rr.n; i++)
            if (!seen[i])
                all = 0;
        CHECK(all, "roundtrip: no corner lost");
    }
    CHECK(nphs == 5 && phs_ok, "roundtrip: all 5 phases to 1e-6");

    /* the doc's own example lines parse */
    {
        static const char *ex[] = {
            "VER spec=1 w=1920 h=1080 pitch_mm=0.2745 Td_ms=16.667",
            "FRM cam=0 k=1234 t=5182.30412",
            "CTR cam=0 k=1234 count=48213 conf=0.94",
            "PHS cam=0 k=1234 patch=2 g=0.371",
            "CRN cam=0 k=1234 id=57 u=812.431 v=402.117",
            "BND cam=0 k=87 pitch=289.3 dir=down",
            "DSP k=309882 t=5182.29561"
        };
        int i, ok = 1;
        double x;
        for (i = 0; i < 7; i++)
            if (mv_session_parse(&rec, ex[i]) != MV_OK)
                ok = 0;
        /* last one: DSP numeric fields */
        if (ok)
            ok = mv_session_num(&rec, "k", &x) == MV_OK && x == 309882.0
              && mv_session_num(&rec, "t", &x) == MV_OK
              && fabs(x - 5182.29561) < 1e-9;
        /* dir=down is a string field */
        if (ok) {
            const char *d;
            mv_session_parse(&rec, ex[5]);
            d = mv_session_field(&rec, "dir");
            ok = d && strcmp(d, "down") == 0
              && mv_session_num(&rec, "dir", &x) == MV_ERR;
        }
        CHECK(ok, "roundtrip: doc example lines parse");
        CHECK(mv_session_parse(&rec, "no_equals_here foo") == MV_ERR,
              "roundtrip: malformed field rejected");
        CHECK(mv_session_parse(&rec, "# comment\n") == MV_OK
              && rec.type[0] == '\0' && rec.nf == 0,
              "roundtrip: comment line is typeless");
    }
}

/* deterministic gaussian: LCG + Box-Muller */
static unsigned long long lcg_state = 0x243f6a8885a308d3ull;
static double lcg_uniform(void)
{
    lcg_state = lcg_state * 6364136223846793005ull + 1442695040888963407ull;
    return ((double)(lcg_state >> 11) + 0.5) / 9007199254740992.0;
}
static double gauss(void)
{
    double u1 = lcg_uniform(), u2 = lcg_uniform();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}

static void test_clockfit(void)
{
    enum { N = 120 };
    static double t[N], k[N];
    const double rate = 58.8, T_true = 1.0 / 58.8, phi_true = 12.34;
    double phi, T;
    int i;

    /* counters spread over ~60 s of capture */
    for (i = 0; i < N; i++) {
        k[i] = floor(i * 29.4 + 0.5);
        t[i] = phi_true + k[i] * T_true;
    }

    /* exact recovery, no noise */
    CHECK(mv_clock_fit(&phi, &T, t, k, N) == MV_OK,
          "clockfit: fits clean data");
    CHECK(fabs(T - T_true) < 1e-12 && fabs(phi - phi_true) < 1e-12,
          "clockfit: exact (1e-12) at zero noise");

    /* jitter 2 ms + 10 gross +-2 s outliers */
    lcg_state = 0x243f6a8885a308d3ull;
    for (i = 0; i < N; i++)
        t[i] += 0.002 * gauss();
    for (i = 0; i < 10; i++)
        t[7 + i * 11] += (i % 2) ? 2.0 : -2.0;
    CHECK(mv_clock_fit(&phi, &T, t, k, N) == MV_OK,
          "clockfit: fits contaminated data");
    printf("      rate %.6f Hz (true %.6f), offset %.6f s (true %.6f)\n",
           1.0 / T, rate, phi, phi_true);
    CHECK(fabs(1.0 / T - rate) < 0.01,
          "clockfit: rate within 0.01 Hz under outliers");
    CHECK(fabs(phi - phi_true) < 0.005,
          "clockfit: offset within 5 ms under outliers");

    /* degenerate inputs */
    CHECK(mv_clock_fit(&phi, &T, t, k, 1) == MV_ERR,
          "clockfit: rejects n=1");
    k[0] = k[1] = 5.0;
    t[0] = t[1] = 1.0;
    CHECK(mv_clock_fit(&phi, &T, t, k, 2) == MV_ERR,
          "clockfit: rejects constant counter");
}

int main(void)
{
    test_samplers();
    test_roundtrip();
    test_clockfit();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall session tests passed\n");
    return 0;
}
