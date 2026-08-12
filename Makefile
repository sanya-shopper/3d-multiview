CC      ?= cc
CFLAGS  ?= -std=c99 -pedantic -Wall -Wextra -O2
CPPFLAGS += -Iinclude
LDLIBS   = -lm -lpthread

SRC = src/mat.c src/cam.c src/epipolar.c src/triangulate.c src/rectify.c \
      src/stereo.c src/img.c src/cloud.c src/target.c src/graycode.c \
      src/calib.c src/pattern.c src/render.c src/reader.c src/tsdf.c \
      src/refine.c src/optimal.c src/feat.c src/session.c src/photo.c src/bundle.c src/rot.c
OBJ = $(SRC:.c=.o)

all: libmv.a demo_synthetic demo_calibrate demo_track demo_diagnose demo_insects demo_lightlog demo_people demo_patternsim demo_tsdf demo_room test_mv

libmv.a: $(OBJ)
	ar rcs $@ $(OBJ)

demo_synthetic: demo/synthetic.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/synthetic.o libmv.a $(LDLIBS)

demo_calibrate: demo/calibrate.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/calibrate.o libmv.a $(LDLIBS)

demo_tsdf: demo/tsdfsim.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/tsdfsim.o libmv.a $(LDLIBS)

demo_room: demo/roomsim.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/roomsim.o libmv.a $(LDLIBS)

demo_slight: demo/slightsim.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/slightsim.o libmv.a $(LDLIBS)

calibreal: tools/calibreal.o libmv.a
	$(CC) $(CFLAGS) -o $@ tools/calibreal.o libmv.a $(LDLIBS)

scenecloud: tools/scenecloud.o libmv.a
	$(CC) $(CFLAGS) -o $@ tools/scenecloud.o libmv.a $(LDLIBS)

rigcalib: tools/rigcalib.o libmv.a
	$(CC) $(CFLAGS) -o $@ tools/rigcalib.o libmv.a $(LDLIBS)

annotate: tools/annotate.o libmv.a
	$(CC) $(CFLAGS) -o $@ tools/annotate.o libmv.a $(LDLIBS)

slreal: tools/slreal.o libmv.a
	$(CC) $(CFLAGS) -o $@ tools/slreal.o libmv.a $(LDLIBS)

replaycam: tools/replaycam.o libmv.a
	$(CC) $(CFLAGS) -o $@ tools/replaycam.o libmv.a $(LDLIBS)

hubengine: tools/livehub.o tools/hub_clock.o tools/hub_solve.o tools/hub_pair.o libmv.a
	$(CC) $(CFLAGS) -o $@ tools/livehub.o tools/hub_clock.o \
	  tools/hub_solve.o tools/hub_pair.o libmv.a $(LDLIBS)

nettest: tools/nettest.o libmv.a
	$(CC) $(CFLAGS) -o $@ tools/nettest.o libmv.a $(LDLIBS)

genframes: tools/genframes.o libmv.a
	$(CC) $(CFLAGS) -o $@ tools/genframes.o libmv.a $(LDLIBS)

stream_cam_v4l2: tools/stream_cam_v4l2.o
	$(CC) $(CFLAGS) -o $@ tools/stream_cam_v4l2.o $(LDLIBS)

demo_patternsim: demo/patternsim.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/patternsim.o libmv.a $(LDLIBS)

demo_insects: demo/track_insects.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/track_insects.o libmv.a $(LDLIBS)

demo_lightlog: demo/lightlog.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/lightlog.o libmv.a $(LDLIBS)

demo_people: demo/track_people.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/track_people.o libmv.a $(LDLIBS)

demo_diagnose: demo/diagnose.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/diagnose.o libmv.a $(LDLIBS)

demo_track: demo/track_robot.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/track_robot.o libmv.a $(LDLIBS)

test_mv: tests/test_mv.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_mv.o libmv.a $(LDLIBS)

test_refine: tests/test_refine.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_refine.o libmv.a $(LDLIBS)

test_optimal: tests/test_optimal.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_optimal.o libmv.a $(LDLIBS)

test_feat: tests/test_feat.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_feat.o libmv.a $(LDLIBS)

test_session: tests/test_session.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_session.o libmv.a $(LDLIBS)

test_photo: tests/test_photo.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_photo.o libmv.a $(LDLIBS)

test_bundle: tests/test_bundle.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_bundle.o libmv.a $(LDLIBS)

test_mux: tests/test_mux.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_mux.o libmv.a $(LDLIBS)

test_reader_speed: tests/test_reader_speed.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_reader_speed.o libmv.a $(LDLIBS)

test_clock_sync: tests/test_clock_sync.o tools/hub_clock.o
	$(CC) $(CFLAGS) -o $@ tests/test_clock_sync.o tools/hub_clock.o $(LDLIBS)

test_hub_solve: tests/test_hub_solve.o tools/hub_solve.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_hub_solve.o tools/hub_solve.o libmv.a $(LDLIBS)

test_hub_pair: tests/test_hub_pair.o tools/hub_pair.o
	$(CC) $(CFLAGS) -o $@ tests/test_hub_pair.o tools/hub_pair.o $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

check: test_mv test_refine test_optimal test_feat test_session test_photo test_bundle test_mux test_reader_speed test_clock_sync test_hub_solve test_hub_pair
	./test_mv
	./test_refine
	./test_optimal
	./test_feat
	./test_session
	./test_photo
	./test_bundle
	./test_mux
	./test_reader_speed
	./test_clock_sync
	./test_hub_solve
	./test_hub_pair
	python3 tests/check_targets.py
	$(MAKE) doc/multiview.aux
	python3 tests/check_bib.py
	$(MAKE) checkweb

checkbib:
	python3 tests/check_bib.py

# The web companion's math layer (web/model.js) is unit-tested headlessly.
# node is present on the GitHub CI runners and dev machines; if it is
# genuinely absent the check reports the skip loudly rather than passing
# silently.
checkweb:
	@if command -v node >/dev/null 2>&1; then \
	  node tests/test_web_model.js && node tests/test_rig_ui.js; \
	else \
	  echo "SKIPPED: web tests (node not installed)"; \
	fi
	python3 tests/check_web_version.py

# Regenerate the visible "page updated" stamp shown in the web footers.
# Run whenever web/ content changes; check_web_version.py enforces it.
webstamp:
	@printf '/* generated by `make webstamp` -- do not edit by hand */\nvar MV_BUILD = { stamp: "%s" };\nif (typeof document !== "undefined") {\n  var _bs = document.getElementById("buildstamp");\n  if (_bs) _bs.textContent = "page updated " + MV_BUILD.stamp;\n}\n' \
	  "$$(TZ=America/Los_Angeles date '+%Y-%m-%d %H:%M %Z')" > web/version.js
	@echo "stamped: $$(grep -o '[0-9][^\"]*' web/version.js | head -1)"

# libFuzzer targets (clang, Linux CI has the runtime; macOS Xcode may not)
FUZZ_SRC = $(filter-out src/%.o,$(SRC))
fuzz_reader: fuzz/fuzz_reader.c $(SRC)
	clang -std=c99 -g -O1 -fsanitize=fuzzer,address,undefined -Iinclude \
	  fuzz/fuzz_reader.c $(SRC) -lm -o $@
fuzz_pgm: fuzz/fuzz_pgm.c $(SRC)
	clang -std=c99 -g -O1 -fsanitize=fuzzer,address,undefined -Iinclude \
	  fuzz/fuzz_pgm.c $(SRC) -lm -o $@

# portable (no-libFuzzer-runtime) standalone drivers for local use
drv_reader: fuzz/standalone.c fuzz/fuzz_reader.c $(SRC)
	$(CC) -std=c99 -O1 -g -fsanitize=address,undefined -Iinclude \
	  fuzz/standalone.c fuzz/fuzz_reader.c $(SRC) -lm -o $@
drv_pgm: fuzz/standalone.c fuzz/fuzz_pgm.c $(SRC)
	$(CC) -std=c99 -O1 -g -fsanitize=address,undefined -Iinclude \
	  fuzz/standalone.c fuzz/fuzz_pgm.c $(SRC) -lm -o $@

demo: demo_synthetic demo_calibrate demo_track demo_diagnose \
      demo_insects demo_lightlog demo_people demo_patternsim demo_tsdf
	./demo_synthetic
	./demo_calibrate
	./demo_track
	./demo_diagnose
	./demo_insects
	./demo_lightlog
	./demo_people
	./demo_patternsim
	./demo_tsdf
	./demo_room

doc: doc/multiview.pdf doc/multiview.aux

# tests/check_bib.py reads doc/multiview.aux, a latex byproduct that
# `make clean` removes while the versioned pdf stays. The pdf is then
# newer than its sources, so the rule below never fires and `make doc`
# cannot restore the aux -- which made `make check` fail after a clean
# with advice ("run make doc first") that did not work. Keying a rule
# on the aux itself forces the one rebuild that regenerates it.
doc/multiview.aux:
	$(MAKE) -B doc/multiview.pdf

doc/multiview.pdf: doc/multiview.tex doc/refs.bib
	cd doc && pdflatex -interaction=nonstopmode multiview >/dev/null \
	&& bibtex multiview >/dev/null \
	&& pdflatex -interaction=nonstopmode multiview >/dev/null \
	&& pdflatex -interaction=nonstopmode multiview >/dev/null

clean:
	rm -f $(OBJ) demo/synthetic.o demo/calibrate.o demo/track_robot.o demo/diagnose.o demo/track_insects.o demo/lightlog.o demo/track_people.o demo/patternsim.o demo/tsdfsim.o demo/roomsim.o demo/slightsim.o tools/calibreal.o tools/scenecloud.o tools/rigcalib.o tools/annotate.o tools/slreal.o tools/replaycam.o tools/livehub.o tools/nettest.o tools/genframes.o tools/stream_cam_v4l2.o tests/test_mv.o tests/test_refine.o tests/test_optimal.o tests/test_feat.o tests/test_session.o tests/test_photo.o tests/test_bundle.o tests/test_mux.o libmv.a \
	      demo_synthetic demo_calibrate demo_track demo_diagnose demo_insects demo_lightlog demo_people demo_patternsim demo_tsdf demo_room test_mv out_cloud.ply out_track.ply \
	      target_letter.pgm \
	      doc/*.aux doc/*.log doc/*.bbl doc/*.blg doc/*.out doc/*.toc doc/*.brf

.PHONY: all check demo doc clean
