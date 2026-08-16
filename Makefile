# Build output goes to the disposable tree (CLAUDE.md T2/T4); the pdf stays
# in the repo (T9). BUILD_TARGET_PREFIX comes from ~/.zshenv.
BUILD_TARGET_PREFIX ?= /Users/thv/Claude/Projects
OUT = $(BUILD_TARGET_PREFIX)/_buildoutput/3d-multiview

CC      ?= cc
CFLAGS  ?= -std=c99 -pedantic -Wall -Wextra -O2
CPPFLAGS += -Iinclude
LDLIBS   = -lm -lpthread

SRC = src/mat.c src/cam.c src/epipolar.c src/triangulate.c src/rectify.c \
      src/stereo.c src/img.c src/cloud.c src/target.c src/graycode.c \
      src/calib.c src/pattern.c src/render.c src/reader.c src/tsdf.c \
      src/refine.c src/optimal.c src/feat.c src/session.c src/photo.c src/bundle.c src/rot.c \
      src/plane.c src/sync.c src/track.c
OBJ = $(SRC:%.c=$(OUT)/%.o)

all: $(OUT)/libmv.a $(OUT)/demo_synthetic $(OUT)/demo_calibrate $(OUT)/demo_track $(OUT)/demo_diagnose $(OUT)/demo_insects $(OUT)/demo_lightlog $(OUT)/demo_people $(OUT)/demo_patternsim $(OUT)/demo_tsdf $(OUT)/demo_room $(OUT)/test_mv

$(OUT)/libmv.a: $(OBJ)
	ar rcs $@ $(OBJ)

$(OUT)/demo_synthetic: $(OUT)/demo/synthetic.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/synthetic.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/demo_calibrate: $(OUT)/demo/calibrate.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/calibrate.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/demo_tsdf: $(OUT)/demo/tsdfsim.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/tsdfsim.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/demo_room: $(OUT)/demo/roomsim.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/roomsim.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/demo_slight: $(OUT)/demo/slightsim.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/slightsim.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/calibreal: $(OUT)/tools/calibreal.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/calibreal.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/scenecloud: $(OUT)/tools/scenecloud.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/scenecloud.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/rigcalib: $(OUT)/tools/rigcalib.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/rigcalib.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/annotate: $(OUT)/tools/annotate.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/annotate.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/slreal: $(OUT)/tools/slreal.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/slreal.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/replaycam: $(OUT)/tools/replaycam.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/replaycam.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/hubengine: $(OUT)/tools/livehub.o $(OUT)/tools/hub_clock.o $(OUT)/tools/hub_solve.o $(OUT)/tools/hub_pair.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/livehub.o $(OUT)/tools/hub_clock.o \
	  $(OUT)/tools/hub_solve.o $(OUT)/tools/hub_pair.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/nettest: $(OUT)/tools/nettest.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/nettest.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/genframes: $(OUT)/tools/genframes.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/genframes.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/densereal: $(OUT)/tools/densereal.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/densereal.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/stream_cam_v4l2: $(OUT)/tools/stream_cam_v4l2.o
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tools/stream_cam_v4l2.o $(LDLIBS)

$(OUT)/demo_patternsim: $(OUT)/demo/patternsim.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/patternsim.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/demo_insects: $(OUT)/demo/track_insects.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/track_insects.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/demo_lightlog: $(OUT)/demo/lightlog.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/lightlog.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/demo_people: $(OUT)/demo/track_people.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/track_people.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/demo_diagnose: $(OUT)/demo/diagnose.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/diagnose.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/demo_track: $(OUT)/demo/track_robot.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/demo/track_robot.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_mv: $(OUT)/tests/test_mv.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_mv.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_refine: $(OUT)/tests/test_refine.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_refine.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_optimal: $(OUT)/tests/test_optimal.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_optimal.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_feat: $(OUT)/tests/test_feat.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_feat.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_session: $(OUT)/tests/test_session.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_session.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_photo: $(OUT)/tests/test_photo.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_photo.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_bundle: $(OUT)/tests/test_bundle.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_bundle.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_mux: $(OUT)/tests/test_mux.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_mux.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_reader_speed: $(OUT)/tests/test_reader_speed.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_reader_speed.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_clock_sync: $(OUT)/tests/test_clock_sync.o $(OUT)/tools/hub_clock.o
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_clock_sync.o $(OUT)/tools/hub_clock.o $(LDLIBS)

$(OUT)/test_hub_solve: $(OUT)/tests/test_hub_solve.o $(OUT)/tools/hub_solve.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_hub_solve.o $(OUT)/tools/hub_solve.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_hub_pair: $(OUT)/tests/test_hub_pair.o $(OUT)/tools/hub_pair.o
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_hub_pair.o $(OUT)/tools/hub_pair.o $(LDLIBS)

$(OUT)/test_tsdffast: $(OUT)/tests/test_tsdffast.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_tsdffast.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_plane: $(OUT)/tests/test_plane.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_plane.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_sync: $(OUT)/tests/test_sync.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_sync.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/test_track: $(OUT)/tests/test_track.o $(OUT)/libmv.a
	@mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $(OUT)/tests/test_track.o $(OUT)/libmv.a $(LDLIBS)

$(OUT)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

check: $(OUT)/test_mv $(OUT)/test_refine $(OUT)/test_optimal $(OUT)/test_feat $(OUT)/test_session $(OUT)/test_photo $(OUT)/test_bundle $(OUT)/test_mux $(OUT)/test_reader_speed $(OUT)/test_clock_sync $(OUT)/test_hub_solve $(OUT)/test_hub_pair $(OUT)/test_tsdffast $(OUT)/test_plane $(OUT)/test_sync $(OUT)/test_track
	$(OUT)/test_mv
	$(OUT)/test_refine
	$(OUT)/test_optimal
	$(OUT)/test_feat
	$(OUT)/test_session
	$(OUT)/test_photo
	$(OUT)/test_bundle
	$(OUT)/test_mux
	$(OUT)/test_reader_speed
	$(OUT)/test_clock_sync
	$(OUT)/test_hub_solve
	$(OUT)/test_hub_pair
	$(OUT)/test_tsdffast
	$(OUT)/test_plane
	$(OUT)/test_sync
	$(OUT)/test_track
	python3 tests/check_targets.py
	$(MAKE) $(OUT)/doc/multiview.aux
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

demo: $(OUT)/demo_synthetic $(OUT)/demo_calibrate $(OUT)/demo_track $(OUT)/demo_diagnose \
      $(OUT)/demo_insects $(OUT)/demo_lightlog $(OUT)/demo_people $(OUT)/demo_patternsim $(OUT)/demo_tsdf \
      $(OUT)/demo_room
	$(OUT)/demo_synthetic
	$(OUT)/demo_calibrate
	$(OUT)/demo_track
	$(OUT)/demo_diagnose
	$(OUT)/demo_insects
	$(OUT)/demo_lightlog
	$(OUT)/demo_people
	$(OUT)/demo_patternsim
	$(OUT)/demo_tsdf
	$(OUT)/demo_room

doc: doc/multiview.pdf $(OUT)/doc/multiview.aux

# tests/check_bib.py reads doc/multiview.aux, a latex byproduct that
# `make clean` removes while the versioned pdf stays. The pdf is then
# newer than its sources, so the rule below never fires and `make doc`
# cannot restore the aux -- which made `make check` fail after a clean
# with advice ("run make doc first") that did not work. Keying a rule
# on the aux itself forces the one rebuild that regenerates it.
$(OUT)/doc/multiview.aux:
	$(MAKE) -B doc/multiview.pdf

# bibtex refuses to write outside the current directory (TeX Live ships
# openout_any = p), so the run happens inside $(OUT)/doc with the sources found
# through TEXINPUTS/BIBINPUTS. The pdf is copied back into the repo (T9).
doc/multiview.pdf: doc/multiview.tex doc/refs.bib
	@mkdir -p $(OUT)/doc
	cd $(OUT)/doc \
	&& TEXINPUTS=$(CURDIR)/doc//: BIBINPUTS=$(CURDIR)/doc//: \
	   pdflatex -interaction=nonstopmode $(CURDIR)/doc/multiview.tex >/dev/null \
	&& BIBINPUTS=$(CURDIR)/doc//: bibtex multiview >/dev/null \
	&& TEXINPUTS=$(CURDIR)/doc//: BIBINPUTS=$(CURDIR)/doc//: \
	   pdflatex -interaction=nonstopmode $(CURDIR)/doc/multiview.tex >/dev/null \
	&& TEXINPUTS=$(CURDIR)/doc//: BIBINPUTS=$(CURDIR)/doc//: \
	   pdflatex -interaction=nonstopmode $(CURDIR)/doc/multiview.tex >/dev/null
	cp $(OUT)/doc/multiview.pdf doc/multiview.pdf

clean:
	rm -rf $(OUT)

.PHONY: all check demo doc clean
