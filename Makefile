CC      ?= cc
CFLAGS  ?= -std=c99 -pedantic -Wall -Wextra -O2
CPPFLAGS += -Iinclude
LDLIBS   = -lm

SRC = src/mat.c src/cam.c src/epipolar.c src/triangulate.c src/rectify.c \
      src/stereo.c src/img.c src/cloud.c src/target.c src/graycode.c \
      src/calib.c src/pattern.c src/render.c src/reader.c src/tsdf.c \
      src/refine.c src/optimal.c src/feat.c src/session.c src/photo.c src/bundle.c
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

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

check: test_mv test_refine test_optimal test_feat test_session test_photo test_bundle test_mux
	./test_mv
	./test_refine
	./test_optimal
	./test_feat
	./test_session
	./test_photo
	./test_bundle
	./test_mux
	python3 tests/check_bib.py

checkbib:
	python3 tests/check_bib.py

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

doc: doc/multiview.pdf

doc/multiview.pdf: doc/multiview.tex doc/refs.bib
	cd doc && pdflatex -interaction=nonstopmode multiview >/dev/null \
	&& bibtex multiview >/dev/null \
	&& pdflatex -interaction=nonstopmode multiview >/dev/null \
	&& pdflatex -interaction=nonstopmode multiview >/dev/null

clean:
	rm -f $(OBJ) demo/synthetic.o demo/calibrate.o demo/track_robot.o demo/diagnose.o demo/track_insects.o demo/lightlog.o demo/track_people.o demo/patternsim.o demo/tsdfsim.o demo/roomsim.o demo/slightsim.o tools/calibreal.o tools/scenecloud.o tools/rigcalib.o tools/annotate.o tools/slreal.o tests/test_mv.o tests/test_refine.o tests/test_optimal.o tests/test_feat.o tests/test_session.o tests/test_photo.o tests/test_bundle.o tests/test_mux.o libmv.a \
	      demo_synthetic demo_calibrate demo_track demo_diagnose demo_insects demo_lightlog demo_people demo_patternsim demo_tsdf demo_room test_mv out_cloud.ply out_track.ply \
	      target_letter.pgm \
	      doc/*.aux doc/*.log doc/*.bbl doc/*.blg doc/*.out doc/*.toc doc/*.brf

.PHONY: all check demo doc clean
