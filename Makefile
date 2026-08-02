CC      ?= cc
CFLAGS  ?= -std=c99 -pedantic -Wall -Wextra -O2
CPPFLAGS += -Iinclude
LDLIBS   = -lm

SRC = src/mat.c src/cam.c src/epipolar.c src/triangulate.c src/rectify.c \
      src/stereo.c src/img.c src/cloud.c src/target.c src/graycode.c \
      src/calib.c
OBJ = $(SRC:.c=.o)

all: libmv.a demo_synthetic demo_calibrate test_mv

libmv.a: $(OBJ)
	ar rcs $@ $(OBJ)

demo_synthetic: demo/synthetic.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/synthetic.o libmv.a $(LDLIBS)

demo_calibrate: demo/calibrate.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/calibrate.o libmv.a $(LDLIBS)

test_mv: tests/test_mv.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_mv.o libmv.a $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

check: test_mv
	./test_mv

demo: demo_synthetic demo_calibrate
	./demo_synthetic
	./demo_calibrate

doc: doc/multiview.pdf

doc/multiview.pdf: doc/multiview.tex doc/refs.bib
	cd doc && pdflatex -interaction=nonstopmode multiview >/dev/null \
	&& bibtex multiview >/dev/null \
	&& pdflatex -interaction=nonstopmode multiview >/dev/null \
	&& pdflatex -interaction=nonstopmode multiview >/dev/null

clean:
	rm -f $(OBJ) demo/synthetic.o demo/calibrate.o tests/test_mv.o libmv.a \
	      demo_synthetic demo_calibrate test_mv out_cloud.ply \
	      target_letter.pgm \
	      doc/*.aux doc/*.log doc/*.bbl doc/*.blg doc/*.out doc/*.toc

.PHONY: all check demo doc clean
