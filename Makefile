CC      ?= cc
CFLAGS  ?= -std=c99 -pedantic -Wall -Wextra -O2
CPPFLAGS += -Iinclude
LDLIBS   = -lm

SRC = src/mat.c src/cam.c src/epipolar.c src/triangulate.c src/rectify.c \
      src/stereo.c src/img.c src/cloud.c
OBJ = $(SRC:.c=.o)

all: libmv.a demo_synthetic test_mv

libmv.a: $(OBJ)
	ar rcs $@ $(OBJ)

demo_synthetic: demo/synthetic.o libmv.a
	$(CC) $(CFLAGS) -o $@ demo/synthetic.o libmv.a $(LDLIBS)

test_mv: tests/test_mv.o libmv.a
	$(CC) $(CFLAGS) -o $@ tests/test_mv.o libmv.a $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

check: test_mv
	./test_mv

demo: demo_synthetic
	./demo_synthetic

doc: doc/multiview.pdf

doc/multiview.pdf: doc/multiview.tex doc/refs.bib
	cd doc && pdflatex -interaction=nonstopmode multiview >/dev/null \
	&& bibtex multiview >/dev/null \
	&& pdflatex -interaction=nonstopmode multiview >/dev/null \
	&& pdflatex -interaction=nonstopmode multiview >/dev/null

clean:
	rm -f $(OBJ) demo/synthetic.o tests/test_mv.o libmv.a \
	      demo_synthetic test_mv out_cloud.ply \
	      doc/*.aux doc/*.log doc/*.bbl doc/*.blg doc/*.out doc/*.toc

.PHONY: all check demo doc clean
