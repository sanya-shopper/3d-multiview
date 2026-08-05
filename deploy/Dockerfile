# Ubuntu validation image for the multiview repo.
#
# Purpose: reproduce the owner's ultimate deployment target (Ubuntu) so
# Linux portability problems surface before any real Ubuntu machine
# exists.  Build with the REPO ROOT as context:
#
#     docker build -f deploy/Dockerfile -t multiview-linux-validate .
#     docker run --rm multiview-linux-validate
#
# (deploy/validate.sh does both and summarizes the result.)
#
# python3 is installed for completeness (tests/check_bib.py), but the
# in-container script deliberately skips the bib check: it is a doc
# check, not a Linux-portability check, and the doc toolchain (LaTeX)
# is not part of this image.

FROM ubuntu:24.04

RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential python3 \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /multiview
COPY . /multiview

# Start from a pristine tree even if the host context contained
# build artifacts, then hand off to the staged build/test script.
RUN make clean >/dev/null 2>&1 || true

CMD ["/bin/sh", "deploy/build-linux.sh"]
