#!/usr/bin/env python3
"""Build-consistency check: every make target named by CI or the deploy
scripts must exist in the Makefile, and every suite in `make check`
must also run on Linux CI.

Why this exists: the hub's make target was renamed livehub -> hubengine
and the CI workflow, deploy/build-linux.sh, stress-live.sh and
valgrind-hub.sh kept calling the old name.  Every Linux job failed with
"No rule to make target 'livehub'" for days -- the sanitizer, valgrind
and scan-build gates were silently not running, which is worse than a
loud failure.  A parallel build then added three test suites to
`make check` that Linux CI never learned about, so they were never
exercised under ASan/valgrind either.  Both are name drift between
files that no compiler checks; this test is the check.
"""

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read(rel):
    with open(os.path.join(ROOT, rel), encoding="utf-8") as f:
        return f.read()


def strip_prose(text):
    """Drop comment lines and CI template expressions before scanning:
    both contain English words that read like make targets."""
    text = re.sub(r"\$\{\{.*?\}\}", "$X", text, flags=re.S)
    out = []
    for line in text.splitlines():
        s = line.lstrip()
        # comments, and YAML step names (human prose, never commands --
        # "make targets ... resolve" is a title, not an invocation)
        if s.startswith("#") or re.match(r"^-?\s*name:\s", s):
            continue
        out.append(line.split(" #")[0])
    return "\n".join(out)


def makefile_targets(text):
    """Explicit target names (pattern rules and variables excluded)."""
    out = set()
    for line in text.splitlines():
        if line.startswith("\t") or line.lstrip().startswith("#"):
            continue
        m = re.match(r"^(\$\(OUT\)/)?([A-Za-z0-9_./+-]+)\s*:(?!=)", line)
        if m:
            # Binaries are built into $(OUT) (CLAUDE.md T2), so a rule reads
            # "$(OUT)/hubengine:". The scripts still invoke ./hubengine, and
            # this check is about a target existing, not about where it lands.
            out.add(m.group(2))
            if m.group(1):
                out.add(m.group(1) + m.group(2))
    return out


def is_target_word(w):
    """A word in a make invocation that names a target, not a flag,
    a variable assignment, or a shell expansion."""
    if not w or w.startswith("-") or "=" in w or "$" in w:
        return False
    return re.match(r"^[A-Za-z0-9_./+-]+$", w) is not None


def invoked_targets(text):
    """Targets named on `make ...` command lines."""
    out = set()
    for m in re.finditer(r"\bmake\s+([^\n|&;]*)", text):
        for w in m.group(1).split():
            if is_target_word(w):
                out.add(w)
    return out


def for_loop_targets(text):
    """Targets listed in `for t in ... ; do ... make "$t"` loops --
    build-linux.sh builds one target per stage this way."""
    out = set()
    for m in re.finditer(r"for\s+\w+\s+in\s+(.*?);?\s*do", text,
                         re.DOTALL):
        body = m.group(1).replace("\\\n", " ")
        for w in body.split():
            if is_target_word(w):
                out.add(w)
    return out


def check_suite_binaries(text):
    """Binaries the `check` target runs, i.e. ./test_foo lines."""
    out = set()
    m = re.search(r"^check:.*?(?=^\S)", text, re.M | re.S)
    if m:
        for line in m.group(0).splitlines():
            b = re.match(r"^\t\./(\S+)", line)
            if b:
                out.add(b.group(1))
    return out


def main():
    fails = []
    mk = read("Makefile")
    targets = makefile_targets(mk)

    # 1. every referenced make target exists
    sources = {
        ".github/workflows/ci.yml": invoked_targets,
        "deploy/build-linux.sh": lambda t: (invoked_targets(t)
                                            | for_loop_targets(t)),
        "deploy/stress-live.sh": invoked_targets,
        "deploy/valgrind-hub.sh": invoked_targets,
    }
    for rel, extract in sources.items():
        text = strip_prose(read(rel))
        for t in sorted(extract(text)):
            if t not in targets:
                fails.append("%s: make target %r does not exist in the "
                             "Makefile" % (rel, t))

    # 2. binaries the scripts execute directly must be built by a target
    #    of the same name (./hubengine, ./replaycam, ...)
    for rel in ("deploy/stress-live.sh", "deploy/valgrind-hub.sh"):
        text = strip_prose(read(rel))
        for m in re.finditer(r"\./([A-Za-z0-9_]+)\s", text):
            b = m.group(1)
            if b.startswith("test_") or b in targets:
                continue
            if b in ("configure",):
                continue
            fails.append("%s: runs ./%s but no Makefile target builds "
                         "it" % (rel, b))

    # 3. every suite in `make check` also runs on Linux CI
    suites = check_suite_binaries(mk)
    linux = for_loop_targets(strip_prose(read("deploy/build-linux.sh")))
    for s in sorted(suites - linux):
        fails.append("deploy/build-linux.sh: test suite %r is in "
                     "`make check` but never built/run on Linux CI" % s)

    # 4. suites must be gitignored, so a build never dirties the tree
    ignored = set(read(".gitignore").split())
    for s in sorted(suites - ignored):
        fails.append(".gitignore: test binary %r is not ignored" % s)

    if fails:
        print("BUILD CONSISTENCY CHECK FAILED:")
        for f in fails:
            print("  - " + f)
        return 1
    print("build consistency ok: %d make targets, %d suites, all "
          "references resolve" % (len(targets), len(suites)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
