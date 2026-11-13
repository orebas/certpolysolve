#!/usr/bin/env python3
"""Check that the guides quote their example programs as they are.

Every code block of a guide that holds a whole program (one with `int main(`) must appear verbatim in one of the
guide's example files, which the build compiles and the suite runs, and every example file must be quoted in its
guide. Shorter blocks are illustrations and are not checked. An example file may begin with a comment that the
guide leaves out.

- docs/quickstart-cpp.md quotes examples/quickstart/*.cpp (the quickstart_* tests), and docs/examples-cpp.md quotes
  examples/topics/*.cpp (the topic_* tests).
- docs/api.md quotes examples/guide/*.c and *.cpp (guide_first_solve_c and guide_first_solve_cpp).
- README.md quotes examples/readme/first_solve.cpp (readme_first_solve).

A C block of docs/api.md that begins "/* Declared in <certpolysolve/certpolysolve.h>. */" lists prototypes: each must be
a declaration of the header, compared with CPS_API removed and runs of white space made one space.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# The guide, the languages of its fenced blocks, the folder of its examples and their file patterns.
GUIDES = [
    ("docs/quickstart-cpp.md", ("cpp",), "examples/quickstart", ("*.cpp",)),
    ("docs/api.md", ("c", "cpp"), "examples/guide", ("*.c", "*.cpp")),
    ("docs/api-spec-0.5.md", ("c", "cpp"), "examples/guide", ("first_solve.c", "first_solve.cpp")),
    ("README.md", ("cpp",), "examples/readme", ("*.cpp",)),
    ("docs/examples-cpp.md", ("cpp",), "examples/topics", ("*.cpp",)),
]


def check(guide, languages, folder, patterns):
    text = (ROOT / guide).read_text(encoding="utf-8")
    fence = "|".join(languages)
    blocks = [b for b in re.findall(rf"```(?:{fence})\n(.*?)```", text, flags=re.S) if "int main(" in b]
    examples = sorted(path for pattern in patterns for path in (ROOT / folder).glob(pattern))
    sources = {path.name: path.read_text(encoding="utf-8") for path in examples}
    failures = []
    quoted = set()
    for i, block in enumerate(blocks, 1):
        matches = [name for name, source in sources.items() if block in source]
        if not matches:
            first = block.strip().splitlines()[0]
            failures.append(f"{guide}: program {i} (starting {first!r}) is not in any file of {folder}")
        quoted.update(matches)
    for name in sorted(set(sources) - quoted):
        failures.append(f"{folder}/{name} is not quoted in {guide}")
    return len(blocks), failures


PROTOTYPES = "/* Declared in <certpolysolve/certpolysolve.h>. */"


def normalized(text):
    return " ".join(text.replace("CPS_API ", "").split())


def check_prototypes(guide="docs/api.md", header="include/certpolysolve/certpolysolve.h"):
    """Every prototype the guide lists is one of the header's declarations, word for word."""
    source = re.sub(r"/\*.*?\*/", " ", (ROOT / header).read_text(encoding="utf-8"), flags=re.S)
    source = "\n".join(line for line in source.splitlines() if not line.lstrip().startswith("#"))  # no #define, #if
    declarations = {normalized(d) + ";" for d in source.split(";")}
    text = (ROOT / guide).read_text(encoding="utf-8")
    blocks = [b for b in re.findall(r"```c\n(.*?)```", text, flags=re.S) if b.startswith(PROTOTYPES)]
    failures, count = [], 0
    for block in blocks:
        for line in block.splitlines()[1:]:
            if line.strip():
                count += 1
                if normalized(line) not in declarations:
                    failures.append(f"{guide}: {line.strip()!r} is not declared so in {header}")
    return count, failures


def main():
    status = 0
    count, failures = check_prototypes()
    for failure in failures:
        print(failure)
    if failures or not count:
        status = 1
    print(f"docs/api.md: {count} prototypes as the header declares them")
    for guide, languages, folder, patterns in GUIDES:
        count, failures = check(guide, languages, folder, patterns)
        for failure in failures:
            print(failure)
        if failures:
            status = 1
        else:
            print(f"{guide}: {count} programs quoted as they are")
    return status


if __name__ == "__main__":
    sys.exit(main())
