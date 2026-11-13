#!/usr/bin/env python3
"""The pages a user reads carry no residue of the project's work (docs/WRITING.md): no dates of the development, no
commit hashes, no gate or workstream numbers, no session names, no quotations of the project's discussions. The
check is a list of patterns over a list of pages; a page joins the list when it has been written for users.
Usage: test_doc_residue.py (from anywhere); exit 1 on the first page with residue."""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# The style sheet (docs/WRITING.md) and llms.txt name the records in order to set them apart, so they are not listed.
PAGES = ["README.md", "docs/quickstart-cpp.md", "docs/examples-cpp.md", "docs/how-it-works.md", "docs/glossary.md", "docs/api-spec-0.5.md",
         "docs/api.md", "docs/mathematics.md", "docs/production-methods.md", "docs/resource-policy.md", "docs/testing.md",
         "examples/README.md"]
MONTHS = "January|February|March|April|May|June|July|August|September|October|November|December"
PATTERNS = [
    (rf"\b\d{{1,2}} (?:{MONTHS})\b", "a date"),
    (rf"\b(?:{MONTHS}) 20\d\d\b", "a date"),
    (r"\b20\d\d-\d\d-\d\d\b", "a date"),
    (r"`(?=[0-9a-f]*[a-f])[0-9a-f]{7,40}`", "a commit hash"),  # a hex word with a letter; 998244353 is a prime
    (r"\bGate \d+\b", "a gate number"),
    (r"\bW\d\b", "a workstream number"),
    (r"(?<!-)\b(?!C99\b|C11\b|C17\b|C23\b)(?:S|C|O)\d{1,2}\b(?! API)", "an item number of the plans"),  # not C99, not -O2
    (r"\bOren\b", "a name from the project's discussions"),
    # The words of the process in prose; `handoff()`, `src/session.cpp` and the like are names in code.
    (r"(?<![`/\w])(?:handoff|HANDOFF|AGENTS\.md|session)\b(?![`(.\w])", "a word of the project's own process"),
    (r"\bbenchmark-results/", "a path into the project's records"),
]


def main():
    status = 0
    for page in PAGES:
        path = ROOT / page
        if not path.exists():
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            for pattern, what in PATTERNS:
                found = re.search(pattern, line)
                if found:
                    print(f"{page}:{number}: {what} ({found.group(0)!r})")
                    status = 1
    if status == 0:
        print(f"{len(PAGES)} user pages carry no residue")
    return status


if __name__ == "__main__":
    sys.exit(main())
