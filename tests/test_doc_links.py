#!/usr/bin/env python3
"""Check local Markdown links, and that docs/README.md lists every user page and dev/docs/README.md every record.

Relative links must resolve inside the repository; a trailing ":LINE" or
":FIRST-LAST" source-line suffix is allowed. Links into ignored local data
trees (benchmark results, papers, build trees) and absolute links to other
machine-local trees are skipped, because a fresh clone does not contain them.
Verbatim archives under dev/docs/history/ keep their original root-relative
links and are not checked. In a release archive, which leaves out what .gitattributes marks export-ignore (the
records), links into those paths are not checked.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
LOCAL_DATA = ("build",)
SCANNED = ["README.md",
           "src/README.md",
           "examples/README.md"]


def left_out():
    """What a release archive leaves out (the export-ignore lines of .gitattributes), as path prefixes."""
    try:
        lines = (ROOT / ".gitattributes").read_text(encoding="utf-8").splitlines()
    except OSError:
        return ()
    return tuple(line.split()[0] for line in lines
                 if not line.lstrip().startswith("#") and line.split()[-1:] == ["export-ignore"])


def documents():
    for name in SCANNED:
        path = ROOT / name
        if path.exists():
            yield path
    for path in sorted((ROOT / "docs").rglob("*.md")):
        yield path
    for path in sorted((ROOT / "dev" / "docs").rglob("*.md")):
        if "history" not in path.relative_to(ROOT / "dev" / "docs").parts:
            yield path


def targets(path):
    text = path.read_text(encoding="utf-8")
    # Ignore fenced code, where bracketed text is not a link.
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    for match in LINK.finditer(text):
        target = match.group(1)
        if re.match(r"^[a-z][a-z0-9+.-]*:", target) or target.startswith("#"):
            continue
        yield target.split("#", 1)[0]


LEFT_OUT = left_out()


def main():
    missing = []
    checked = 0
    for path in documents():
        for target in targets(path):
            if target.startswith(("/", "~")) and not target.startswith(str(ROOT)):
                continue  # machine-local reference outside the repository
            resolved = (path.parent / target).resolve()
            if not resolved.exists():
                stripped = re.sub(r":\d+(-\d+)?$", "", target)
                resolved = (path.parent / stripped).resolve()
            try:
                relative = resolved.relative_to(ROOT).as_posix()
            except ValueError:
                missing.append(f"{path.relative_to(ROOT)} -> {target} (outside repository)")
                continue
            if relative.startswith(LOCAL_DATA):
                continue
            # In a release archive the records are left out on purpose (the archive is not a
            # target worth a second, unlinked copy of the pages); in the repository they are checked.
            if any(relative.startswith(p) and not (ROOT / p).exists() for p in LEFT_OUT):
                continue
            checked += 1
            if not resolved.exists():
                missing.append(f"{path.relative_to(ROOT)} -> {target}")
    unindexed = []
    for folder in ("docs", "dev/docs"):
        if not (ROOT / folder).is_dir():
            continue  # a release archive carries no records: dev/ is export-ignore in .gitattributes
        index = (ROOT / folder / "README.md").read_text(encoding="utf-8")
        unindexed += [f"{folder}/{p.name}" for p in sorted((ROOT / folder).glob("*.md"))
                      if p.name != "README.md" and f"({p.name}" not in index]
    for line in missing:
        print("missing link:", line)
    for name in unindexed:
        print("not listed in its folder's README.md:", name)
    print(f"{checked} local links checked; {len(missing)} missing; {len(unindexed)} unindexed docs")
    return 1 if missing or unindexed else 0


if __name__ == "__main__":
    sys.exit(main())
