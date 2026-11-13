#!/usr/bin/env python3
"""data/time-by-degree.csv from a gate's timed run (benchmarks/gate/check.sh's folder): one line per Auto solve of
the test set at the default 53 bits that met its goal, as the 24 October data was drawn: the candidate's all-roots,
accuracy and real-roots jobs (ids cand-auto--, cand-auto-accuracy--, cand-real--), the degree from the answer, the
process's wall-clock seconds from the run's record. Only the answer's first bytes are read (the largest answers are
hundreds of megabytes). Prints the table of docs/how-it-works.md ("How fast"): per band of degree, the median and the
ninetieth percentile.

    time_by_degree.py CHECK_DIR [OUT_CSV]     (default OUT_CSV: data/time-by-degree.csv beside this script)
"""
import csv
import json
import re
import statistics
import sys
from pathlib import Path

QUERIES = {"cand-auto--": "all", "cand-auto-accuracy--": "accuracy", "cand-real--": "real"}
BANDS = [(0, 256, "below 256"), (256, 1024, "256 to 1,023"), (1024, 4096, "1,024 to 4,095"),
         (4096, 16384, "4,096 to 16,383"), (16384, None, "16,384 and up")]
HEAD_BYTES = 4096


def query_of(job_id):
    for prefix, query in QUERIES.items():
        if job_id.startswith(prefix):
            return query
    return None


def answer_head(path):
    """degree, complete and goal_met from the answer's first bytes, or None."""
    try:
        with open(path, "rb") as f:
            head = f.read(HEAD_BYTES).decode("utf-8", "replace")
    except OSError:
        return None
    degree = re.search(r'"degree":(\d+)', head)
    complete = re.search(r'"complete":(true|false)', head)
    met = re.search(r'"goal_met":(true|false)', head)
    if not (degree and complete and met):
        return None
    return int(degree.group(1)), complete.group(1) == "true", met.group(1) == "true"


def percentile(values, share):
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(share * len(ordered)))]


def main():
    check = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parent / "data" / "time-by-degree.csv"
    jobs = {j["id"]: j for j in json.load(open(check / "jobs.json"))}
    rows, skipped = [], 0
    for line in open(check / "run" / "results.jsonl"):
        record = json.loads(line)
        query = query_of(record["id"])
        job = jobs.get(record["id"])
        if query is None or job is None or job["command"][3] != "53":
            continue
        head = answer_head(check / "run" / record["id"] / "stdout.txt") if record.get("returncode") == 0 else None
        if head is None or not (head[1] and head[2]):
            skipped += 1
            continue
        rows.append((head[0], round(record["wall_seconds"], 4), query))
    rows.sort()
    with open(out, "w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["degree", "seconds", "query"])
        writer.writerows(rows)
    print(f"{len(rows)} solves written to {out} ({skipped} not finished or unreadable)")
    print("| Degree | Typical solve | Nine in ten under |")
    print("|---|---:|---:|")
    for low, high, name in BANDS:
        times = [s for d, s, _ in rows if d >= low and (high is None or d < high)]
        if times:
            print(f"| {name} | {statistics.median(times):.3g} s | {percentile(times, 0.9):.3g} s |   ({len(times)} solves)")


if __name__ == "__main__":
    main()
