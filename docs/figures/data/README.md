# The figures' data

Every file here is the solver's own output, saved as it came; `../make_figures.py` draws the figures from them and
nothing else. To regenerate after a change to the solver, run the commands below from the repository root with a
built `certpolysolve_cli` and run the script again.

- `degree24-coefficients.txt`, `degree24-roots.txt`, `degree24-echo.txt`: a polynomial of degree 24 with small
  integer coefficients (constant term first in the file), solved with the defaults:
  `build/dev/certpolysolve_cli $(cat docs/figures/data/degree24-coefficients.txt) > degree24-roots.txt 2> degree24-echo.txt`.
  (Asked for boxes of radius at most 0.05 instead, the solver still proves boxes of 10⁻¹³ and less, because a tight
  box is cheap here; so the figure shows the roots as dots.)
- `pair-coefficients.txt`, `pair-clusters.txt`, `pair-accuracy.txt`: (x − 1)(x − 1 − 2⁻⁶⁰)(x + 1), whose
  coefficients are exact in the `m@e` form (m · 2^e), solved twice:
  `build/dev/certpolysolve_cli --clusters 1e-6 $(cat docs/figures/data/pair-coefficients.txt) > pair-clusters.txt 2>&1`
  (boxes of radius at most 10⁻⁶, each with its count) and
  `build/dev/certpolysolve_cli $(cat docs/figures/data/pair-coefficients.txt) > pair-accuracy.txt 2>&1` (the default:
  each root alone).
- `time-by-degree.csv`: one line per solve of the test set's Auto queries that finished, from a timed run of the
  whole set on one core each (degree, wall-clock seconds, the kind of query). The columns come from the run's
  records: `python3 docs/figures/time_by_degree.py CHECK_DIR` (a timed run's check folder, `
  gate-records/checks/NAME`) writes the file and prints the table of how-it-works.md; the run itself is not
  reproducible from this folder.

- `compare-small-dense.csv`: from the comparison of CertPolySolve with MPSolve, numpy.roots and PolynomialRoots.jl
  (`benchmarks/headline/`; the record is `dev/docs/headline-benchmark-20261025.md`). Per degree and solver, over
  five dense polynomials with coefficients drawn from [-1, 1], the median time of one solve call on one core (the
  median of three in-process repetitions), how many of the five each solver answered, and the median accuracy its
  roots reached in bits (53 for the solvers that prove it). Made by
  `python3 benchmarks/headline/report.py CAMPAIGN --readme-csv docs/figures/data/compare-small-dense.csv`.
