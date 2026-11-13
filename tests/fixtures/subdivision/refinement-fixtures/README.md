# Exact refinement controls

Four square-free Gaussian-integer quadratics with irrational roots:

| File | Polynomial |
| --- | --- |
| `sqrt2.ps` | z²−2 |
| `sqrt_i.ps` | z²−i |
| `translated_sqrt2.ps` | (z−2⁸⁰)²−2 |
| `near_zero.ps` | 2⁸⁰z²−2⁸¹z−1 |

Their nonzero discriminants establish square-freeness. None has zero as a
root. The translated and near-zero examples exercise the scale dependence of
relative accuracy. `test_subdivision_refinement` checks the original exact
polynomial with the rational Sturm–Cauchy rectangle counter after isolation
and after each refinement request. It also checks the exported accuracy
inequality using rational arithmetic.
