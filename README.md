# TKMPMeasure.jl

This package tries to hide the complexity inherent in measure relaxations of Lasserre's hierarchies.

Define a relaxation of a measure on X as
```
@variable(m, μ, TKMPMeasure, on=X, maxdegree=4)
```
the appropriate moment and localizing matrices are created automatically. `μ` now behaves as a de facto measure and can be used in SOS constraints:
```
@constraint(m, ∫(p, μ) >= 0, domain=Y, maxdegree=4)
```
while, in reality, `μ` is only a TKMP-constrained vector masquerading as a measure.