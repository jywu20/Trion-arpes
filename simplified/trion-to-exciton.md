ARPES spectrum of trion, with discrete excitonic final states
===========

# With a single exciton state 

Run [`eeh-heatmap-prototype-discrete-final-single-exciton.jl`](eeh-heatmap-prototype-discrete-final-single-exciton.jl),
and we get the follows:

![](eeh-heatmap-prototype-discrete-final-single-exciton.png)

Note the function invocation takes the following form:
```julia
Akω_total = trion_ARPES_eeh(trion, P, Ak1k2, IntraValley2DExciton, [exciton_direct], [Avck_A1s_bright], [rk], k1_list, ω_list, broaden)
```
Here we have a caveat related to the language features of Julia.
The array `[exciton_direct]` can be enriched by other types of excitons,
and they make up an array of things with a shared abstract type only.
This in turn means we have to pass the shared type (here it's a concrete type `IntraValley2DExciton`, but it can be an abstract type) to the function.

To make sure the exciton wave function is read properly,
we also draw what we read from the exciton eigenvectors file:

![](exciton-patch.png)

# The linear-parabolic splitting

Run [`exciton-band-gamma-Qiu-data.jl`](exciton-band-gamma-Qiu-data.jl),
we get 

![](exciton-bands.png)

This means the solver works as expected.