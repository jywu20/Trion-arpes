# The contributions of $\vb{k}_1$s to the final ARPES heatmap

Running [`ehh-heatmap-prototype.jl`](ehh-heatmap-prototype.jl), we get 
![Trion ARPES with $\vb{P} = 1.2 \vb{w}$](ehh-P-1.2.png)

In this figure the trion momentum is $\vb{P} = 1.2 \vb{w}$.
It can be observed that the ARPES signature includes a plurity of dispersion relations.

This figure is also produced in the rightmost panel in the output of [`ehh-k1-contribution-benchmark.jl`](ehh-k1-contribution-benchmark.jl),
which is to [benchmark the output of the new program with the old program](benchmarks.md#comparison-with-old-codes).

# The structure of the contribution of each $\vb{k}_1$

In the ARPES formula we have a dispersion relation and a wave function structural factor from the trion wave function.
Running [`ehh-k1-wfn-and-dispersion-prototype.jl`](ehh-k1-wfn-and-dispersion-prototype.jl),
we get the following plot:

![The dispersion relation and wave function structural factor of a single $\vb{k}_1$](ehh-ring-and-peak-contribution-example.png)

It can be observed that 
- In the plane of $\vb{k}_1$, the dispersion relation is like a ring,
while the wave function structural factor has two peaks.
- The peaks are close to (but not identical to) where $\vb{k}_1 = 0$ or $\vb{k}_2 = 0$, which are represented by crosses. The reason is when $\vb{k}_1 = 0$, $\vb{k}_2$ is suboptimal w.r.t. the wave function of another momentum degree of freedom.
- The center of the dispersion relation circle is the same as the center of the image of the wave function structural factor. This is a consequence of equivalence of $\vb{k}_{1,2}$.

The image above fixes $\omega$, so that we can have a 2D representation of the dispersion relation.
When we move $\omega$, the dispersion relation circle changes its size (but not position),
and correspondingly the final contribution from one $\vb{k}_1$ changes its strength.
Preferrably, the dispersion relation circle should go close to the peaks of the wave function structural factor.
