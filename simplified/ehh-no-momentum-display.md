[`ehh-no-momentum-display.jl`](ehh-no-momentum-display.jl) is based on [`ehh-heatmap-prototype.jl`](ehh-k1-contribution.md#the-contributions-of-s-to-the-final-arpes-heatmap).
The output is 

![](ehh-no-momentum-display.png)

Note that the vertical gray dotted lines correspond to the momentum of the intensity plots (b, d, e).

We can also plot the panels in a more compact form.
This is done by [`ehh-no-momentum-display-compact.jl`](ehh-no-momentum-display-compact.jl).
The result is shown below:

![](ehh-no-momentum-display-compact.png)

Here we also include [`ehh-shifted-momentum-display.md`](ehh-shifted-momentum-display.md).

We can actually incorporate the [finite-Q calculation, i.e. `ehh-shifted-momentum-display.md`](ehh-shifted-momentum-display.md) into this plot.
This is done in [`ehh-display-very-compact.jl`](ehh-display-very-compact.jl).
The output is shown below:

![](ehh-display-very-compact.png)

To make the frequency positions of the signatures more clear,
we use [`ehh-display-very-compact-2.jl`](ehh-display-very-compact-2.jl) to plot 

![](ehh-display-very-compact-2.png)

Finally, we may want to incorporate the interactive corrections within the figure,
and this is done using [`ehh-display-very-compact-2-residue-energy.jl`](ehh-display-very-compact-2-residue-energy.jl):

![](ehh-display-very-compact-2-residue-energy.png)

Note: all figures except the last are plotted using the 3D trion wave function
which is theoretically wrong although it doesn't make things qualitatively wrong.

Run [`ehh-display-very-compact-2-residue-energy-ionization.jl`](ehh-display-very-compact-2-residue-energy-ionization.jl) to get 

![](ehh-display-very-compact-2-residue-energy-ionization.png)

If we choose to not use a hard coded Ehh but a $Z(E)$ correction factor to the spectral weight factor to represent how Coulomb repulsion between the holes scatters the two holes to the free-hole-pair final state,
we can run [`ehh-display-very-compact-2-residue-energy-distribution-ionization.jl`](ehh-display-very-compact-2-residue-energy-distribution-ionization.jl) and get 

![](ehh-display-very-compact-2-residue-energy-distribution-ionization.png)

The downward shift of the trion signature due to the Z(E) factor, 0.197 eV, is also reported when running the script.