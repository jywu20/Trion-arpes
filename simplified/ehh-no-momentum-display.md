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