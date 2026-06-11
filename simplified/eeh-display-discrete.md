ARPES imaging of eeh trions, with excitons forming in the residue state
===========

Running [`eeh-heatmap-prototype-discrete-final.jl`](eeh-heatmap-prototype-discrete-final.jl),
we get 

![](eeh-heatmap-prototype-discrete-final.png)

To make the figure look more "elaborated", run [`eeh-heatmap-prototype-discrete-final-zoom-in-display.jl`](eeh-heatmap-prototype-discrete-final-zoom-in-display.jl),
and we get 

![](eeh-heatmap-prototype-discrete-final-zoom-in-display.png)

Or we can try to break the y axis:
run [`eeh-heatmap-prototype-discrete-final-zoom-in-display-fold.jl`](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold.jl), and we get 

![](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold.png)

Run [`eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-fully-labeled.jl`](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-fully-labeled.jl) to have a fully labeled version with info like "this is the A series unlike-spin 1s exciton" included into the figure.

We note that the intensity of the 2s peak does not vanish.
This is not due to visual illusions.
If we plot the exciton wave functions with respect to $\abs{\vb{k}}$,
this is what we get for the 1s state:

![](exciton-0-patch-2.png)

and this is what we get for the 2s state:

![](exciton-0-patch-10.png)

And it can be observed that the absolute intensity of the negative center of the 2s wave function happens to be higher than that of the positive center of the 1s wave function.
After computing the overlap between the exciton wave functions and the trion wave function,
here is what we get for the 1s state (the third panel is the trion wave function):

![](exciton-0-patch-overlap-2.png)

and here is what we get for the 2s state (again, the third panel is the trion wave function):

![](exciton-0-patch-overlap-10.png)

Note that the structure of the trion wave function means that for the 2s state, it is the negative central part that contributes most to the overlap matrix element;
the contribution from the positive part of the 2s exciton wave function being weak 
means that there is almost no cancellation from the influences of the positive and negative parts of the 2s exciton wave function,
which, together with the fact that the central part of the 2s exciton wave function has a higher intensity,
means the overlap between the exciton and the trion wave functions should not vanish,
although it's still small compared with the 1s peak.

Some of the images below are produced using an incorrect exciton phase correction function, and you may find the 2s peak to be higher than the 1s peak.
This is not physical, and the figures used in the paper are based on the correct implementation.

A detailed look of the 1s and 2s signals: [`eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-K.jl`](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-K.jl) and [`eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-Kp.jl`](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-Kp.jl)

![](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-K.png)

![](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-Kp.png)

These plots only show the first couple exciton bands;
when we include the first 200 bands into the calculation,
by running [`eeh-heatmap-prototype-discrete-final-lots-of-bands.jl`](eeh-heatmap-prototype-discrete-final-lots-of-bands.jl),
we get 

![](eeh-heatmap-prototype-discrete-final-lots-of-bands.png)

Increasing the number of bands from 200 to 300 only makes the high-energy continuum wider.

Run [`eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-K-lots-of-bands.jl`](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-K-lots-of-bands.jl) and [`eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-Kp-lots-of-bands.jl`](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-Kp-lots-of-bands.jl) to get the following plots:

![](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-K-lots-of-bands.png)

![](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-linecut-Kp-lots-of-bands.png)

Note that there are subtleties in how the peak heights are calculated, because the exciton wave function involves two bands in the output of QE and the phase correction is non-trivial.
Run [`exciton-patch-overlap-0.jl`](exciton-patch-overlap-0.jl) to see the overlap matrix elements visualized.
The above two images are generated using the correct exciton wave function phase correction scheme.

We note that the gray frequency-dependent figures in all figures above are *not* momentum-integrated.
For momentum-integrated EDC, run [`eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-K-lots-of-bands-integrated.jl`](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-K-lots-of-bands-integrated.jl) to get 

![eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-K-lots-of-bands-integrated.png](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-K-lots-of-bands-integrated.png)

and run [`eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-Kp-lots-of-bands-integrated.jl`](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-Kp-lots-of-bands-integrated.jl) to get 

![](eeh-heatmap-prototype-discrete-final-zoom-in-display-fold-Kp-lots-of-bands-integrated.png)

Note that the linear exciton bands give the momentum integrated EDC an asymmetric line shape.