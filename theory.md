# Discussion 

- Trion energy: Berkelboob's paper
- Ignore the wave function for now 
- Two electrons, or two holes?
- What happens when we have series tetron?



# Exciton 

$$
\begin{aligned}
P(t) &= \sum_{S, \vb{Q}} \rho_{S \vb{Q}} \abs{A_{\vb{k}'}^{S \vb{Q}}}^2 \abs{M^{fc}_{\vb{k} \vb{k}'}}^2 
\int_{t_0}^t \dd{t_1} \int_{t_0}^t \dd{t_2}
\ee^{-\ii (E_{S \vb{Q}} + \epsilon_{v \vb{k}' - \vb{Q}} - \omega ) (t_1 - t_2)} 
s(t_1) s(t_2)
 \\ 
&= \abs{M^{fc}_{\vb{k} \vb{k}'}}^2 
\sum_{S, \vb{Q}} \rho_{S \vb{Q}} \abs{A_{\vb{k}'}^{S \vb{Q}}}^2 
\cdot \mathrm{broaden}(E_{S \vb{Q}} + \epsilon_{v \vb{k}' - \vb{Q}} - \omega )
.
\end{aligned}
$$

# Trion: one electron, two holes 

## Brief analysis without considering the Hamiltonian

Again consider the case of 2H TMD. We have 
$$
\epsv{\vb{k}} = - \frac{\vb{k}^2}{2 \mh}, \quad 
\epsc{\vb{k}} = \Eg + \frac{(\vb{k} - \vb{w})^2}{2 \me} .
$$
The bare trion energy, which is the sum of the energies of electrons and holes, is 
$$
E_{\text{trion, bare}} = \epsc{\ke} - \epsv{\khi{1}} - \epsv{\khi{2}}.

$$
Here the sign notations are:
- $\mh > 0$
- $\kh$ are the momenta read from the band diagram; hence $\vb{P} = \ke - \khi{1} - \khi{2}$;

The bare trion energy should be minimized to see around which $\ke$ we have highest ARPES intensity.
Thus we have the following problem:
$$
\min E_{\text{trion, bare}} =  
\Eg + \frac{(\ke - \vb{w})^2}{2 \me} + \frac{\khi{1}^2}{2 \mh} + \frac{\khi{2}^2}{2 \mh}
\quad \text{s.t.} \quad 
\vb{P} = \ke - \khi{1} - \khi{2}.
$$
Here only $\khi{1, 2}$ are variables; 
$\vb{P}$ is controlled by the initial state, 
and $\ke$ is controlled by where we are in the AREPS heatmap.

Lagrangian multiplier method tells us that we should let $\khi{1} = \khi{2}$ 
to maximize the bare trion energy.
The exact value of the minimal bare exciton energy is not important 
because it doesn't appear in the ARPES spectrum anyway.
Thus we find that the trion ARPES signature should roughly follow the following dispersion relation:
$$
\omega = E_{S \vb{P}} - 2 \cdot \epsv{(\ke - \vb{P}) / 2}.
$$
Two immediate consequences:
- A curvature different from that of an exciton ARPES signature 
- Different center positions 

Next question is the $\vb{P}$ of trion modes with the lowest energies.
It seems the best idea is to let $\ke = \vb{w}$, and $\khi{1, 2} = 0$,
and in this way $\vb{P} = \vb{w}$ when the trion energy is minimized.
So for a "low-temperature" signature, we have 
$$
\omega = E_{S, \vb{P} = \vb{w}} - 2 \cdot \epsv{(\ke - \vb{w}) / 2}.
$$
That's to say, a signature at *half of* $\vb{w}$ appears; 
that's in sharp contrast with the case of excitons, 
where the ARPES signature appears at $\ke = \vb{w}$.

## Trion Hamiltonian in the two-band model

When we build the effective trion Hamiltonian it seems to be a good idea 
to turn the the sign convention that treats holes as "real" quasiparticles; 
thus, in this section, 
a hole with momentum $\kh$ means a hole placed at $- \kh$ on the band plot.
The momentum conservation relation now is 
$$
\vb{P} = \ke + \khi{1} + \khi{2}. \label{P-conservation}
$$

The Hamiltonian is 
$$
H = \frac{(\ke - \vb{w})^2}{2 \me} + \frac{\khi{1}^2}{2 \mh} + \frac{\khi{2}^2}{2 \mh} 
+ V(\rhi{1} - \rhi{2}) - V(\rhi{1} - \re) - V(\rhi{2} - \re).
\label{negative-trion-H}
$$
There should be no controversy regarding 
$$
\vb{r}_1 = \rhi{1} - \re , \quad 
\vb{r}_2 = \rhi{2} - \re.
$$
We need to identify the "external" degrees of freedom of the trion
and fix them when calculating the ARPES signature of one trion mode
and thus we want a good definition of $\vb{R}$ (and hence $\vb{P}$).
The definition however is definite dependent to the Hamiltonian.
Whether we still have 
$$
\vb{R} = \frac{
    \me \re + \mh (\rhi{1} + \rhi{2})
}{
    M
} \label{R-def-direct}
$$
is a question. Another way to ask the question is why we have $\eqref{R-def-direct}$ when $\vb{w} = 0$.
We may use the criterion 
$$
\dv{\vb{R}}{t} = \frac{1}{\ii \hbar} \comm{\vb{R}}{H} = \frac{\vb{P}}{\underbrace{\me + 2\mh}_M}, \label{R-EOM-P-relation}
$$
where $\vb{P}$ is given by $\eqref{P-conservation}$.
$\eqref{R-EOM-P-relation}$ is correct when $\vb{w} = 0$.
Here we assume that $\comm{x}{p} = \ii \hbar$ is always correct
-- as this is what is used for the semi-classical EOM 
$$
\dv{\vb{r}}{t} = \pdv{E_{\vb{k}}}{\vb{k}} , \quad 
\dv{\vb{k}}{t} = - \pdv{V_{\text{external}}}{\vb{r}} = - \pdv{(- q \vb{r} \cdot \vb{E})}{\vb{r}} = q \vb{E}
$$
of wave packets.
Note that $\comm{x}{p} = \ii \hbar$ is equivalent to $\comm{x}{p - \const.} = \ii \hbar$

When $\vb{w} \neq 0$, we have 
$$
\frac{1}{\ii \hbar} \comm{
    \frac{
        \me \re + \mh (\rhi{1} + \rhi{2})
    }{
        M
    }
}{
    H
} = \frac{\ke - \vb{w} + \khi{1} + \khi{2}}{M} = \frac{\vb{P} - \vb{w}}{M}.
$$
It seems if we choose $\eqref{R-def-direct}$, then we have 
$$
\dv{\vb{R}}{t} = \frac{\vb{P} - \vb{w}}{M}. \label{R-EOM-P-relation-with-w}
$$
This seems reasonable for me: recall that when $\vb{P} = \vb{w}$, the trion energy is minimized.
Thus we should expect to see something like this:
$$
H = \frac{(\vb{P} - \vb{w})^2}{2M} + \cdots,
$$
and this is consistent with $\eqref{R-EOM-P-relation-with-w}$.

Since $\eqref{R-def-direct}$ is still true even when $\vb{w} \neq 0$,
Eqs. (6-9) in Variationally optimized orbital approach to trions in two-dimensional materials
are still correct.
Thus we have 
$$
\khi{1} = \vb{k}_1 + \frac{\mh}{M} \vb{P}, \quad 
\khi{2} = \vb{k}_2 + \frac{\mh}{M} \vb{P}, \quad 
\ke     = \frac{\me}{M} \vb{P} - \vb{k}_1 - \vb{k}_2,
$$
and the total Hamiltonian $\eqref{negative-trion-H}$ can be verified to be (TODO: double check) 
$$
H = \frac{(\vb{P} - \vb{w})^2}{2 M} + \Eg 
+ \frac{\vb{k}_1^2}{2\mh} + \frac{\vb{k}_2^2}{2\mh}
- \frac{\vb{w}^2}{2M} + \frac{(\vb{k}_1 + \vb{k}_2 + \vb{w})^2}{2\me}
- V(\vb{r}_1) - V(\vb{r}_2) + V(\vb{r}_1 - \vb{r}_2).
$$
It can be seen that when $\vb{w}$ is chosen to be zero, 
this Hamiltonian is two direct-gap exciton Hamiltonians plus one interaction term;
but when $\vb{w} \neq 0$, 
first, there is a constant term containing $\vb{w}$, 
and second, the exciton Hamiltonian seems to be 
$$
H = \frac{\vb{k}_1^2}{2\mh} + \frac{(\vb{k}_1 + \vb{w})^2}{2\me} - V(\vb{r}_1).
$$
It can be seen that  
