# Units 

1. Unit of length: Å. Thus the length of wave vector is 1/Å.
2. Unit of energy: eV. 
3. Unit of time: by setting $\hbar = 1$. 
   Thus 1fs × $n$ eV / $\hbar$ = $n$ × the time quantity corresponding to 1fs, 
   and we find 1fs is equivalent to 
   $$
   10^{-15} \cdot \frac{e}{\hbar} = 1.5193 \text{ time unit here}.
   $$
4. Unit of mass: the electron mass. 
   Thus an effective mass dispersion relation, in SI, is 
   $$
   E = \frac{\hbar^2 \vb{k}^2}{2\me^*} 
   = \underbrace{\frac{\hbar^2 (1 Å^{-1})^2}{2 \me \cdot 1 \mathrm{eV}}}_{3.80998 \mathrm{eV}} \cdot \frac{(\vb{k} / 1 Å^{-1})^2}{2 \me^* / \me},
   $$
   and therefore when the unit of $\vb{k}$ is set to $Å^{-1}$ and the unit of the effective mass is set to $\me$,
   the dispersion relation, in eV, becomes 
   $$
   E = 3.80998 \frac{\vb{k}^2}{2\me^*}.
   $$ 

# Discussion 

- Trion energy: Berkelboob's paper
- Ignore the wave function for now 
- Two electrons, or two holes?
- What happens when we have series tetron?

# Outline of the project 

1. The peak positions are controlled by $\vb{w}$ and $\vb{P}$
2. Decisive evidence for existence of trion: the position of the intensity peak.
   The difference can be seen when $\vb{w} \neq 0$ or when $\vb{P} \neq 0$.
3. The most robust criterion: the position of the center of the overall signature, 
   (which is just the center of the trion mode with the lowest energy )
   when the gap between the predominant valley and the predominant peak is indirect.
4. Even when the 

# Formalism of time-resolved ARPES 

1. Most general theory: transition rate $P \propto G^<$
2. Assuming that the pump stops before the probe starts 
3. Ignoring "self-driven Floquet" and things like that:
  the initial state of probing is assumed to be a real stationary state, 
  or at least $\rho = \sum_n \rho_n \dyad{{\Psi_n}}{{\Psi_n}}$

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

# Trion: one electron, two holes on the same band 

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
Of course now the relation between $\vb{r}$ and $\vb{P}$ is not linear, 
but this is just the consequence of the EOM given by an arbitrary $\epsilon_{\vb{k}}$. 

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
\label{exciton-ham-in-trion-finite-w} 
$$
This should be double checked; it doesn't look very elegant because of the additional $\vb{w}^2 / 2\me$ term.
The sign of $\vb{k}_1$ shouldn't be a particularly big problem 
since we can always add a minus sign before $\vb{r}_{1, 2}$.

Anyway if this Hamiltonian is correct, 
what does it mean to the intensity peak of $\phi_{\vb{k}_{1, 2}}$?
We can again try to minimize $\eqref{exciton-ham-in-trion-finite-w}$:
the result is 
$$
\underbrace{
    \left( \frac{1}{\me} + \frac{1}{\mh}\right) 
}_{1 / \mu} \vb{k}_1 + \frac{\vb{w}}{\me} = 0, \quad 
\vb{k}_1 = - \frac{\mh}{\me + \mh} \vb{w},
$$
and therefore at the brightest spot in the ARPES signature of a single $\vb{P}$ mode, we have 
$$
\ke = \frac{\me}{M} \vb{P} + \frac{2 \mh}{\me + \mh} \vb{w}.
$$
The center of the overall signature is at 
$$
\ke = \frac{\me}{M} \vb{w} + \frac{2 \mh}{\me + \mh} \vb{w} .
$$
Anyway it's not the same as $\vb{w}$.


This prediction disagrees with the [[theory#Brief analysis without considering the Hamiltonian|naive analysis]]:
the reason is we have an additional $\vb{k}_1 \cdot \vb{k}_2 / \me$ term, 
which is of the same order of magnitude of $H(\vb{k}_1)$ and $H(\vb{k}_2)$.
Therefore the "perturbation" is not small at all; 
this might be seen as evidence for the correctness of [[theory#Brief analysis without considering the Hamiltonian|naive analysis]];
what analysis is closer to the truth is not clear.

The kinetic part of the Hamiltonian looks like 
$$
H_0 = \frac{\vb{k}_1^2}{2\mu} + \frac{\vb{k}_1 \cdot \vb{w}}{\me} + \frac{\vb{k}_2^2}{2\mu} + \frac{\vb{k}_2 \cdot \vb{w}}{\me} 
+ \frac{\vb{k}_1 \cdot \vb{k}_2}{\me} + \const
$$
we want to rewrite this into something like 
$$
\frac{(\vb{k}_1 - \vb{u})^2}{\const} + \frac{(\vb{k}_2 - \vb{u})^2}{\const} + \frac{(\vb{k}_1 - \vb{u}) \cdot (\vb{k}_2 - \vb{u})}{\const} 
$$
It seems we can have 
$$
H_0 = \frac{(\vb{k}_1 - \vb{u})^2}{2\mu} + \frac{(\vb{k}_2 - \vb{u})^2}{2\mu} + \frac{(\vb{k}_1 - \vb{u}) \cdot (\vb{k}_2 - \vb{u})}{\me} + \const, \quad 
\vb{u} = - \frac{\vb{w}}{1 + \frac{\me}{\mu}}.
$$
The energy minimum therefore appears at 
$$
\vb{k}_1 = \vb{k}_2 = \vb{u}, 
$$
which leads to 
$$
\ke = \frac{\me}{M} \vb{P} + \frac{2 \mh}{M} \vb{w}.
$$
It seems this is also the result we obtain from 
$$
\min E_{\text{trion, bare}} =  
\Eg + \frac{(\ke - \vb{w})^2}{2 \me} + \frac{\khi{1}^2}{2 \mh} + \frac{\khi{2}^2}{2 \mh}
\quad \text{s.t.} \quad 
\vb{P} = \ke - \khi{1} - \khi{2}.
$$
The problem can be alternatively written as 
$$
\min \frac{(\ke - \vb{w})^2}{2 \me} + \frac{(\ke - \vb{P})^2}{4 \mh}.
$$
The maximum is achieved when 
$$
\ke = \frac{\frac{\vb{w}}{\me} + \frac{\vb{P}}{2\mh}}{\frac{1}{\me} + \frac{1}{2\mh}} = \frac{\me}{M} \vb{P} + \frac{2 \mh}{M} \vb{w},
$$
the same as the one we obtained in the first approach.

The problem however is whether $\phi(\vb{k}_1 - \vb{u}, \vb{k}_2 - \vb{u})$ 
really has its peak at $\vb{k}_{1, 2} = \vb{u}$.
This however is a clear, well-defined question.

# Trion, one electron and two holes; two holes are on two peaks 

$$
H_{\text{free}} = \frac{\ke^2}{2 \me} + \Eg + \frac{\khi{1}^2}{2\mh} + \frac{(\khi{2} - \vb{w})^2}{2 \mh} 
$$
The constraint is (here to keep $\epsv{}$ consistent with the form in band diagrams,
no additional minus sign is placed before $\vb{k}$ of hole)
$$
\vb{P} = \ke - \khi{1} - \khi{2}.
$$
Again, minimizing the energy (this energy minimization condition should be verified more carefully: maybe by really solving a trion Hamiltonian?), we find the peak intensity in the signature of one exciton mode appears when  
$$
\frac{\ke}{\me} + \frac{\khi{1}}{\mh} = 0, \quad
\frac{\ke}{\me} + \frac{\khi{2} - \vb{w}}{\mh} = 0,
$$
and therefore 
$$
\ke = \frac{\me (\vb{P} + \vb{w})}{\me + 2 \mh}, \quad 
\khi{1} = - \frac{\mh (\vb{P} + \vb{w})}{\me + 2 \mh}, \quad 
\khi{2} = - \frac{\mh \vb{P} - \me \vb{w} - \mh \vb{w}}{\me + 2 \mh}.
$$
Since the energy is globally minimized at $\vb{P} = - \vb{w}$,  
the overall position of the light spot should be 
$$
\ke  = 0. 
$$

Does this mean that in order to observe obvious distinction between excitons and trions, 
there should be at least some sort of indirect band gap in the system?



# Trion, two electrons 

Assuming that the two electrons are in different valleys, 
the interference between the two is non-existent, 
and we find 
$$
P(t) \propto \int_0^{t} \dd{t_1} \int_0^{t} \dd{t_2} s(t_1) s(t_2) \left(
    \ee^{- \ii (t_1 - t_2) (E_n - \epsc{2} + \epsv{} - \omega) \abs{A_{\vb{k}' c_2 v}}} + 
    \ee^{- \ii (t_1 - t_2) (E_n - \epsc{1} + \epsv{} - \omega) \abs{A_{c_1 \vb{k}'  v}}}
\right).
$$
And 

# Challenges 

- [ ] commutation relation for holes - still $[x, p] = \ii \hbar$?
- [ ] The 1BZ: what to do in the hexagonal case?
- [ ] Real space or momentum space? 

Current procedure: 
1. Prove $\vb{P} = \ke + \khi{1} + \khi{2}$ is a constant (always the case from the form of Coulomb interaction)
2. Define $\vb{R}$ that is the corresponding coordinate variable of $\vb{P}$ - 
   here is a subtlety: if $\comm{x}{p}.= \ii \hbar$, then so is $\comm{x + \const}{p}$.
   But I believe 
   $$
   \vb{R} = \frac{\mh \rhi{1} + \mh \rhi{2} + \me \re}{2\mh + \me}
   $$
   is always a good choice, regardless of the positions of the maxima and minima of the valleys and peaks: 
   we always have 
   $$
   \dv{\vb{R}}{t} = \frac{1}{\ii \hbar} \frac{\vb{P} - \sum \vb{w}}{\sum m},
   $$
   which is desirable.
   Note that since $\partial_{\vb{R}}$ is not well-defined when the other variables are not specified, 
   we can explicitly *define* that $\vb{P} = -\ii \hbar \partial_{\vb{R}}$; 
   this equation can be used below to finally settle down the definitions of the other variables.
3. Use $\vb{r}_{1, 2} = \rhi{1, 2} - \re$ as coordinates
4. Write $\ke, \khi{1}, \khi{2}$ in terms of $\partial_{\vb{r}_{1, 2}}$ and $\partial_{\vb{R}}$ - the constant in $\vb{R}$ should disappear here; note that the momentum variables corresponding to $\vb{r}_{1, 2}$ may differ from $\partial_{1, 2}$ by constants.
   Now since the form of $\vb{R}$ doesn't depend on the positions of the valleys and peaks, 
   we *always* have 
   $$
   \ke = \ii \hbar \partial_{\vb{r}_1} + \ii \hbar \partial_{\vb{r}_2} + \frac{\me}{M} \vb{P}, \quad 
   \khi{1, 2} = - \ii \hbar \partial_{\vb{r}_{1, 2}} + \frac{\mh}{M} \vb{P}.
   $$
5. Write the kinetic energy in terms of $\partial_{1, 2}$; decide the definition of $\vb{p}_{1, 2}$, that is, the constant differences between them and $\partial_{\vb{r}_{1, 2}}$
6. The resulting problem about $\vb{r}_{1, 2}$ and $\vb{p}_{1, 2}$ should have nothing different from a problem about two coupled "ordinary" electrons (i.e. electrons whose energy minimum appears at $\vb{k} = 0$)

For convenience, we use $\vb{k}_{1, 2}$ to refer to $- \ii \hbar \partial_{\vb{r}_{1, 2}}$.
We want to find the minimum of 
$$
H_0 = \frac{(\vb{P} - \vb{w})^2}{2 M} + \Eg 
+ \frac{\vb{k}_1^2}{2\mh} + \frac{\vb{k}_2^2}{2\mh}
- \frac{\vb{w}^2}{2M} + \frac{(\vb{k}_1 + \vb{k}_2 + \vb{w})^2}{2\me}
$$
w.r.t. $\vb{k}_{1, 2}$ so that we can redefine $\vb{k}_{1, 2}$.
Then we see indeed that the peaks in the ARPES signature can be found by minimizing terms containing $\vb{k}_1$ in $H_0$ and terms containing $\vb{k}_2$ independently. 

# Note on applying Lagrangian multiplier 

The optimization problem 
$$
\min f(x) + g(y) \quad \text{  s.t. } \abs{h(x, y)} = 0
$$
can be solved by 
$$
\partial_{x, y, \lambda} (f(x) + g(y) - \lambda h(x, y)^2) = 0.
$$
Thus we get 
$$
\partial_x f - 2 \lambda h \partial_h = 0, \quad \partial_y g - 2 \lambda h \partial_y h = 0, \quad h = 0.
$$
The problem is, the expression $\partial_x f - 2 \lambda h \partial_h = 0$
is intended to mean that the gradient of $f$ and $h^2$ are parallel, 
and when the gradient of $h^2$ is zero, 
the condition that the gradient of $f$ and $h^2$ are parallel is trivially true, 
but the equation $\partial_x f - 2 \lambda h \partial_h h = 0$ is no longer true; 
we may want to put $\lambda$ before $\partial_x f$ in this case.