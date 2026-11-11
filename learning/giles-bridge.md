# From Giles's MLMC review to quantum mean-field dynamics

[Learning route](README.md) · [Scientific overview](../science/README.md) · [Exact model and result](../science/model-and-results.md)

The starting point is an expectation that is useful to know but expensive to compute. A simulation can produce one approximate sample. Improving the simulation reduces its systematic error, while averaging samples reduces statistical uncertainty. Multilevel Monte Carlo separates these two tasks and spends most of its effort on inexpensive corrections.

This project adds an obstacle: the equation being simulated depends on an unknown curve of expectations. The algorithm must first construct that curve accurately enough to simulate the requested final observable. The bridge below follows the cost of doing so, from the classical multilevel allocation to the paid quantum reconstruction.

## 1. The selected source and the relevant route

The educational anchor is Michael B. Giles, **Multilevel Monte Carlo methods**, *Acta Numerica* **24**, 259–328 (2015), [DOI: 10.1017/S096249291500001X](https://doi.org/10.1017/S096249291500001X). Use the [author's corrected 70-page PDF](https://people.maths.ox.ac.uk/gilesm/files/acta15.pdf), identified by the correction notice on the [author's article page](https://people.maths.ox.ac.uk/gilesm/acta/). Its first-page header says 2018; the journal publication is the 2015 article. The section numbers below refer to that corrected file.

| Source section | Idea to take into this project | Project destination |
|---|---|---|
| [§§1.2–1.3, pp. 3–4](https://people.maths.ox.ac.uk/gilesm/files/acta15.pdf#page=3) | Couple a costly approximation to a cheap one; estimate a telescoping sum. | Fine/coarse path corrections |
| [§2.1, Theorem 1, pp. 6–8](https://people.maths.ox.ac.uk/gilesm/files/acta15.pdf#page=6) | Separate bias, correction variance and sample cost. | The exponents and coefficient cost below |
| [§3.1, pp. 21–22](https://people.maths.ox.ac.uk/gilesm/files/acta15.pdf#page=21) | Allocate samples and decide how far to refine. | Certified allocation and confidence budgets |
| [§5.1, starting p. 29](https://people.maths.ox.ac.uk/gilesm/files/acta15.pdf#page=29) | Construct coupled Brownian increments for Euler paths. | Additive-noise path simulation |
| [§5.2, starting p. 33](https://people.maths.ox.ac.uk/gilesm/files/acta15.pdf#page=33) | Relate strong convergence to correction variance; recognize improved classical estimators. | Strong-order-one coupling and classical challenges |

Giles supplies the classical foundation. Quantum mean estimation, reversible computation and the nonlinear feedback construction require additional arguments, supplied in the linked scientific pages. Their primary citations are supporting research references, not additional educational anchors.

The [milestone map](milestones.md) connects these foundations to the full construction, resource accounting, application and comparison.

## 2. One expectation, two errors

Let $`P`$ be a random output, such as a scattering observable of a simulated chain. The target is the number $`\mathbb E[P]`$. A finer numerical simulation produces $`P_\ell`$, where the level $`\ell`$ increases as its time step decreases. Even infinitely many samples of $`P_\ell`$ would recover only $`\mathbb E[P_\ell]`$.

For an estimator $`\widehat p`$,

```math
\mathbb E[(\widehat p-\mathbb E[P])^2]
=\mathrm{Var}(\widehat p)
+\bigl(\mathbb E[\widehat p]-\mathbb E[P]\bigr)^2.
```

The first term is statistical uncertainty. The second is squared bias. More samples reduce the first term; a more accurate simulation is needed to reduce the second. Paying for high accuracy on every sample can be unnecessarily expensive.

Now use the identity

```math
\mathbb E[P_L]
=\mathbb E[P_0]
+\sum_{\ell=1}^{L}\mathbb E[P_\ell-P_{\ell-1}].
```

Define $`D_0=P_0`$ and $`D_\ell=P_\ell-P_{\ell-1}`$ for $`\ell\geq1`$. Each difference is calculated using two approximations driven by the **same** randomness. Independent fine and coarse paths could fluctuate substantially even when their expectations are close. Sharing their Brownian increments makes the actual difference small.

If level $`\ell`$ uses $`N_\ell`$ independent samples of $`D_\ell`$, its contribution to variance is $`V_\ell/N_\ell`$, where $`V_\ell=\mathrm{Var}(D_\ell)`$. If one coupled sample costs $`C_\ell`$, its contribution to work is $`N_\ell C_\ell`$. Samples on different levels can be independent; the important coupling is inside each fine/coarse pair.

For a target sampling variance $`v`$, minimizing work with a real-valued sample allocation gives

```math
N_\ell
=\frac1v\left(\sum_i\sqrt{V_iC_i}\right)
\sqrt{\frac{V_\ell}{C_\ell}},
\qquad
W=\frac1v\left(\sum_i\sqrt{V_iC_i}\right)^2.
```

This follows by differentiating $`\sum_\ell N_\ell C_\ell+\lambda(\sum_\ell V_\ell/N_\ell-v)`$ with respect to each $`N_\ell`$. Enforcing integer counts and at least one sample on each retained level increases this real-valued allocation's work by at most $`\sum_\ell C_\ell`$. That term remains relevant later.

## 3. Translate the notation before comparing rates

There are two different refinements in this project. The path mesh controls how accurately a diffusion is simulated. The temporal support controls how much of the unknown mean curve one coefficient describes. They need not have the same size.

| Symbol | Meaning in this bridge |
|---|---|
| $`\alpha_{\rm ML}`$ | Giles's weak-bias exponent: $`\lvert\mathbb E[P_\ell-P]\rvert=O(2^{-\alpha_{\rm ML}\ell})`$. |
| $`\beta`$ | Correction-variance exponent: $`V_\ell=O(2^{-\beta\ell})`$. This is the variance exponent, not the standard-deviation exponent. |
| $`\gamma`$ | Cost exponent of one sample: $`C_\ell=O(2^{\gamma\ell})`$. It must reflect the actual implementation. |
| $`h_\ell=T2^{-\ell}`$ | Simulation time step; the index $`\ell`$ refines a stochastic path. |
| $`H_j=T2^{-j}`$ | Width of one temporal support; the index $`j`$ refines the representation of the mean curve. |
| $`\delta`$ | Allowed additive error in one scalar temporal coefficient. |
| $`\eta`$ | Allowed error in one complete approximate mean-curve update. |
| $`\epsilon`$ | Allowed error in the final requested classical expectation. |
| $`K`$ | Number of feedback features $`g_1,\ldots,g_K`$; $`d`$ is the state dimension. |
| $`p_{\rm fail}`$ | Failure probability used here to avoid confusing it with $`\alpha_{\rm ML}`$. The main theorem denotes its failure budget by $`\alpha`$. |
| $`L_{\rm path}`$ and $`J`$ | Finest simulation level and finest temporal-reconstruction level, respectively. Neither is the number of beads in the polymer application. |

Giles's Theorem 1 assumes positive exponents, $`\alpha_{\rm ML}\geq\tfrac12\min\lbrace \beta,\gamma\rbrace `$, and the stated bias, variance and expected-cost bounds. Its classical RMS-error work conclusion has the three regimes

```math
W_{\rm MLMC}
=
\begin{cases}
O(\epsilon^{-2}),&\beta\gt \gamma,\\
O(\epsilon^{-2}\log^2(1/\epsilon)),&\beta=\gamma,\\
O(\epsilon^{-2-(\gamma-\beta)/\alpha_{\rm ML}}),&\beta\lt \gamma.
\end{cases}
```

The comparisons express whether finer corrections become cheap enough in variance to offset their higher simulation cost. Applying the theorem requires proving the relevant bounds for the actual estimator. An observed slope on a short numerical plot does not supply a certified constant or a uniform theorem.

Giles's §3.1 implementation estimates level variances and uses observed weak differences to choose further refinement. The theorem in this repository instead chooses its cutoffs and budgets from supplied certified bounds. Empirical allocation is useful in an implementation, but importing its stopping decision would require a separate error guarantee.

The RMS guarantee and a success-probability guarantee are different. This project returns an estimate within $`\epsilon`$ with probability at least $`1-p_{\rm fail}`$. For its classical comparator, a constant-success estimator can be repeated independently and combined by a median. Quantum mean estimation has its own confidence amplification. Both overheads are charged.

## 4. Why additive noise changes the path estimate

A stochastic differential equation specifies a deterministic drift and a random forcing. With a known feedback curve $`v`$, this project simulates

```math
dX_t^v=a(X_t^v,v(t))\,dt+\Sigma\,dW_t,
```

where $`W_t`$ is Brownian motion and $`\Sigma`$ is constant. Brownian increments over disjoint intervals are independent, with covariance equal to the interval length. The Euler update is

```math
X_{n+1}^{v,h}
=X_n^{v,h}+a(X_n^{v,h},v(nh))h+\Sigma\Delta W_n.
```

A coarse Brownian increment is the sum of the corresponding fine increments. This is the concrete coupling used in Giles's §5.1.

Generic Euler simulation has strong order one half under the usual conditions. Strong order describes closeness of coupled random states, not just their averages. The Milstein discussion in §5.2 shows how strong order one can improve correction variance. In the scalar Milstein formula, the extra term contains a spatial derivative of the diffusion coefficient. That term vanishes for constant diffusion.

For the present additive-noise assumptions, the repository proves directly that

```math
\sup_{0\leq t\leq T}
\left(\mathbb E\lVert X_t^v-X_t^{v,h}\rVert^2\right)^{1/2}
\leq C_Eh.
```

The sufficient conditions include two bounded state derivatives of the drift and a uniformly Lipschitz supplied curve $`v`$. A linear-confining extension uses fourth moments. The displayed supremum is over fixed-time mean-square errors; it is not a pathwise supremum inside the expectation. See [the additive-noise Euler estimate](../science/dynkin-details.md#additive-noise-euler-error-with-an-explicit-constant).

For a Lipschitz terminal observable $`f`$, the shared-path difference consequently has $`L^2`$ size $`O(h_\ell)`$ and variance $`O(h_\ell^2)`$. A full path costs $`O(h_\ell^{-1})`$ steps. At fixed model parameters this gives the sufficient classical exponents

```math
\alpha_{\rm ML}=1,\qquad\beta=2,\qquad\gamma=1.
```

Thus a known-feedback expectation already has a classical $`O(\epsilon^{-2})`$ MLMC upper bound. Comparing quantum computation only with a cubic-cost single-level Euler calculation would miss this improvement.

## 5. A complete random-time detail example

Consider the especially simple diffusion

```math
dX_t=dW_t,\qquad X_0=0,\qquad g(x)=\cos x.
```

Its mean curve is known exactly:

```math
\phi(t)=\mathbb E[\cos W_t]=e^{-t/2}.
```

The midpoint detail on $`[0,1]`$ is the signed deviation of the midpoint from the line joining the endpoints:

```math
c=\phi(1/2)-\frac{\phi(0)+\phi(1)}2
=e^{-1/4}-\frac{1+e^{-1/2}}2
\approx-0.02446455.
```

A generator describes the instantaneous expected change of an observable. Here it is $`\mathcal A=\tfrac12\partial_x^2`$, so

```math
\mathcal Ag(x)=-\tfrac12\cos x.
```

Itô's formula gives

```math
\cos W_t
=1-\frac12\int_0^t\cos W_s\,ds
-\int_0^t\sin W_s\,dW_s.
```

The final stochastic integral has mean zero. Taking expectations therefore expresses changes in $`\phi`$ as integrals of its generator expectation. Subtracting the two half-interval integrals produces the midpoint detail.

Let $`\tau`$ be independent and uniform on $`[0,1]`$, and define

```math
s(t)=
\begin{cases}
+1,&0\leq t\lt 1/2,\\
-1,&1/2\leq t\leq1.
\end{cases}
\qquad
Z=-\frac14s(\tau)\cos W_\tau.
```

The factor $`1/4`$ combines the midpoint-detail factor $`1/2`$ with the generator factor $`1/2`$. Direct integration verifies the answer:

```math
\begin{aligned}
\mathbb E[Z]
&=-\frac14\left(
\int_0^{1/2}e^{-t/2}\,dt
-\int_{1/2}^1e^{-t/2}\,dt
\right)\\
&=e^{-1/4}-\frac{1+e^{-1/2}}2=c.
\end{aligned}
```

One sample can be produced by drawing $`\tau`$ and an independent standard normal $`G`$, then evaluating $`Z=-s(\tau)\cos(\sqrt\tau G)/4`$. The output is a signed scalar bounded by $`1/4`$.

For an interval of width $`H`$, the corresponding multiplier is $`-H/4`$, so the second moment is at most $`H^2/16`$. This is the support-width factor the general construction needs. It does not require the Brownian path to be differentiable.

This example is analytically easy. Its exact expectation is already written above, and Brownian motion can be sampled at one time without stepping through its past. Euler has no simulation bias for this process. It illustrates the estimator and its factors, not an advantage instance or a nonzero multilevel-error experiment. A nonlinear drift generally requires the full path prefix from time zero to $`\tau`$, even when the detail occupies a very short interval.

## 6. Put the unknown feedback back into the problem

The actual task is

```math
dX_t=a(X_t,m(t))\,dt+\Sigma\,dW_t,
\qquad m(t)=\mathbb E[g(X_t)],
\qquad z=\mathbb E[f(X_T)].
```

The input contains coefficient programs and a known initial state or an explicitly supplied sampler. It does not contain $`m(t)`$. To make progress, fix a candidate curve $`v`$, simulate the ordinary diffusion $`X^v`$, and define its mean curve

```math
\Phi(v)(t)=\mathbb E[g(X_t^v)].
```

The correct curve satisfies $`m=\Phi(m)`$. This is the fixed-point problem solved by Picard iteration: compute an approximate $`\Phi(v)`$, use it as the next candidate, and repeat.

For component $`g_k`$, write

```math
F_k(x,v)
=a(x,v)\cdot\nabla g_k(x)
+\frac12\mathrm{tr}
\left(\Sigma\Sigma^{\mathsf T}\nabla^2g_k(x)\right).
```

For $`I=[u,u+H]`$, let $`s_I`$ be positive on its first half and negative on its second. The same calculation as in the example gives

```math
\begin{aligned}
c_{k,I}(v)
&=\Phi_k(v)(u+H/2)
-\frac{\Phi_k(v)(u)+\Phi_k(v)(u+H)}2,\\
Y_{k,I}
&=\frac H2s_I(\tau)F_k(X_\tau^v,v(\tau)),
\qquad \tau\sim\mathrm{Unif}(I),\\
c_{k,I}(v)&=\mathbb E[Y_{k,I}].
\end{aligned}
```

The curve also needs an endpoint difference. Its estimator has multiplier $`T`$, rather than $`T/2`$:

```math
\Phi_k(v)(T)-\mathbb E[g_k(X_0)]
=\mathbb E\left[T F_k(X_\tau^v,v(\tau))\right],
\qquad \tau\sim\mathrm{Unif}[0,T].
```

Both have multiplier bounded by their support width. The initial means are known through the charged input-evaluation contract.

In the numerical sampler, simulate to $`\tau`$ with a partial final Euler step. Its last feedback value can be frozen at the preceding grid point. Couple fine and coarse paths with the same endpoint, initial sample and Brownian increments. The generator's Lipschitz bound, the strong order-one estimate and the $`H`$ multiplier yield

```math
\left\lvert\mathbb E[Y_{H,h}]-\mathbb E[Y_H]\right\rvert
\leq C Hh,
\qquad
\lVert Y_{H,h}-Y_{H,2h}\rVert_{L^2}\leq C' Hh.
```

The bounded theorem obtains the needed generator regularity from three bounded derivatives of $`g_k`$. Linear confinement instead uses a weighted generator bound and fourth moments. These are [proved sampling estimates](../science/dynkin-details.md), not assumptions about an already prepared distribution.

## 7. One coefficient: the classical and quantum bills

Normalize $`T`$ to one for the scaling calculations in this section. Its effect on model constants remains charged. For a detail of width $`H`$,

```math
\mathrm{bias}=O(Hh_\ell),\qquad
V_\ell=O(H^2h_\ell^2),\qquad
C_\ell=O(h_\ell^{-1}).
```

The base level has second moment $`O(H^2)`$. Every sample includes the complete path prefix. The required finest step is proportional to $`\min\lbrace 1,\delta/H\rbrace `$, with constants chosen for a bias allowance such as $`\delta/4`$. A small support can therefore use a coarse path grid even when $`h_\ell\gt H`$.

### Classical allocation

Set $`S=\sum_{\ell=0}^{L_{\rm path}}\sqrt{V_\ell C_\ell}=O(H)`$. One conservative integer allocation is

```math
N_\ell=
\max\left\lbrace 1,
\left\lceil
16\delta^{-2}S\sqrt{V_\ell/C_\ell}
\right\rceil\right\rbrace .
```

It gives sampling variance at most $`\delta^2/16`$. Chebyshev's inequality then bounds the probability of statistical error exceeding $`\delta/2`$ by $`1/4`$. The simulation bias and finite-bit errors have separate allowances. Independent repetition and a median supply the desired higher confidence. The minimum of one sample also covers a potentially nonzero deterministic correction.

Before logarithmic confidence factors, the work is

```math
W_{\rm cl}(H,\delta)
=O\left(1+H^2\delta^{-2}+H\delta^{-1}\right).
```

The quadratic term is the allocated sampling work. The remaining terms pay for the retained levels, including the finest sample. This is the classical derivation in [Classical comparison](../science/classical-comparison.md#2-a-direct-classical-analogue-under-the-same-analytic-assumptions).

### Quantum estimation

The quantum ingredient is a mean estimator for a reversibly computed random variable $`D`$ with a supplied second-moment bound $`\mathbb E[D^2]\leq\sigma^2`$. Up to logarithms, accuracy $`e`$ costs $`O(1+\sigma/e)`$ uses of the sample circuit and its inverse. These circuits compute over uniform random bits, and the final answer is an ordinary classical number. The preparation circuit itself has a nonzero gate cost.

Here the correction norm is $`\sigma_\ell=O(Hh_\ell)`$. Assigning error $`e_\ell`$ of order $`\delta/(L_{\rm path}+1)`$ gives level work proportional to

```math
h_\ell^{-1}\left(1+\frac{Hh_\ell}{e_\ell}\right).
```

The $`h_\ell`$ factors cancel in the second term. Summing levels and retaining the mandatory finest-level cost gives

```math
W_{\rm q}(H,\delta)
=\widetilde O(1+H/\delta).
```

The tilde hides logarithmic accuracy and confidence factors. This is quantum MLMC applied to the particular support-scaled hierarchy. The mean-estimation and quantum-MLMC ingredients are attributed to Montanaro and An et al. in the [supporting bibliography](../references/related-work.md) and [circuit account](../science/finite-precision-and-circuits.md).

Crucially, $`O(Hh_\ell)`$ is a second-moment scale. It is not a bound on every possible sample. Dividing a correction by that scale and treating it as a probability would be invalid. The variance-aware quantum estimator and its finite-range/precision treatment are part of the construction.

## 8. From coefficient savings to a complete mean curve

A dyadic piecewise-linear curve can be built from its endpoints and midpoint details. Each detail multiplies a triangular function supported on its interval and equal to one at its midpoint. Adding finer triangles corrects new midpoint values without changing previous nodes.

At scale $`H_j`$, there are $`O(H_j^{-1})`$ intervals. At any one time, only one triangle from that scale is nonzero. Consequently coefficient errors accumulate by the **number of scales**, not by the total number of coefficients. With $`J`$ scales and $`K`$ components, a sufficient choice is

```math
\delta=\frac{\eta}{4\sqrt K(J+1)},
\qquad
J=O(\log(1/\eta))
```

at fixed model parameters. Lipschitz interpolation uses $`O(\eta^{-1})`$ nodes. The reconstruction also controls slopes and projects values into the allowed moment range. Thus every successful next iterate remains an admissible Lipschitz curve for the preceding sampling estimates. See [Feedback reconstruction](../science/feedback-reconstruction.md).

For quantum estimation, the total at one scale is

```math
O(H_j^{-1})\,\widetilde O(1+H_j/\delta)
=\widetilde O(H_j^{-1}+\delta^{-1}).
```

Summing all scales gives $`\widetilde O(\eta^{-1})`$ at fixed $`K`$ and model constants. The first term pays for the actual number of coefficients. The second uses the fact that all supports at a scale have a fixed total width.

Applying the classical coefficient bill to the same representation instead gives

```math
O\left(2^J+\delta^{-2}\sum_j H_j+(J+1)\delta^{-1}\right)
=\widetilde O(\eta^{-2}).
```

The geometric sum of support widths is bounded. This produces a matched classical upper bound without adding a particle-number approximation.

For the nonlinear feedback, a Volterra stability estimate controls successive approximate Picard updates. The error after $`R`$ updates has the form

```math
\lVert v_R-m\rVert_\infty
\leq E_0\frac{A_T^R}{R!}
+C_{\rm stab}\eta,
```

where $`E_0`$ bounds the initial curve error and the other constants depend explicitly on the model and horizon. The factorial term decreases as the iteration count grows. This allows a logarithmic sufficient iteration count at fixed parameters without assuming that a single update is a short-time contraction. Large stability constants still increase the work.

Finally, the computed classical curve $`v_R`$ is fixed while estimating $`\mathbb E[f(X_T^{v_R})]`$. Its feedback error and terminal simulation/estimation error receive separate portions of the final $`\epsilon`$ budget. Conditional confidence bounds are valid for every admissible previously measured curve and are combined over all calls. The result is one scalar with the stated success probability, rather than a quantum state containing an inaccessible answer.

## 9. Preparation and memory remain part of the result

The direct quantum implementation has $`\widetilde O(\epsilon^{-1})`$ logical-gate dependence at fixed model parameters, under the supplied bit-efficient evaluation contract. This includes finite random bits, Gaussian generation, drift and generator arithmetic, controls, the inverse sampler, classical curve reconstruction and circuit generation. Feedback values at deterministic grid times are compiled into the schedule. They are not granted a free coherent table lookup.

There is a distinct choice of coherent memory:

| Implementation | Logical gates at fixed model parameters | Coherent workspace |
|---|---|---|
| Direct retained random tape and history | $`\widetilde O(\epsilon^{-1})`$ | Can be $`\widetilde O(\epsilon^{-1})`$ |
| Space-bounded PRG with ternary reversible checkpoints | $`\widetilde O(\epsilon^{-\log_3 5})`$ | Polynomial in logarithmic accuracy and confidence parameters at fixed model |

The second row replaces the long random tape by an explicitly analyzed short seed and recomputes intermediate states. It preserves both the first and second moments needed by MLMC. A ternary checkpoint bridge uses five recursive subcomputations per three segments, producing $`\log_3 5\approx1.465`$. Its smaller coherent workspace has a higher gate exponent. The two rows must not be combined into an inverse-linear, small-workspace claim.

Both implementations separately store the classical feedback curve, with $`\widetilde O(\epsilon^{-1})`$ bits at fixed model parameters. Dimension, noise norms, stiffness, derivative bounds, parameter precision and the actual coefficient program still affect costs. None of these gate estimates is a physical hardware runtime forecast. The complete accounting is in [Finite precision and circuits](../science/finite-precision-and-circuits.md) and [Memory and recomputation](../science/memory-and-recomputation.md).

## 10. What the bridge establishes, and where the comparison can fail

The inherited classical lesson is to pay separately for approximation bias, coupled statistical uncertainty and sample construction. The added feedback mechanism retains the temporal width in both bias and second moments, chooses a different simulation cutoff for each support, and keeps the measured mean curves regular enough for the next update. These steps explain how the unknown-feedback reconstruction fits the quantum precision bound. Global curve parameterization and quantum mean estimation alone are not the new claim; their earlier uses are identified in [Related work](../references/related-work.md).

The classical $`\widetilde O(\epsilon^{-2})`$ construction is an upper bound. It is not a lower bound on arbitrary classical algorithms given explicit coefficient programs. In particular:

- Low-dimensional quadrature or a PDE solver can be cheaper. The cited smooth one-dimensional mean-field quadrature regime already has essentially inverse-linear cost.
- QMC, multilevel QMC and Brownian-coordinate choices can improve classical sampling. Giles's §5.2 itself discusses such gains on coarse levels; a coarse-level IID variance argument does not exclude them.
- Exact conditional expectations, learned control variates, Gaussian closures or low-rank representations can remove much of the stochastic uncertainty if their construction and errors are inexpensive enough.
- For one fixed analytic model, a sufficiently effective deterministic spectral method may eventually have a better accuracy exponent, despite a large dimension-dependent prefactor.

The [finite-bead polymer task](../science/polymer-dynamics.md) supplies a concrete input, explicit Gaussian initial law and useful density/scattering outputs. Its finite Fourier interaction does not give an exact finite moment closure. That observation motivates the computational question but proves no classical hardness. The illustrative classical calculations and the full [classical comparison](../science/classical-comparison.md) delimit the claim.

The supported conclusion is a constructive quantum accuracy upper bound, with a separate coherent-memory tradeoff, and a conditional route to advantage where stronger classical structure does not resolve the same output more cheaply. The exact assumptions and quantified constants remain in [Model and results](../science/model-and-results.md).
