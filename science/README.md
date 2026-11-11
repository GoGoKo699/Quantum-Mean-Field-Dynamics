# Scientific overview

[Home](../README.md) · [Exact result](model-and-results.md) · [Classical comparison](classical-comparison.md) · [Related work](../references/related-work.md)

## An expectation that determines its own dynamics

A mean-field model describes a representative random system whose force depends on population averages. For a polymer chain, those averages can be a few Fourier density modes. For other models, they may be specified bounded features of the state. The computational input provides the coefficient formulas and initial distribution; the unknown population averages must still be found.

The output in this project is one additive-error estimate of a bounded observable at a finite time. That choice matters. The algorithm need not recover a probability density, print a trajectory ensemble or reconstruct a quantum state. It does need to obtain the feedback used by the dynamics from the supplied classical input.

Quantum estimation can reduce the sample dependence of an ordinary expectation, but substituting it into every sequential time update can repeatedly rebuild long path histories. A useful algorithm must control that accumulated cost. The construction here instead computes a classical approximation to the entire feedback curve and then uses it as the known coefficient of an ordinary diffusion. The previous curve is always classical when a new quantum estimation problem begins.

## The central identity

Let $v(t)$ be a known candidate mean curve and let $X^v$ solve the resulting ordinary diffusion. Its next Picard curve is $\Phi(v)(t)=\mathbb E[g(X_t^v)]$. For one scalar feature $g_k$, define its generator observable

```math
F_k(x,m)=a(x,m)\cdot\nabla g_k(x)
+\frac12\operatorname{tr}(\Sigma\Sigma^\mathsf T\nabla^2g_k(x)).
```

On an interval $I=[u,u+H]$, the midpoint detail of $\Phi_k(v)$ is

```math
c_{k,I}=\Phi_k(v)(u+H/2)
-\frac{\Phi_k(v)(u)+\Phi_k(v)(u+H)}2.
```

Dynkin's formula turns this into

```math
c_{k,I}=\mathbb E\!\left[
\frac H2 s_I(\tau)F_k(X_\tau^v,v(\tau))
\right],\qquad \tau\sim\operatorname{Unif}(I),
```

where $s_I$ is positive on the first half and negative on the second. This identity estimates a curve detail without separately estimating three nearby expectations and subtracting them. Its $H$ multiplier is retained through the full path simulation analysis.

With additive noise, a continuous Euler approximation has fixed-time strong error $O(h)$ under the stated finite regularity. The fine and coarse paths share both the random terminal time and all Brownian increments, including their final partial step. The payoff bias and the $L^2$ norm of their difference are consequently both $O(Hh)$. [Dynkin details](dynkin-details.md) proves the constants, the coupling, and the linear-confining extension.

## Why the cost of the whole curve stays controlled

At coefficient tolerance $\delta$, the support-width factor permits the finest simulation mesh to scale as $h\asymp\delta/H$. A wide support needs a fine path; a short support can use a coarser one. Each call still simulates its prefix from time zero, and the cost of every call is included.

Quantum multilevel estimation uses the second-moment bound on the coupled corrections, rather than incorrectly treating their standard deviation as a maximum possible value. The cost of one coefficient is $\widetilde O(1+H/\delta)$. At a temporal scale of width $H$, there are $O(H^{-1})$ coefficients. Their total cost is therefore $\widetilde O(H^{-1}+\delta^{-1})$. The temporal resolution and coefficient tolerance are both linked to the requested curve error, so summing all scales retains one inverse power of accuracy.

Only one dyadic hat at each scale is nonzero at any fixed time. This keeps the reconstruction error proportional to the number of scales, rather than to the total number of stored coefficients. The same observation bounds the slopes of the reconstructed curve. Projection of its nodes onto the allowed moment range preserves that regularity bound.

Global Picard iteration supplies the remaining nonlinear argument. Differences propagate through a Volterra integral operator; its repeated powers have a factorial bound. The number of iterations grows only logarithmically with inverse requested accuracy at fixed model parameters. The finite-time stability factors remain visible in the theorem. [Feedback reconstruction](feedback-reconstruction.md) gives the complete error and conditional-confidence calculation.

## From ordinary bits to a quantum circuit

The input contains finite programs for the actual coefficients and generator observables, together with certified bounds. The latter are genuine input promises. Smoothness by itself is not an efficient derivative-evaluation oracle.

Uniform bits generate approximate Gaussian increments using an explicit fixed-clock transformation. The analysis controls normal tails, arithmetic errors, terminal-time discretization and the failure of separately quantized Gaussian marginals to telescope exactly. Known mean-curve values at deterministic grid times are computed classically and embedded as circuit constants. There is no coherent lookup into an uncharged feedback table.

The direct reversible implementation retains path history and has the near inverse-linear gate bound. A second implementation uses a small-space pseudorandom generator and reversible checkpointing. The pseudorandom generator is required to preserve the first and second moments of the classical sample program; it is not claimed to fool a quantum computer. Quantum estimation then acts on the resulting explicit seed distribution. The construction trades additional gates for much smaller coherent workspace, while retaining the separately charged classical curve memory.

These claims depend on actual circuit-generation and arithmetic costs. [Finite precision and circuits](finite-precision-and-circuits.md) and [Memory and recomputation](memory-and-recomputation.md) state and prove that interface.

## What the advantage statement means

The same temporal reconstruction and couplings are available to classical computation. Ordinary MLMC gives $\widetilde O(\epsilon^{-2})$ work under the same analytic assumptions. Thus the comparison is already against a method that uses the new reconstruction, rather than against a naive particle simulation. Known particle antithetic methods give compatible strong baselines under their additional assumptions.

Explicit coefficient programs also permit structure-exploiting classical methods. One-dimensional Gauss quadrature can match the inverse-linear precision exponent. Analytic closures, PDE solvers, tensor representations, regression control variates and effective QMC may be better on other subclasses. A positive coarse-level variance only explains the cost of a sampling estimator; it proves no lower bound against these alternatives.

The supported contribution is therefore a complete quantum upper bound and a conditional route to advantage for a useful class of nonlinear stochastic tasks. Its significance is the removal of an additional feedback-construction cost, with preparation and output included. The [classical comparison](classical-comparison.md) identifies counterregimes, and the [polymer task](polymer-dynamics.md) supplies a concrete model and an initial classical stress test. An unconditional separation or a fault-tolerant hardware crossover would require further evidence.

## Proof and reproduction map

| Question | Authoritative account | What depends on it |
|---|---|---|
| What is the precise computational task? | [Model and results](model-and-results.md) | All guarantees and resource comparisons |
| Why do temporal details have $O(Hh)$ corrections? | [Dynkin details](dynkin-details.md) | The support-dependent multilevel allocation |
| Why do measured curves converge and remain admissible? | [Feedback reconstruction](feedback-reconstruction.md) | Complete nonlinear error and success probability |
| Are random bits, arithmetic and circuit construction affordable? | [Finite precision and circuits](finite-precision-and-circuits.md) | The logical gate upper bound |
| Can coherent workspace be reduced? | [Memory and recomputation](memory-and-recomputation.md) | The separate time-space tradeoff |
| What is a useful explicit model? | [Polymer dynamics](polymer-dynamics.md) | Physical input, observable and model limits |
| Which classical algorithms challenge the claim? | [Classical comparison](classical-comparison.md) | Scope of the conditional advantage |
| Which claims have been checked numerically? | [Verification](../verification/README.md) | Reproducible examples, not theorem certification |

