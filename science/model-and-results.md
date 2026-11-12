# Model, guarantees and resource scope

[Home](../README.md) · [Dynkin details](dynkin-details.md) · [Feedback reconstruction](feedback-reconstruction.md) · [Finite precision and circuits](finite-precision-and-circuits.md) · [Memory and recomputation](memory-and-recomputation.md)

## Computational task

The input describes a diffusion whose drift depends on finitely many expectations of its own state. The requested output is one classical expectation at a prescribed finite time:

```math
dX_t=a(X_t,m(t))\,dt+\Sigma\,dW_t,\qquad
m(t)=\mathbb E[g(X_t)],\qquad
X_0=x_0,\qquad
z=\mathbb E[f(X_T)].
```

Here $`X_t\in\mathbb R^d`$, $`W_t\in\mathbb R^q`$ is standard Brownian motion, $`\Sigma\in\mathbb R^{d\times q}`$ is constant, and $`g=(g_1,\ldots,g_K)`$. The horizon $`T\gt 0`$ and the moment count $`K`$ are supplied model parameters. The output specification is

```math
\Pr\lbrace \vert \widehat z-z\vert \leq\epsilon\rbrace \geq1-\alpha,
\qquad 0\lt \epsilon\lt 1,\quad 0\lt \alpha\lt \tfrac12.
```

The input consists of ordinary classical descriptions of $`x_0,\Sigma,a,g,f`$ and certified bounds/evaluation procedures specified below. It contains neither a prepared distribution of $`X_t`$ nor a table of its unknown moments.

The main result is an upper bound for this task. Its precision exponent is compared with particular classical algorithm families in the accompanying comparison; it is not a lower bound on every classical algorithm receiving the same coefficient descriptions.

## Norms and assumptions

Vector norms are Euclidean. Matrix norms are operator norms unless a Frobenius subscript is present. The norm of a multilinear derivative is its induced operator norm. Put

```math
D=\prod_{k=1}^K[-G_k,G_k],\qquad
G=\left(\sum_{k=1}^K G_k^2\right)^{1/2},
\qquad
\nu=\Vert \Sigma\Vert _{\mathrm F}.
```

The bounded-drift theorem uses the following assumptions.

| Quantity | Supplied bound or condition |
|---|---|
| Initial state | A known deterministic $`x_0\in\mathbb R^d`$. Its finite-bit approximation is charged. |
| Drift magnitude | $`\lvert a(x,m)\rvert\leq B_a`$ for every $`x\in\mathbb R^d`$ and $`m\in D`$. |
| State dependence | $`a(\cdot,m)`$ is twice continuously differentiable, $`\lVert D_xa\rVert\leq L_x`$, and $`\lVert D_x^2a\rVert\leq M_a`$, uniformly in $`m`$. |
| Moment dependence | $`\lvert a(x,m)-a(x,n)\rvert\leq L_m\lvert m-n\rvert`$, uniformly in $`x`$. The spatial derivatives are jointly continuous in $`(x,m)`$. Differentiability in $`m`$ is not required. |
| Mean observables | $`g_k\in C_b^3(\mathbb R^d)`$, $`\lvert g_k\rvert\leq G_k`$, and $`\lVert D^jg_k\rVert\leq D_{j,k}`$ for $`j=1,2,3`$. |
| Requested observable | $`\lvert f\rvert\leq F_0`$ and $`\lvert f(x)-f(y)\rvert\leq L_f\lvert x-y\rvert`$. |
| Arithmetic | Uniform finite-bit programs evaluate the initial moment vector, drift, requested observable and generator observables to certified errors. Their gate, workspace and classical compilation costs are charged as in [Finite precision and circuits](finite-precision-and-circuits.md). |

The generator observables are

```math
F_k(x,m)
=
a(x,m)\cdot\nabla g_k(x)
+\frac12\mathrm{tr}\!\left(
\Sigma\Sigma^\mathsf T\nabla^2g_k(x)
\right).
```

Evaluating $`F_k`$ includes evaluating the required derivatives of $`g_k`$. These evaluations are not free derivative oracles. A program for an explicit polynomial, trigonometric or other suitably approximable coefficient family can supply them; an arbitrary numerical black box for $`g_k`$ does not automatically do so.

The assumptions are dimension sensitive. For example, $`\nu^2`$ sums noise power over all $`q`$ driving coordinates, and derivative bounds can depend on $`d`$ and $`K`$. No dimension independence is inferred from the notation.

### Constants used in the proof

Define

```math
\begin{aligned}
L_g
&=\left(\sum_{k=1}^K D_{1,k}^2\right)^{1/2},\\
B_{F,k}
&=B_aD_{1,k}+\frac12\nu^2D_{2,k},\\
L_{F,x,k}
&=L_xD_{1,k}+B_aD_{2,k}+\frac12\nu^2D_{3,k},\\
L_{F,m,k}
&=L_mD_{1,k},\\
B_F
&=\left(\sum_{k=1}^K B_{F,k}^2\right)^{1/2},\\
L_*
&=2(B_F+1).
\end{aligned}
```

Thus $`\vert F_k\vert \leq B_{F,k}`$ and

```math
\vert F_k(x,m)-F_k(y,n)\vert 
\leq L_{F,x,k}\vert x-y\vert +L_{F,m,k}\vert m-n\vert .
```

For a moment curve $`v:[0,T]\to D`$ with Lipschitz constant at most $`L_*`$, let $`X^v`$ be the ordinary diffusion with drift $`a(x,v(t))`$. The continuous Euler approximation, including its final partial step, satisfies

```math
\sup_{0\leq t\leq T}
\Vert X_t^v-X_t^{v,h}\Vert _{L^2}
\leq C_Eh,
```

with the explicit sufficient constant

```math
\begin{aligned}
B_c&=L_mL_*+L_xB_a+\frac12M_a\nu^2,\\
B_M&=L_x\nu,\\
C_E&=e^{L_xT}\left(TB_c+\sqrt T\,B_M\right).
\end{aligned}
```

The proof, including Lipschitz time dependence, is given in [Dynkin details](dynkin-details.md#additive-noise-euler-error-with-an-explicit-constant).

For later resource bounds put

```math
\begin{aligned}
B_k&=L_{F,x,k}C_E+L_{F,m,k}L_*,\\
S_k&=3L_{F,x,k}C_E+L_{F,m,k}L_*,\\
\overline S_k&=\max\lbrace B_{F,k}/T,S_k\rbrace .
\end{aligned}
```

These constants account for bias, coupled second moments and the initial sampling level, respectively.

## Main theorem

**Theorem 1: finite-moment feedback with additive noise.** Under the assumptions above, there is an algorithm that returns $`\widehat z`$ with the specified additive error and success probability. It uses classical approximate Picard iteration, a dyadic representation of the mean curve, random-time generator estimators, and variance-aware quantum multilevel estimation.

At fixed model parameters, when the supplied evaluators and classical circuit generators have costs polynomial in the working bit precision and state-range bit length as specified in the circuit input contract, its logical gate count has precision dependence

```math
\widetilde O(\epsilon^{-1}).
```

The suppressed factors are logarithms in inverse accuracy and failure probability, together with the explicitly displayed model, arithmetic and stability factors. This statement does not suppress a dimension-dependent derivative bound, a large diffusion norm, or a long-time stability factor as a universal constant.

A direct history implementation can require $`\widetilde O(\epsilon^{-1})`$ coherent workspace at fixed model size. [Memory and recomputation](memory-and-recomputation.md) gives a separate tradeoff: for a specified fixed $`\kappa\gt 0`$ and a suitable small-space classical sample program, gate dependence $`\widetilde O(\epsilon^{-1-\kappa})`$ can be combined with $`\mathrm{poly}(\log(1/\epsilon))`$ coherent workspace at fixed model parameters. The dependence on the actual sample program and its workspace is retained there. The known classical mean curve remains separately stored and charged.

### Explicit accuracy parameters

Let

```math
E_x(T)=
\begin{cases}
(e^{L_xT}-1)/L_x,&L_x\gt 0,\\
T,&L_x=0,
\end{cases}
\qquad
D_f=L_fL_mE_x(T),
\qquad
C_{\mathrm{out}}=\max\lbrace 1,D_f\rbrace .
```

Set

```math
A=L_gL_m,\qquad
C_P=e^{(L_x+A)T},\qquad
\zeta=\frac{\epsilon}{2C_{\mathrm{out}}},
```

and choose

```math
0\lt \eta\leq
\min\left\lbrace 
\frac{\zeta}{2C_P},\,
1,\,
2T(B_F+1)
\right\rbrace .
```

Take an integer $`R\geq1`$ such that

```math
2G\,e^{L_xT}\frac{(AT)^R}{R!}\leq\frac{\zeta}{2}.
```

Such an $`R`$ exists. For $`A\gt 0`$ and $`G\gt 0`$, a sufficient explicit choice is any integer satisfying

```math
R\geq
\max\left\lbrace 
1,\,
2eAT,\,
\log_2^+\!\left(\frac{4G e^{L_xT}}{\zeta}\right)
\right\rbrace ,
\qquad
\log_2^+(u)=\max\lbrace 0,\log_2u\rbrace .
```

If $`A=0`$ or $`G=0`$, the corresponding vanishing initial-error term is interpreted directly and $`R=1`$ is sufficient for that term.

For each approximate Picard map use

```math
J=
\max\left\lbrace 
1,\,
\left\lceil
\log_2\frac{2T(B_F+1)}{\eta}
\right\rceil
\right\rbrace ,
\qquad
h_{\mathrm{time}}=T2^{-J},
\qquad
\delta=\frac{\eta}{4\sqrt K(J+1)}.
```

There are exactly $`2^J`$ estimated coefficients per moment: one endpoint coefficient and $`2^J-1`$ midpoint details. Each coefficient is estimated to absolute error at most $`\delta`$. Its allowed ideal-circuit failure probability is

```math
p_{\mathrm{coef}}=\frac{\alpha}{4RK2^J}.
```

Compute the initial vector $`\mu_0=g(x_0)`$ to Euclidean error at most $`\eta/4`$, project it onto $`D`$, and reuse that same approximation as the constant term in every reconstructed curve. Its evaluation cost is charged once. Coefficient error, initial-vector error and temporal interpolation then give a uniform approximate-map error at most $`3\eta/4\lt \eta`$.

The terminal expectation, after the final moment curve has been computed, is estimated to error $`\epsilon/2`$ with ideal-circuit failure probability at most $`\alpha/4`$. The coefficient calls together use an ideal-circuit failure budget $`\alpha/4`$. A further budget $`\alpha/2`$ bounds the total change in output probabilities caused by synthesis of the required non-Boolean quantum gates, as proved in [Finite precision and circuits](finite-precision-and-circuits.md).

These choices give a final moment-curve error at most $`\zeta`$, so its effect on the requested output is at most $`D_f\zeta\leq\epsilon/2`$. The two output errors sum to at most $`\epsilon`$. The conditional union bound over the ideal coefficient and terminal calls, followed by the compiled-circuit probability bound, gives total failure probability at most $`\alpha`$. The complete argument is in [Feedback reconstruction](feedback-reconstruction.md).

## A resource bound with the model factors exposed

For a coefficient of moment $`k`$ on a support of width $`H\leq T`$, let

```math
L_{k,H}
=
\max\left\lbrace 
0,\,
\left\lceil
\log_2\frac{4HB_kT}{\delta}
\right\rceil
\right\rbrace .
```

When $`B_k=0`$, take $`L_{k,H}=0`$. This is a simulation refinement level, not a temporal interpolation level.

Suppose a binary sample circuit on $`n`$ Euler steps, including its inverse, controls and deterministic coefficient schedule, has cost at most $`C_s(b)(n+1)`$ at working precision $`b`$. The symbol $`C_s`$ includes the state dimension, noise dimension, moment count and actual coefficient/generator programs. Let $`\Lambda`$ be the maximum polylogarithmic overhead of the variance-aware mean estimator over the finite collection of calls, including their confidence amplification and range/precision parameters. Then one coefficient has the sufficient cost bound

```math
W_{k,H}
\leq
C\,C_s(b)\Lambda
\left[
1+
\frac{HT}{\delta}
\left\lbrace 
B_k+\overline S_k(L_{k,H}+1)^2
\right\rbrace 
\right],
```

for a universal constant $`C`$, with the fixed finite-bit error allocations in [Finite precision and circuits](finite-precision-and-circuits.md). The proof does not rescale a correction by its standard deviation as though that were a pointwise bound.

Put $`L_{\max}=\max_k L_{k,T}`$. Summing the actual coefficient supports gives, for one Picard map,

```math
W_{\Phi}
\leq
C\,C_s(b)\Lambda
\left[
K2^J+
\frac{T^2(J+1)}{\delta}
\sum_{k=1}^K
\left\lbrace 
B_k+\overline S_k(L_{\max}+1)^2
\right\rbrace 
\right].
```

This uses the exact fact that the total support width at each temporal scale is $`T`$. Since

```math
2^J\leq\frac{4T(B_F+1)}{\eta},
\qquad
\delta^{-1}=\frac{4\sqrt K(J+1)}{\eta},
```

the displayed model factors multiply $`\eta^{-1}`$ and logarithms.

For the requested observable, define

```math
B_f=L_fC_E,\qquad
\overline S_f=\max\lbrace F_0/T,3L_fC_E\rbrace ,
\qquad
L_f^{\mathrm{sim}}
=
\max\left\lbrace 
0,\left\lceil
\log_2\frac{8B_fT}{\epsilon}
\right\rceil
\right\rbrace ,
```

again taking level zero if $`B_f=0`$. A sufficient terminal cost is

```math
W_f
\leq
C\,C_s(b)\Lambda
\left[
1+
\frac{2T}{\epsilon}
\left\lbrace 
B_f+
\overline S_f(L_f^{\mathrm{sim}}+1)^2
\right\rbrace 
\right].
```

The total cost is $`RW_\Phi+W_f`$, plus input reading, the certified initial-moment evaluation, and the explicitly charged classical reconstruction/compilation. Computing and storing one curve takes $`O(K2^J)`$ values and at most $`O(KJ2^J)`$ elementary reconstruction operations. Evaluating the known curve at the deterministic times used by a sampler costs a logarithmic classical search per point or a monotone scan. The resulting constants are embedded in the sampler; no coherent random-access memory is assumed.

The finite-bit construction makes $`b`$ logarithmic in inverse requested error, failure probability and the actual model/count bounds under its uniform evaluation assumptions. If a supplied coefficient family does not admit that bit complexity, its cost must replace $`C_s(b)`$; the theorem does not turn an expensive coefficient evaluator into a cheap one.

## Linear-confining drift: a separate corollary

The bounded-drift assumption can be replaced by a specified linear-confining model

```math
a(x,m)=-Ax+c(x,m),\qquad A=A^\mathsf T\succeq0,
```

where $`c`$ is bounded and satisfies the analogous state-derivative and moment-Lipschitz assumptions. The observables remain bounded with the stated derivative bounds.

**Corollary 2: linear confinement.** The conclusion of Theorem 1 continues to hold with the moment and weighted-generator constants derived in [Dynkin details](dynkin-details.md#linear-confining-drift). In particular, a uniform fourth-moment bound and an $`L^4`$ Euler error bound replace the bounded generator argument. The precision exponent is unchanged.

The constants depend explicitly on $`\Vert A\Vert `$, $`\Sigma`$, $`x_0`$, $`T`$ and the nonlinear coefficient bounds. A finite-bead polymer or other finite-dimensional model can be an instance of this corollary after its coefficients and observable are checked. Neither a bead-continuum limit nor an uncharged stiff limit follows from it.

### Supplied random initial laws

Both results allow a known random initial state independent of the Brownian motion when its sampling and finite-bit approximation costs are charged, fine and coarse paths share the same initial draw, and its initial moment vector $`\mu_0=\mathbb E[g(X_0)]`$ has a certified evaluation procedure. The linear-confining case also requires the stated initial fourth-moment and approximation bounds. Under the same polynomial bit-complexity contract for these extra procedures, the precision exponent is unchanged. The detailed argument is in [Known random initial laws](dynkin-details.md#known-random-initial-laws); the [polymer example](polymer-dynamics.md) supplies an explicit Gaussian sampler and initial Fourier moments.

## Attribution and boundaries

The constituent tools have established antecedents. Quantum multilevel Monte Carlo is developed by An, Linden, Liu, Montanaro, Shao and Wang in [Quantum 5, 481 (2021), Theorem 2; Theorem 3 specializes to SDEs](https://quantum-journal.org/papers/q-2021-06-24-481/). Miyamoto's [McKean–Vlasov algorithm, EPJ Quantum Technology 13, 3 (2026)](https://link.springer.com/article/10.1140/epjqt/s40507-025-00453-z) formulates the finite-moment feedback task and identifies the cost of repeatedly rebuilding paths in sequential time updates. Classical iterative multilevel particle approximation is developed by Szpruch, Tan and Tse in [Annals of Applied Probability 29 (2019)](https://arxiv.org/abs/1706.00907).

The contribution here is the complete feedback construction: represent temporal details by a randomly sampled generator value, retain the support-width factor in both bias and coupled second moments, use a different finest simulation grid for each support width, and reconstruct the unknown feedback curve with controlled error and regularity. Dynkin's identity, Picard iteration, dyadic interpolation and quantum mean estimation are credited components of that construction.

The result concerns fixed finite-time input models under the stated controls and accounting. Multiplicative diffusion, unrestricted long-time behavior, relative-error rare-event outputs, finite physical hardware crossover, and unconditional classical lower bounds are outside the theorem.

