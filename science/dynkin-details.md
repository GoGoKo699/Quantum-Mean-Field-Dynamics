# Random-time Dynkin details

[Home](../README.md) · [Model and results](model-and-results.md) · [Feedback reconstruction](feedback-reconstruction.md) · [Finite precision and circuits](finite-precision-and-circuits.md)

This page proves the sampling estimates used by the feedback algorithm. All norms and input constants have the meanings in [Model and results](model-and-results.md). The moment curve $`v:[0,T]\to D`$ is fixed throughout each statement. In the algorithm it is the classical output of a preceding iteration; the bounds are uniform over all such curves with Lipschitz constant at most $`L_*`$.

## The generator and a temporal detail

For the ordinary diffusion

```math
dX_t^v=a(X_t^v,v(t))\,dt+\Sigma\,dW_t,\qquad X_0^v=x_0,
```

Itô's formula and boundedness of the derivatives give

```math
g_k(X_t^v)
=
g_k(x_0)
+
\int_0^t F_k(X_s^v,v(s))\,ds
+
\int_0^t \nabla g_k(X_s^v)^\mathsf T\Sigma\,dW_s.
```

The stochastic integral is square-integrable and has zero mean. Thus, writing $`\Phi_k(v)(t)=\mathbb E[g_k(X_t^v)]`$,

```math
\Phi_k(v)(t)
=
g_k(x_0)+\int_0^t
\mathbb E[F_k(X_s^v,v(s))]\,ds.
```

In particular, $`\Phi(v)`$ is Lipschitz with constant at most $`B_F`$, independently of the derivative of $`v`$.

Let $`I=[u,u+H]\subset[0,T]`$ and let $`u_{1/2}=u+H/2`$. Define its midpoint detail by

```math
c_{k,I}(v)
=
\Phi_k(v)(u_{1/2})
-\frac{\Phi_k(v)(u)+\Phi_k(v)(u+H)}2.
```

Let $`\tau`$ be independent of the Brownian motion and uniform on $`I`$, and put

```math
s_I(t)=
\begin{cases}
+1,&u\leq t\lt u_{1/2},\\
-1,&u_{1/2}\leq t\leq u+H.
\end{cases}
```

**Lemma 1: random-time representation.**

```math
c_{k,I}(v)=\mathbb E[Y_{k,I}],
\qquad
Y_{k,I}
=
\frac H2\,s_I(\tau)F_k(X_\tau^v,v(\tau)).
```

**Proof.** Subtract the three generator-integral formulas above. The portion before $`u`$ cancels, leaving

```math
c_{k,I}(v)
=
\frac12\int_u^{u+H/2}\mathbb E[F_k(X_s^v,v(s))]\,ds
-
\frac12\int_{u+H/2}^{u+H}\mathbb E[F_k(X_s^v,v(s))]\,ds.
```

Integration against the uniform law on $`I`$ gives the claimed expectation. All exchanges are justified by the bounded generator. ∎

The endpoint coefficient has the related representation

```math
\Phi_k(v)(T)-g_k(x_0)
=
\mathbb E\!\left[T F_k(X_\tau^v,v(\tau))\right],
\qquad \tau\sim\mathrm{Unif}[0,T].
```

For common bounds, both cases can be written as $`Y=Hw(\tau)F_k(X_\tau^v,v(\tau))`$, with $`\vert w\vert \leq1`$. For a midpoint detail, $`w=s_I/2`$; for the endpoint coefficient, $`H=T`$ and $`w=1`$.

The sampled variable is a generator value at one random time. The estimator does not replace the integral by a pathwise quadrature sum. This distinction is responsible for the support-width factor retained below.

## Additive-noise Euler error with an explicit constant

For a uniform mesh $`h=T/n`$, let $`\pi_h(t)=h\lfloor t/h\rfloor`$ for $`t\lt T`$, with the endpoint convention chosen consistently. The continuous Euler process is

```math
X_t^{v,h}
=
x_0+
\int_0^t a(X_{\pi_h(s)}^{v,h},v(\pi_h(s)))\,ds
+\Sigma W_t.
```

At an off-grid terminal time $`t`$, this definition includes the partial step

```math
X_t^{v,h}
=
X_{\pi_h(t)}^{v,h}
+
a(X_{\pi_h(t)}^{v,h},v(\pi_h(t)))(t-\pi_h(t))
+
\Sigma(W_t-W_{\pi_h(t)}).
```

The driving Brownian motion is the same as for $`X^v`$.

**Lemma 2: uniform fixed-time strong error.** If $`v`$ is $`L_*`$-Lipschitz, then

```math
\sup_{0\leq t\leq T}
\Vert X_t^v-X_t^{v,h}\Vert _{L^2}
\leq C_Eh,
```

where

```math
C_E=e^{L_xT}
\left[
T\left(L_mL_*+L_xB_a+\frac12M_a\nu^2\right)
+\sqrt T\,L_x\nu
\right].
```

This is a supremum of fixed-time $`L^2`$ errors. A pathwise supremum inside the expectation is not asserted or needed.

**Proof.** Write $`b(t,x)=a(x,v(t))`$. It is $`L_x`$-Lipschitz in $`x`$ and $`L_mL_*`$-Lipschitz in $`t`$. First suppose it is smooth in time. Along the exact process,

```math
db(t,X_t^v)
=
c(t,X_t^v)\,dt+M(t,X_t^v)\,dW_t,
```

with

```math
\begin{aligned}
c&=\partial_t b+D_xb\,b
+\frac12\sum_{j=1}^qD_x^2b[\Sigma e_j,\Sigma e_j],\\
M&=D_xb\,\Sigma.
\end{aligned}
```

Therefore $`\vert c\vert \leq B_c`$ and $`\Vert M\Vert _{\mathrm F}\leq B_M`$, where

```math
B_c=L_mL_*+L_xB_a+\tfrac12M_a\nu^2,
\qquad
B_M=L_x\nu.
```

Decompose the error as

```math
X_t^v-X_t^{v,h}
=
\int_0^t
\left[
b(\pi_h(s),X_{\pi_h(s)}^v)
-b(\pi_h(s),X_{\pi_h(s)}^{v,h})
\right]ds
+
R_t,
```

where

```math
R_t
=
\int_0^t
\left[
b(s,X_s^v)-b(\pi_h(s),X_{\pi_h(s)}^v)
\right]ds.
```

The drift part of $`R_t`$ has norm at most $`B_cTh`$. By stochastic Fubini, its martingale part is

```math
\int_0^t k_{t,h}(r)M(r,X_r^v)\,dW_r,
\qquad
k_{t,h}(r)=\min\lbrace t,\pi_h(r)+h\rbrace -r.
```

The deterministic kernel satisfies $`0\leq k_{t,h}\leq h`$. Itô isometry consequently bounds this martingale's $`L^2`$ norm by $`B_Mh\sqrt T`$. Hence

```math
\Vert R_t\Vert _{L^2}\leq h(B_cT+B_M\sqrt T)
```

uniformly in $`t`$. Let $`e(t)=\sup_{u\leq t}\Vert X_u^v-X_u^{v,h}\Vert _{L^2}`$. Minkowski's inequality gives

```math
e(t)
\leq
L_x\int_0^t e(s)\,ds+h(B_cT+B_M\sqrt T).
```

Gronwall's inequality proves the displayed constant.

For a merely Lipschitz curve $`v`$, the coefficient $`b(t,x)`$ is Lipschitz in time. Extend it constantly beyond the interval and mollify only in time. The mollifications preserve the same time-Lipschitz, state-derivative and magnitude bounds. The preceding argument applies uniformly. Their coefficients converge uniformly to $`b`$; the exact and Euler solutions converge in fixed-time $`L^2`$ by the usual Lipschitz stability inequality. Passing to the limit proves the same estimate without differentiability in $`m`$ or a twice differentiable mean curve. ∎

If $`\tau`$ is independent and uniform on any subinterval $`I`$, then

```math
\Vert X_\tau^v-X_\tau^{v,h}\Vert _{L^2}^2
=
\frac1H\int_I
\mathbb E\vert X_t^v-X_t^{v,h}\vert ^2\,dt
\leq C_E^2h^2.
```

This is the reason a random terminal time does not change the strong-error order.

## Coupled random-time estimators

Define the implementable ideal-real-arithmetic variable

```math
Y_{k,I,h}
=
Hw(\tau)
F_k\!\left(X_\tau^{v,h},v(\pi_h(\tau))\right).
```

The endpoint moment value is frozen at the last grid time, just as the Euler drift is. In a circuit, that value is one of the known classical constants belonging to deterministic grid steps.

**Lemma 3: support-scaled bias and second moment.** For every $`h=T/n`$,

```math
\Vert Y_{k,I,h}-Y_{k,I}\Vert _{L^2}\leq HB_kh,
\qquad
\vert \mathbb E Y_{k,I,h}-\mathbb E Y_{k,I}\vert \leq HB_kh,
```

where $`B_k=L_{F,x,k}C_E+L_{F,m,k}L_*`$. If $`h`$ and $`2h`$ are nested meshes and both approximations use the same $`\tau`$ and the same Brownian motion, then

```math
\Vert Y_{k,I,h}-Y_{k,I,2h}\Vert _{L^2}
\leq HS_kh,
\qquad
S_k=3L_{F,x,k}C_E+L_{F,m,k}L_*.
```

Thus the coupled correction has second moment at most $`H^2S_k^2h^2`$. This is stronger than a variance bound.

**Proof.** The generator Lipschitz inequality and Lemma 2 give

```math
\begin{aligned}
\Vert Y_{k,I,h}-Y_{k,I}\Vert _{L^2}
&\leq
H\left[
L_{F,x,k}\Vert X_\tau^{v,h}-X_\tau^v\Vert _{L^2}
+
L_{F,m,k}L_*h
\right]\\
&\leq HB_kh.
\end{aligned}
```

The expectation bound follows from Cauchy–Schwarz.

For coupled grids, the triangle inequality through the exact process gives

```math
\Vert X_\tau^{v,h}-X_\tau^{v,2h}\Vert _{L^2}
\leq C_E(h+2h)=3C_Eh.
```

The two frozen grid times differ by at most $`h`$, so their moment values differ by at most $`L_*h`$. Apply the same generator Lipschitz inequality. ∎

The small quantity $`HS_kh`$ is an $`L^2`$ bound, not a pointwise range bound. In the bounded-drift model the correction's pointwise range is only bounded by $`2HB_{F,k}`$. A variance-aware or second-moment-aware mean estimator is essential; directly amplitude-encoding a correction after division by $`Hh`$ would generally be invalid.

### A coupled partial step

One can generate the fine Brownian increments up to the common $`\tau`$. Full fine increments have covariance $`hI_q`$; the last partial increment has covariance $`(\tau-\pi_h(\tau))I_q`$. A coarse full or partial increment is the sum of its constituent fine increments. Conditional on $`\tau`$, these sums have exactly the Brownian covariance and independence required by the coarse path. Both paths can be simulated in a fixed deterministic loop, masking their updates after $`\tau`$.

The support interval $`I`$ need not be aligned to the simulation mesh, and the simulation mesh need not resolve $`H`$. The final partial increment represents the chosen time separately. Simulation error remains $`O(h)`$, while the payoff multiplies it by $`H`$.

The finite-bit version, including independent Gaussian seeds, partial-time arithmetic, controls, workspace and inverse calls, is proved in [Finite precision and circuits](finite-precision-and-circuits.md).

## Quantum multilevel cost for one coefficient

Let $`h_\ell=T2^{-\ell}`$. For a support width $`H`$ define

```math
Z_0=Y_{k,I,h_0},
\qquad
Z_\ell=Y_{k,I,h_\ell}-Y_{k,I,h_{\ell-1}}
\quad(\ell\geq1).
```

The expectation telescopes:

```math
\sum_{\ell=0}^L\mathbb E Z_\ell
=
\mathbb E Y_{k,I,h_L}.
```

The different $`Z_\ell`$ need not be sampled jointly. Only the two paths inside each correction must have the stated coupling.

Put $`\overline S_k=\max\lbrace B_{F,k}/T,S_k\rbrace `$. Then, for all levels including zero,

```math
\Vert Z_\ell\Vert _{L^2}
\leq H\overline S_kh_\ell.
```

For coefficient tolerance $`\delta\gt 0`$, choose

```math
L=
\max\left\lbrace 
0,\left\lceil
\log_2\frac{4HB_kT}{\delta}
\right\rceil
\right\rbrace ,
```

with $`L=0`$ when $`B_k=0`$. Lemma 3 bounds the terminal discretization bias by $`\delta/4`$.

Estimate each level mean to statistical error at most $`\delta/[4(L+1)]`$ and ideal-circuit failure probability at most $`p_{\mathrm{coef}}/(L+1)`$. The sum of the statistical errors is at most $`\delta/4`$, and the union bound gives ideal coefficient failure probability at most $`p_{\mathrm{coef}}`$. Give the finite-bit implementation a correction perturbation budget at most $`\delta/[16(L+1)]`$, as detailed in the circuit account. These errors sum to at most $`\delta/16`$. The total coefficient error is therefore strictly below $`\delta`$, leaving slack for the specified rounding and optional pseudorandom-generator allocations. Gate synthesis has a separate global output-probability budget.

A second-moment-aware quantum mean estimator uses, up to its stated logarithmic factors,

```math
O\!\left(
1+\frac{H\overline S_kh_\ell(L+1)}{\delta}
\right)
```

calls at level $`\ell`$. The fixed finite-bit perturbation enlarges this by at most a constant in the same allocation. A call costs $`C_s(b)(1+2^\ell)`$, including the complete prefix from time zero. Hence

```math
W_{k,H}
\leq
C\,C_s(b)\Lambda
\left[
1+
\frac{HT}{\delta}
\left\lbrace 
B_k+\overline S_k(L+1)^2
\right\rbrace 
\right].
```

Indeed, $`\sum_{\ell=0}^L2^\ell=O(1+HB_kT/\delta)`$, while $`2^\ell h_\ell=T`$ makes every variance-dependent summand proportional to $`H\overline S_kT(L+1)/\delta`$.

The required finest simulation grid is coefficient dependent:

```math
h_L=O(\delta/H)
```

when refinement is necessary. Requiring every coefficient to use the finest grid needed for a width-$`T`$ coefficient would discard this advantage.

Quantum multilevel estimation is an established method; see [An et al., Quantum 5, 481 (2021)](https://quantum-journal.org/papers/q-2021-06-24-481/). The second-moment mean-estimation interface follows the bounded-$`L^2`$ route in [Montanaro, arXiv:1504.06987v3, Lemma 4 and Algorithm 2](https://arxiv.org/abs/1504.06987v3). The circuit account states and proves the concrete interface used here, including its confidence and finite-bit assumptions.

## Linear-confining drift

This section proves the changes needed when

```math
a(x,m)=-Ax+c(x,m),\qquad A=A^\mathsf T\succeq0.
```

Assume

```math
\vert c\vert \leq B_0,\qquad
\Vert D_xc\Vert \leq L_c,\qquad
\Vert D_x^2c\Vert \leq M_c,\qquad
\vert c(x,m)-c(x,n)\vert \leq L_m\vert m-n\vert .
```

Write $`a_0=\Vert A\Vert `$, $`L_x=a_0+L_c`$ and keep the same bounded $`C^3`$ moment observables. The deterministic initial condition can be replaced by an independent known initial law with supplied fourth-moment bound

```math
M_{0,4}=\Vert X_0\Vert _{L^4},
```

as explained below.

### Fourth moments

The exact variation-of-constants formula is

```math
X_t^v
=
e^{-At}X_0
+
\int_0^t e^{-A(t-s)}c(X_s^v,v(s))\,ds
+
\int_0^t e^{-A(t-s)}\Sigma\,dW_s.
```

Since $`A\succeq0`$, its semigroup is a contraction. A centered Gaussian vector $`Z`$ with covariance $`C`$ satisfies

```math
\mathbb E\vert Z\vert ^4=(\mathrm{tr}C)^2+2\mathrm{tr}(C^2)
\leq3(\mathrm{tr}C)^2.
```

Minkowski's inequality therefore gives the uniform bound

```math
\sup_{t\leq T}\Vert X_t^v\Vert _{L^4}
\leq
M_4
:=
M_{0,4}+B_0T+3^{1/4}\nu\sqrt T.
```

This bound is uniform over all curves $`v`$ in the permitted moment range.

The generator obeys

```math
\vert F_k(x,m)\vert \leq a_k\vert x\vert +b_k,
\qquad
a_k=a_0D_{1,k},
\qquad
b_k=B_0D_{1,k}+\tfrac12\nu^2D_{2,k}.
```

Consequently $`\Phi(v)`$ has Lipschitz constant at most

```math
B_F^{\mathrm{lin}}
=
\left[
\sum_{k=1}^K(a_kM_4+b_k)^2
\right]^{1/2}.
```

Use $`L_*=2(B_F^{\mathrm{lin}}+1)`$ for this model.

### Fourth-moment Euler error

Apply the proof of Lemma 2 in $`L^4`$. Along the exact process, the drift term in the Itô expansion of $`a(X_t^v,v(t))`$ has $`L^4`$ norm at most

```math
B_{c,4}
=
L_mL_*+
L_x(a_0M_4+B_0)
+\tfrac12M_c\nu^2.
```

The martingale coefficient remains pointwise bounded in Frobenius norm by $`B_M=L_x\nu`$.

For a vector stochastic integral $`N_t=\int_0^t H_s\,dW_s`$ with $`\Vert H_s\Vert _{\mathrm F}\leq r(s)`$ for a deterministic function $`r`$, Itô's formula for $`\vert N_t\vert ^4`$ gives

```math
\mathbb E\vert N_t\vert ^4
\leq
3\left(\int_0^t r(s)^2\,ds\right)^2.
```

To see this, its drift is at most $`6r(s)^2\mathbb E\vert N_s\vert ^2`$, and $`\mathbb E\vert N_s\vert ^2\leq\int_0^s r(u)^2du`$. Integrating yields the displayed inequality. Applying it with $`r(s)=B_Mk_{t,h}(s)`$ bounds the local martingale remainder in $`L^4`$ by $`3^{1/4}B_Mh\sqrt T`$. The same Gronwall argument yields

```math
\sup_{t\leq T}
\Vert X_t^v-X_t^{v,h}\Vert _{L^4}
\leq C_{E,4}h,
```

where

```math
C_{E,4}
=
e^{L_xT}
\left(
TB_{c,4}+3^{1/4}\sqrt T\,B_M
\right).
```

It follows that all used Euler paths with $`h\leq T`$ obey

```math
\sup_{t\leq T}\Vert X_t^{v,h}\Vert _{L^4}
\leq M_4^\star:=M_4+C_{E,4}T.
```

This is a sufficient finite-horizon bound without an unstated Euler stability restriction. It can be pessimistic for stiff $`A`$; its explicit $`a_0`$ and exponential factors remain charged. A stiff-limit improvement would require a separately analyzed integrator.

### Weighted generator estimate

Set

```math
\ell_{0,k}
=
L_xD_{1,k}+B_0D_{2,k}+\tfrac12\nu^2D_{3,k},
\qquad
\ell_{1,k}=a_0D_{2,k}.
```

A direct subtraction gives

```math
\vert F_k(x,m)-F_k(y,n)\vert 
\leq
(\ell_{0,k}+\ell_{1,k}\vert y\vert )\vert x-y\vert 
+
L_mD_{1,k}\vert m-n\vert .
```

Hölder's inequality consequently gives, for the coupled states in this construction,

```math
\Vert F_k(X,m)-F_k(Y,n)\Vert _{L^2}
\leq
(\ell_{0,k}+\ell_{1,k}M_4^\star)
\Vert X-Y\Vert _{L^4}
+
L_mD_{1,k}\vert m-n\vert .
```

Thus Lemma 3 remains valid after the substitutions

```math
\begin{aligned}
L_{F,x,k}&\longleftarrow
\ell_{0,k}+\ell_{1,k}M_4^\star,\\
L_{F,m,k}&\longleftarrow L_mD_{1,k},\\
C_E&\longleftarrow C_{E,4},\\
B_{F,k}\ \text{in the base sampling bound}
&\longleftarrow a_kM_4^\star+b_k,\\
B_F\ \text{in the temporal Lipschitz bound}
&\longleftarrow B_F^{\mathrm{lin}}.
\end{aligned}
```

The distinction between the base sampling bound and the exact mean-curve Lipschitz bound avoids a circular definition of $`L_*`$. Every constant is determined from the input bounds in the order displayed above.

The base variable now has a bounded second moment, rather than a bounded range before truncation. The second-moment quantum estimator and the finite-bit tail allocation apply. This proves the linear-confining corollary's analytic estimates.

## Known random initial laws

The same statements hold when a supplied random initial state $`X_0`$ is independent of $`W`$, provided:

1. Its sampling and finite-bit approximation costs are included in every path call.
2. The fine and coarse paths use the same initial draw.
3. The initial moment vector $`\mu_0=\mathbb E[g(X_0)]`$ is supplied by a certified evaluation procedure to the allocated accuracy.
4. For linear confinement, the fourth-moment and finite-bit tail/approximation bounds required above are supplied.

Replace $`g(x_0)`$ by $`\mu_0`$ in all identities. Itô's formula and the martingale expectation argument are unchanged. The Euler error starts from zero because both paths use the same initial draw. In the bounded-drift model its displayed constant is independent of the initial value; in the linear-confining model the displayed $`M_{0,4}`$ carries the initial-law dependence.

For a Gaussian initial law and Fourier observables, the characteristic function provides $`\mu_0`$ directly. Such an instance is developed in [Polymer dynamics](polymer-dynamics.md). An initial distribution with no accessible sampler or no accounted-for initial moments does not satisfy this input contract.

