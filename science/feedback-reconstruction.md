# Feedback reconstruction and end-to-end error

[Home](../README.md) · [Model and results](model-and-results.md) · [Dynkin details](dynkin-details.md) · [Finite precision and circuits](finite-precision-and-circuits.md)

This page completes the proof of Theorem 1. It gives an invariant class of classical mean curves, an approximate global Picard iteration, and a final scalar readout. The estimates hold uniformly over the measured outputs of earlier iterations. Quantum state preparation is never nested inside another quantum expectation estimator.

## The frozen-feedback map

For a deterministic curve $v:[0,T]\to D$, define

```math
dX_t^v=a(X_t^v,v(t))\,dt+\Sigma\,dW_t,\qquad
\Phi(v)(t)=\mathbb E[g(X_t^v)].
```

The known initial vector is $\mu_0=g(x_0)$. In the known-random-initial-law extension it is $\mu_0=\mathbb E[g(X_0)]$, and the same initial sample is used in each coupling.

The drift is globally Lipschitz in its state, so the frozen equation has a unique strong solution. Because $|g_k|\leq G_k$, the map $\Phi$ takes its values in the box $D$. The generator identity in [Dynkin details](dynkin-details.md#the-generator-and-a-temporal-detail) gives

```math
|\Phi(v)(t)-\Phi(v)(s)|\leq B_F|t-s|.
```

This bound does not involve the Lipschitz constant of $v$. In particular, $\Phi$ maps the class of continuous $D$-valued curves into a uniformly Lipschitz class.

### Stability without a short-time contraction assumption

Couple $X^v$ and $X^w$ with the same initial state and Brownian path. Additive noise cancels from their difference. The state and moment Lipschitz bounds give, pathwise,

```math
|X_t^v-X_t^w|
\leq
L_x\int_0^t|X_s^v-X_s^w|\,ds
+L_m\int_0^t|v(s)-w(s)|\,ds.
```

The inhomogeneous Gronwall inequality therefore implies

```math
|X_t^v-X_t^w|
\leq
L_m\int_0^t e^{L_x(t-s)}|v(s)-w(s)|\,ds.
```

Since $g$ is $L_g$-Lipschitz in Euclidean norm,

```math
|\Phi(v)(t)-\Phi(w)(t)|
\leq
A\int_0^t e^{L_x(t-s)}|v(s)-w(s)|\,ds,
\qquad A=L_gL_m.
```

Define the positive Volterra operator

```math
(\mathcal K u)(t)=A\int_0^t e^{L_x(t-s)}u(s)\,ds
```

on nonnegative functions. Repeated convolution gives, for every integer $r\geq1$,

```math
(\mathcal K^r u)(t)
=
\frac{A^r}{(r-1)!}
\int_0^t e^{L_x(t-s)}(t-s)^{r-1}u(s)\,ds.
```

This identity follows by induction: the inner integral over $s\leq u\leq t$ is $(t-s)^r/r$. Hence

```math
(\mathcal K^r\mathbf 1)(t)
\leq e^{L_xt}\frac{(At)^r}{r!}.
```

**Lemma 4: existence, uniqueness and inexact Picard stability.** The map $\Phi$ has a unique fixed point $m$ in the continuous $D$-valued curves. If $v_0$ is any $D$-valued continuous curve and

```math
\|v_{r+1}-\Phi(v_r)\|_\infty\leq\eta,
\qquad r=0,\ldots,R-1,
```

then

```math
\|v_R-m\|_\infty
\leq
2G\,e^{L_xT}\frac{(AT)^R}{R!}
+\eta e^{(L_x+A)T}.
```

**Proof.** For exact Picard iterates starting at $v_0$, the difference of the first two curves is at most $2G$ pointwise, since both lie in $D$. The stability inequality bounds successive differences after $r$ further maps by $2G\mathcal K^r\mathbf 1$. Their uniform norms have a summable factorial upper bound. The iterates therefore converge uniformly to a continuous $D$-valued curve. Continuity of $\Phi$ under its displayed integral bound makes the limit a fixed point.

Two fixed points have difference bounded by $\mathcal K^r(2G\mathbf1)$ for every $r$. The factorial bound tends to zero, proving uniqueness. The process $X^m$ then solves the specified McKean–Vlasov equation, and any other solution would produce a second fixed-point mean curve and then contradict strong uniqueness of the frozen equation.

For the inexact iteration, write $e_r(t)=|v_r(t)-m(t)|$. Positivity gives

```math
e_R
\leq
\mathcal K^R e_0
+\eta\sum_{j=0}^{R-1}\mathcal K^j\mathbf1.
```

Use $e_0\leq2G$ and

```math
\sum_{j=0}^{R-1}\mathcal K^j\mathbf1
\leq
1+e^{L_xT}(e^{AT}-1)
\leq e^{(L_x+A)T}.
```

This proves the result. ∎

The proof does not require $AT<1$ or a norm-contraction condition on $\Phi$. It does retain the displayed finite-time stability factors. Holding $T$ fixed is substantive when these factors are large.

## A dyadic representation of the whole mean curve

Fix a positive integer $J$, and let

```math
t_i=iT2^{-J},\qquad 0\leq i\leq2^J,\qquad
h_{\mathrm{time}}=T2^{-J}.
```

For temporal scale $j=0,\ldots,J-1$, define intervals

```math
I_{j,r}=[rH_j,(r+1)H_j],\qquad
H_j=T2^{-j},\qquad 0\leq r<2^j.
```

Let $\psi_{j,r}$ be the continuous piecewise-linear hat supported on $I_{j,r}$, equal to one at its midpoint and zero at its endpoints. Its nonzero slopes are $\pm2/H_j$.

For any continuous vector function $u$, its piecewise-linear interpolant on the finest grid has the exact finite identity

```math
\mathcal I_Ju(t)
=
u(0)+\frac tT\bigl(u(T)-u(0)\bigr)
+\sum_{j=0}^{J-1}\sum_{r=0}^{2^j-1}
\left[
u\!\left((r+\tfrac12)H_j\right)
-\frac{u(rH_j)+u((r+1)H_j)}2
\right]\psi_{j,r}(t).
```

To verify it, start with the line between the two endpoints. Adding the hats at scale $j$ corrects the value at each newly introduced midpoint and leaves all earlier grid values unchanged. After $J$ scales, the resulting function is the finest-grid interpolant.

For $u=\Phi(v)$, the endpoint difference and every midpoint detail have the random-time estimators already proved in [Dynkin details](dynkin-details.md). No list of pointwise values of $\Phi(v)$ is assumed available.

### Reconstruction error and a uniform slope bound

Use the parameters

```math
0<\eta\leq\min\{1,\,2T(B_F+1)\},
\qquad
J=\max\left\{1,\left\lceil
\log_2\frac{2T(B_F+1)}{\eta}
\right\rceil\right\},
\qquad
\delta=\frac{\eta}{4\sqrt K(J+1)}.
```

They imply

```math
\frac{\eta}{4(B_F+1)}
\leq h_{\mathrm{time}}
\leq\frac{\eta}{2(B_F+1)},
\qquad
2^J\leq\frac{4T(B_F+1)}{\eta}.
```

Compute a binary approximation $\widehat\mu_0\in D$ with

```math
|\widehat\mu_0-\mu_0|\leq\eta/4.
```

A certified approximation can be projected onto $D$ without increasing this error. It is reused in every iteration. The endpoint coefficient always estimates $\Phi(v)(T)-\mu_0$, through its generator representation; it does not subtract an independently noisy initial estimate.

Estimate each scalar endpoint coefficient and midpoint detail to error at most $\delta$. Let $\widetilde u$ be the Faber–Schauder reconstruction using those estimates and $\widehat\mu_0$. Finally project each finest-grid node of $\widetilde u$ onto $D$ and linearly interpolate the projected nodes. Denote the resulting curve by $\widehat\Phi(v)$.

**Lemma 5: stable reconstruction from approximate details.** If every scalar coefficient estimate satisfies its tolerance, then

```math
\|\widehat\Phi(v)-\Phi(v)\|_\infty\leq3\eta/4<\eta,
\qquad
\operatorname{Lip}(\widehat\Phi(v))\leq L_*=2(B_F+1).
```

The curve is $D$-valued.

**Proof.** At a fixed time, at most one hat at each scale is nonzero, and every hat lies in $[0,1]$. The endpoint factor $t/T$ also lies in $[0,1]$. Combining $K$ scalar errors in Euclidean norm gives

```math
\|\widetilde u-\mathcal I_J\Phi(v)\|_\infty
\leq
|\widehat\mu_0-\mu_0|+(J+1)\sqrt K\,\delta
\leq\eta/2.
```

For a $B_F$-Lipschitz vector function on a grid cell $[a,a+h]$, put $\theta=(t-a)/h$. Its interpolation error is at most

```math
(1-\theta)B_F(t-a)+\theta B_F(a+h-t)
=2B_Fh\theta(1-\theta)
\leq B_Fh/2.
```

Thus $\|\mathcal I_J\Phi(v)-\Phi(v)\|_\infty\leq\eta/4$.

At every grid node the exact vector $\Phi(v)(t_i)$ lies in $D$. Euclidean projection onto this closed convex box is nonexpansive, so projection cannot increase the node error. Inside a grid cell, convex interpolation of the two node errors remains at most their maximum. Adding the preceding interpolation error proves the total bound $3\eta/4$.

For regularity, the exact interpolant has Lipschitz constant at most $B_F$. The initial-vector error is a constant function and has zero slope. The error from the estimated coefficients has, on every open finest-grid cell, slope norm at most

```math
\sqrt K\,\delta
\left(\frac1T+\sum_{j=0}^{J-1}\frac2{H_j}\right)
=
\sqrt K\,\delta\left(\frac2{h_{\mathrm{time}}}-\frac1T\right)
\leq
\frac{2\sqrt K\,\delta}{h_{\mathrm{time}}}
\leq
\frac{2(B_F+1)}{J+1}
\leq B_F+1.
```

Therefore $\operatorname{Lip}(\widetilde u)\leq2B_F+1<L_*$. Projection of neighboring nodes cannot increase the distance between them. The projected linear interpolant consequently has no larger Lipschitz constant. It is $D$-valued because the box is convex. ∎

This proof only needs the Lipschitz regularity of $\Phi(v)$ and bounded $C^3$ mean observables. A second time derivative of the mean curve is not assumed. The price is $O(\eta^{-1})$ stored nodes at fixed model parameters.

The reconstruction itself can use exact binary arithmetic. Coefficient estimates and $\widehat\mu_0$ are binary rationals, and every hat value at a normalized dyadic node is a dyadic rational. Their finite sums therefore need no rounding error beyond the coefficient tolerances. Bounds for the clipping box and slope check can be supplied as certified rational bounds. The resulting bit length grows by only logarithmic factors; its cost is included in the classical compilation account.

## The complete algorithm

Use $C_{\mathrm{out}},\zeta,C_P,\eta,R,J,\delta$ from [Model and results](model-and-results.md#explicit-accuracy-parameters). In particular,

```math
\zeta=\frac{\epsilon}{2C_{\mathrm{out}}},
\qquad
\eta\leq\frac{\zeta}{2C_P},
\qquad
2G e^{L_xT}\frac{(AT)^R}{R!}\leq\zeta/2.
```

The classical controller performs the following steps.

1. Compute $\widehat\mu_0$ to Euclidean error $\eta/4$, project it onto $D$, and initialize $v_0(t)=\widehat\mu_0$.
2. For each $r=0,\ldots,R-1$, estimate the endpoint coefficient and all dyadic details of $\Phi(v_r)$, for all $K$ moments. Each scalar estimate has tolerance $\delta$ and ideal-circuit failure budget $\alpha/(4RK2^J)$.
3. Reconstruct and project the nodes as in Lemma 5, obtaining $v_{r+1}$. Check that the stored piecewise-linear curve has Lipschitz constant at most $L_*$. If this check fails, terminate with a failure flag.
4. With the final classical curve $v_R$, estimate $\mathbb E[f(X_T^{v_R})]$ to error $\epsilon/2$, with ideal-circuit failure budget $\alpha/4$, using the terminal multilevel estimator below.

The slope check is a finite classical check of adjacent stored nodes. On the event used in Lemma 5 it always passes. Its purpose is to ensure that every sampler actually called has an admissible input curve, including branches on which an earlier statistical estimate failed. A failure flag counts against the allowed probability of unsuccessful output.

There are $2^J$ coefficients per moment per iteration: one endpoint and

```math
\sum_{j=0}^{J-1}2^j=2^J-1
```

details. The initial-vector computation is deterministic under the certified evaluation contract and does not consume a statistical failure budget.

### Why adaptivity does not multiply failure probabilities

Condition on all outputs from the previous iterations and on the event that no slope check has failed. The next curve is then an ordinary deterministic $D$-valued $L_*$-Lipschitz input. All the estimates in the Dynkin and circuit analyses hold uniformly over that class. Therefore the stated failure bound applies conditionally to each coefficient call.

Enumerate the at most $RK2^J$ coefficient calls. In the ideal quantum circuit, the probability of a first failed coefficient estimate is at most the sum of its conditional failure bounds, which is $\alpha/4$. No independence assumption between the measured estimates is required. If all coefficient calls succeed, Lemma 5 ensures that all slope checks pass and all approximate-map errors are at most $\eta$.

Lemma 4 then yields

```math
\|v_R-m\|_\infty
\leq
2G e^{L_xT}\frac{(AT)^R}{R!}
+\eta C_P
\leq\zeta.
```

The ideal terminal estimator has failure probability at most $\alpha/4$, conditional on any admissible final curve. The ideal algorithm therefore has total failure probability at most $\alpha/2$. The finite gate implementation is synthesized so that the total variation distance between its output law and that of the ideal circuit is at most $\alpha/2$. The hybrid bound establishing this, including all adaptive branches, is given in [Finite precision and circuits](finite-precision-and-circuits.md). The implemented algorithm consequently fails with probability at most $\alpha$.

## Terminal readout and its resource cost

The stability estimate for the state directly gives

```math
\left|
\mathbb E[f(X_T^v)]-\mathbb E[f(X_T^m)]
\right|
\leq
L_fL_m\int_0^T e^{L_x(T-s)}|v(s)-m(s)|\,ds
\leq D_f\|v-m\|_\infty,
```

where

```math
D_f=L_fL_mE_x(T),\qquad
E_x(T)=
\begin{cases}
(e^{L_xT}-1)/L_x,&L_x>0,\\
T,&L_x=0.
\end{cases}
```

Consequently the feedback approximation contributes at most $D_f\zeta\leq\epsilon/2$.

For a fixed admissible curve $v$, use terminal Euler payoffs

```math
P_\ell=f(X_T^{v,h_\ell}),\qquad
h_\ell=T2^{-\ell}.
```

Use $P_0$ at level zero and $P_\ell-P_{\ell-1}$ at higher levels, coupled by the same Brownian path and initial draw. Lemma 2 and Lipschitz continuity of $f$ give

```math
|\mathbb E P_\ell-\mathbb E f(X_T^v)|
\leq B_fh_\ell,\qquad B_f=L_fC_E,
```

and

```math
\|P_\ell-P_{\ell-1}\|_{L^2}
\leq3L_fC_Eh_\ell,\qquad
\|P_0\|_{L^2}\leq F_0.
```

Thus $\overline S_f=\max\{F_0/T,3L_fC_E\}$ bounds all correction norms by $\overline S_fh_\ell$. With terminal tolerance $\epsilon/2$, take

```math
L_f^{\mathrm{sim}}
=
\max\left\{
0,\left\lceil\log_2\frac{8B_fT}{\epsilon}\right\rceil
\right\},
```

with level zero if $B_f=0$. Apply the same per-level statistical, precision and confidence allocations as for one coefficient. This gives terminal error at most $\epsilon/2$ and the sufficient cost

```math
W_f
\leq
C C_s(b)\Lambda
\left[
1+\frac{2T}{\epsilon}
\left\{
B_f+\overline S_f(L_f^{\mathrm{sim}}+1)^2
\right\}
\right].
```

Together with the feedback contribution, this proves the output error at most $\epsilon$ on the event of probability at least $1-\alpha$. Only a bounded Lipschitz terminal observable is needed; the algorithm never evaluates a generator or a derivative of $f$.

## Summing the cost of a whole feedback curve

The crucial sum is over actual support widths, including complete path prefixes in every sample. For a fixed moment $k$, the coefficient cost from the Dynkin analysis is

```math
W_{k,H}
\leq
C C_s(b)\Lambda
\left[
1+\frac{HT}{\delta}
\left\{
B_k+\overline S_k(L_{k,H}+1)^2
\right\}
\right].
```

At temporal scale $j$ there are $2^j$ intervals, each of width $H_j=T2^{-j}$, and hence

```math
\sum_{r=0}^{2^j-1}H_j=T.
```

The endpoint contributes one further width $T$. Since $L_{k,H}\leq L_{k,T}\leq L_{\max}$, one Picard map costs at most

```math
W_\Phi
\leq
C C_s(b)\Lambda
\left[
K2^J+
\frac{T^2(J+1)}{\delta}
\sum_{k=1}^K
\left\{
B_k+\overline S_k(L_{\max}+1)^2
\right\}
\right].
```

The first term counts the nonzero minimum cost of every individual coefficient estimate. It is necessary when the support is small enough that no simulation refinement is needed.

Substituting the choices of $J$ and $\delta$ yields one inverse power of $\eta$ times logarithms and the displayed model factors. At fixed model parameters, $\eta$ is proportional to $\epsilon$, while $R,J,L_{\max}$ and the confidence amplification cost grow at most logarithmically in inverse precision. The terminal cost has the same inverse-precision power. The complete cost is therefore

```math
RW_\Phi+W_f
=\widetilde O(\epsilon^{-1})
```

under the stated bit-program assumptions, with classical input reading, node reconstruction, memory and compilation included as specified in the resource pages.

This summation would fail if every temporal detail were forced to use the finest simulation grid required by the endpoint coefficient. The finest grid is chosen separately for each support: its bias is proportional to $Hh$, so a short support permits a coarser simulation at the same coefficient tolerance.

## Linear confinement and known random initial laws

For the linear-confining model, use $B_F^{\mathrm{lin}}$ in place of $B_F$ for temporal regularity, and the substituted second-moment and Euler constants in [Dynkin details](dynkin-details.md#linear-confining-drift) for coefficient and terminal costs. The state stability inequality still uses $L_x=\|A\|+L_c$ and $L_m$. Every argument on this page then applies unchanged.

For a supplied random initial law, use its certified initial moment vector $\mu_0$ and a common initial draw in each coupled path. Its sampling, approximation and initial-moment costs remain part of the input contract. The example in [Polymer dynamics](polymer-dynamics.md) has an explicit Gaussian initial law and explicit Fourier initial moments.

The method estimates finitely many coefficients of a classical feedback curve before estimating the requested scalar. It does not output the many-particle distribution, a trajectory ensemble, or a full quantum state of the nonlinear diffusion.
