# Application: finite-bead polymer dynamics in a self-consistent Fourier field

The computational task is to predict a density relaxation amplitude or an intrachain scattering observable of a specified coarse-grained polymer model. The input is a finite list of physical parameters and an explicit Gaussian initial-state sampler. The output is one classical number with an additive error guarantee. This application motivates the algorithm; the numerical example below does not establish a physical quantum advantage.

The model belongs to the dynamic single-chain mean-field family developed by Fredrickson and Orland [FO14], with rigorous polymer McKean–Vlasov limits studied by E and Shen [ES13]. The finite Fourier interaction used here is a restriction we impose on the supplied model. It is not a claim that an arbitrary atomistic interaction can be compressed to fixed rank without error. Bibliographic details and the scope of these sources are in [Related work](../references/related-work.md).

## 1. Supplied model and its normalization

One stochastic state is an entire chain with $L\geq2$ beads,

```math
R=(r_1,\ldots,r_L)\in\mathbb R^{3L},\qquad d=3L.
```

Let $A$ be the open-chain graph Laplacian, acting independently on each spatial coordinate:

```math
(AR)_1=r_1-r_2,\quad
(AR)_j=2r_j-r_{j-1}-r_{j+1}\ (1<j<L),\quad
(AR)_L=r_L-r_{L-1}.
```

The input specifies $T>0$, spring coefficient $\kappa>0$, diffusion coefficient $D>0$, and $R_f$ wavevectors $q_a\in\mathbb R^3$ with real interaction coefficients $v_a$. Define the even pair potential

```math
V(z)=\sum_{a=1}^{R_f}v_a\cos(q_a\cdot z).
```

The number $R_f$ is fixed with respect to the requested numerical accuracy. Increasing it to approximate another interaction potential is a separate model-refinement problem whose error and cost are not included for free.

To specify the mean-field normalization, consider $N$ chains with bead coordinates $r_{i,j}$. The dimensionless energy is

```math
U_N=
\frac{\kappa}{2}\sum_{i=1}^N\sum_{j=1}^{L-1}
|r_{i,j+1}-r_{i,j}|^2
\;+
\frac{1}{2NL}\sum_{i,k=1}^N\sum_{j,\ell=1}^L
V(r_{i,j}-r_{k,\ell}).
```

Because $V$ is even, the interaction force on one bead is

```math
-\frac1{NL}\sum_{k=1}^N\sum_{\ell=1}^L
\nabla V(r_{i,j}-r_{k,\ell}).
```

Thus the interaction per bead remains of order one as $N$ increases. Physical density and interaction-strength conventions can be absorbed into the explicitly supplied $v_a$; they must not be silently varied with $N$ or $L$. Self terms have the same specified scaling and vanish from the tagged-chain mean-field limit at fixed $L$.

Define $K=2R_f$ bounded chain features

```math
g_{c,a}(R)=\frac1L\sum_{j=1}^L\cos(q_a\cdot r_j),\qquad
g_{s,a}(R)=\frac1L\sum_{j=1}^L\sin(q_a\cdot r_j),\qquad
m(t)=\mathbb E[g(R_t)].
```

The task is to solve the specified limiting equation

```math
dr_j(t)=b_j(R_t,m(t))\,dt+\sqrt{2D}\,dW_j(t),
```

```math
b_j(R,m)=-\kappa(AR)_j+
\sum_{a=1}^{R_f}v_aq_a
\left[m_{c,a}\sin(q_a\cdot r_j)-m_{s,a}\cos(q_a\cdot r_j)\right].
```

The $W_j$ are independent three-dimensional Brownian motions. This is an exact finite-rank representation of this chosen interaction. It models a spatially smoothed, weak mean field; it does not by itself reproduce hard-core exclusion, uncrossable chains or entanglement physics.

## 2. Obtainable initial state

An explicit nonequilibrium initial law avoids an unknown equilibrium-preparation problem. Draw a center $C\sim N(c_0,\Sigma_C)$ and independent bond vectors

```math
\zeta_j\sim N(0,(D/\kappa)I_3),\quad j=1,\ldots,L-1.
```

Set $S_1=0$, $S_j=\sum_{k<j}\zeta_k$, and

```math
r_j(0)=C+S_j-\frac1L\sum_{\ell=1}^L S_\ell.
```

This sampler uses $O(L)$ arithmetic and $O(L)$ Gaussian draws. It supplies the harmonic-chain Gaussian law for internal coordinates with a separately specified Gaussian center. It is not the interacting equilibrium law. Turning on the supplied $v_a$ at time zero is a well-defined quench.

Equivalently, in Rouse normal modes the nonzero eigenvalues are

```math
\lambda_p=4\sin^2\!\left(\frac{\pi p}{2L}\right),\qquad p=1,\ldots,L-1,
```

with independent coordinate variances $D/(\kappa\lambda_p)$. The center mode has covariance $L\Sigma_C$. The numerical example uses this equivalent representation to give QMC a useful coordinate ordering.

Both algorithms receive the same sampler. Quantum preparation consists of uniform random-bit registers followed by reversible Gaussian and sampling arithmetic; neither a preloaded probability table nor a prepared interacting quantum state is assumed. Gaussian truncation, precision, arithmetic, and inverse calls are part of the resource accounting. Coupled fine and coarse paths share the initial draw.

This input uses the [known random initial-law extension](dynkin-details.md), with initial mean features evaluated once to the specified accuracy and charged. The Gaussian law supplies the required moment and tail bounds.

The initial mean features can also be computed from Gaussian characteristic functions. If the marginal bead mean and covariance are $\mu_j,C_j$, then

```math
\mathbb E[e^{iq\cdot r_j(0)}]
=\exp(iq\cdot\mu_j-\tfrac12q^TC_jq).
```

For the sampler above, $\mu_j=c_0$ and

```math
C_j=\Sigma_C+a_jI_3,\qquad
a_j=\frac{D}{\kappa L^2}\bigl[S_2(j-1)+S_2(L-j)\bigr],\qquad
S_2(n)=\frac{n(n+1)(2n+1)}6.
```

All initial mean features therefore cost $O(LK)$ arithmetic. The initial scattering expectation is likewise explicit:

```math
\mathbb E\Phi_q(R_0)
=\frac1{L^2}\left[L+2\sum_{s=1}^{L-1}(L-s)
\exp\!\left(-\frac{D|q|^2s}{2\kappa}\right)\right].
```

A localized center gives nonzero density modes. By contrast, a spatially uniform center modulo the Fourier periods can make all retained mean modes zero. Without an external perturbation, that symmetry can reduce the dynamics to the exactly soluble Gaussian Rouse case and is excluded from any proposed advantage demonstration.

## 3. Useful output

The primary output is $m_{c,a}(T)$, or the corresponding complex pair $m_{c,a}(T)+i m_{s,a}(T)$. It describes the relaxation or growth of a prescribed spatial monomer-density mode after a quench. If this coordinate is already part of the reconstructed mean curve, reading its terminal value requires no new expectation-estimation problem.

A second useful scalar is the normalized intrachain scattering intensity

```math
S(q,T)=\mathbb E[\Phi_q(R_T)],\qquad
\Phi_q(R)=\left|\frac1L\sum_j e^{iq\cdot r_j}\right|^2\in[0,1].
```

This quantity reports chain conformation during the density evolution. It is different from the squared population density amplitude $|\mathbb E[L^{-1}\sum_j e^{iq\cdot r_j}]|^2$. Evaluating $\Phi_q$ costs $O(L)$; estimating it after mean-curve reconstruction requires the scalar terminal-estimation stage of the algorithm.

These are numerical predictions of the supplied coarse-grained model. Agreement with a particular material or experiment requires separately justified model parameters and model error; no such empirical calibration is claimed here.

The requested scalar meets the observable contract: in the normalized chain norm below, $\Phi_q$ is bounded by one and has Lipschitz constant at most $2|q|$.

## 4. Finite interaction rank does not close the state dynamics

The $K$ mean equations do not form a closed $K$-dimensional ODE. Applying the generator to a cosine feature already produces products of sine and cosine features at a single bead. For example, $\sin^2(q\cdot r_j)$ generates $\cos(2q\cdot r_j)$. Repeated differentiation creates further harmonics. The spring term additionally introduces mixed expectations such as

```math
\mathbb E\!\left[r_{j+1}\sin(q\cdot r_j)\right].
```

Consequently, the one-bead marginal equation requires conditional forces from neighboring beads. An ordinary three-dimensional one-bead Fokker–Planck equation is not an exact replacement for the joint $3L$-dimensional chain equation. A Gaussian initial law does not remain Gaussian under the nonlinear Fourier force.

This is a nonclosure statement, not a classical hardness argument. Dynamic SCFT, density-functional descriptions and covariance closures can approximate the observable well in suitable regimes. Fredrickson–Orland distinguish exact chain evolution from reductions that introduce constitutive or memory approximations [FO14]; single-chain mean-field particle methods already offer substantial classical computational savings [DM06]. Recent Gaussian closures are an additional serious alternative [G26]. Their error must be assessed for the chosen transient observable and tolerance.

## 5. Drift, generator and precision costs

For each bead, calculate its $R_f$ phases, their sine/cosine values, the spring force and the Fourier force. The full drift and all mean features require $O(LK)$ scalar operations, with no dense $3L$ by $3L$ matrix multiplication. Classical finite-rank particle aggregation likewise costs $O(NLK)$ per step, not $O(N^2)$.

The generator for a known mean curve is

```math
\mathcal A_m f=b(R,m)\cdot\nabla f+D\Delta f.
```

Its action on the required features has an explicit inexpensive form:

```math
\begin{aligned}
\mathcal A_mg_{c,a}
&=-\frac1L\sum_j\sin(q_a\cdot r_j)\,q_a\cdot b_j
  -D|q_a|^2g_{c,a},\\
\mathcal A_mg_{s,a}
&=\frac1L\sum_j\cos(q_a\cdot r_j)\,q_a\cdot b_j
  -D|q_a|^2g_{s,a}.
\end{aligned}
```

All $K$ generator features therefore cost $O(LK)$ too. No generic $O(d^2)$ Hessian computation is needed. For the scattering observable at one wavevector,

```math
\mathcal A_m\Phi_q
=2g_c\mathcal A_mg_c+2g_s\mathcal A_mg_s+\frac{2D|q|^2}{L}.
```

The $O(LK)$ count is an arithmetic count. The logical-gate cost additionally includes coordinate and coefficient bit lengths, range reduction and sine/cosine approximation, Gaussian arithmetic, controlled partial Euler steps, reversibility, and the inverse sample circuit. Increasing $|q_a|$, $|v_a|$, $D$, $T$ or $\kappa$ changes regularity and stability constants; those parameters are not suppressed by calling the result an inverse-linear precision bound.

## 6. Linear growth, dimension normalization and stiffness

The harmonic drift is globally Lipschitz with bounded positive-order derivatives but is not bounded in range. Hence $\mathcal A_mg$ need not be bounded either. The linear-growth version of the algorithm must use moment estimates and a variance-aware expectation estimator, with finite-bit truncation explicitly budgeted.

For fixed parameters and an admissible bounded Lipschitz mean curve, the Gaussian input supplies all finite moments. The relevant sufficient estimates are $L^4$ Euler strong order one, bounded fourth moments, and a weighted Lipschitz generator bound. In the normalized chain norm

```math
\|R\|_L^2=\frac1L\sum_j|r_j|^2,
```

a component of $g$ has Lipschitz constant at most $|q_a|$. Sensitivity of the force to the mean curve can also be bounded independently of $L$ in this norm, for fixed Fourier parameters. One may use, for example,

```math
\|b(R,m)-b(R,n)\|_L
\leq 2\Bigl(\sum_a|v_a|\,|q_a|\Bigr)\|m-n\|_\infty.
```

The state-dependent Fourier force is globally Lipschitz with a constant controlled by $\sum_a|v_a||q_a|^2$. The negative spring Laplacian is dissipative, so continuous-time difference estimates need not replace it by a destabilizing positive exponential. These observations help track dimension factors; they do not certify that every resulting constant is small or dimension-independent.

For fixed $m,n$, a sufficient generator estimate has the form

```math
|F(R,m)-F(S,n)|
\leq C(1+\|R\|_L+\|S\|_L)\|R-S\|_L+C\|m-n\|,
```

where $F=\mathcal A_mg_a$. Hölder's inequality and $L^4$ strong coupling then give the required $L^2$ correction estimate proportional to $Hh$ for a Dynkin detail of temporal width $H$ and simulation mesh $h$. Sharper bounds can use bond-force moments rather than the absolute center position. The finite-bit and statistical assumptions of the main theorem still apply.

Although $\|A\|\leq4$, a numerically stable explicit-Euler treatment of stiff springs should use a coarse-step restriction such as $4\kappa h\lesssim1$. The finite-horizon analytic estimate also permits larger coarse steps, at the cost of potentially very pessimistic constants. If a continuum-chain refinement makes $\kappa$ grow with $L$, that dependence must be charged. An exact or implicit harmonic propagator is a legitimate improvement for either classical or quantum computation, but it requires its own coupling and circuit accounting. The current application fixes the finite bead model and makes no continuum-limit cost claim.

### A conservative bead-count ledger

The dimension factors in this model need not introduce $e^{cLT}$ into the stability bound. To see this directly, use Euclidean coordinates $Y=R/\sqrt L$, equivalent to the normalized chain norm. Keep $K,T,\kappa,D,c_0,\Sigma_C,q_a,v_a$ fixed in this accounting. Sufficient constants obey

```math
\begin{aligned}
L_x&\leq4\kappa+2\sum_a|v_a||q_a|^2,\\
L_m&\leq2\sum_a|v_a||q_a|,\\
L_g&\leq\left(2\sum_a|q_a|^2\right)^{1/2},\\
\|\Sigma_Y\|_{\mathrm F}^2&=6D,\\
\|D_Y^2c_Y\|&\leq2\sqrt L\sum_a|v_a||q_a|^3.
\end{aligned}
```

Here $c_Y$ is the nonlinear part of the drift in the $Y$ variables. Thus the Picard/stability exponent involving $L_x+L_gL_m$ is independent of $L$ at fixed supplied parameters. The exact initial second moment is

```math
\mathbb E|Y_0|^2
=|c_0|^2+\operatorname{tr}\Sigma_C
+\frac D\kappa\frac{L^2-1}{2L}.
```

Gaussian fourth moments give $\|Y_0\|_{L^4}=O(\sqrt L)$. The linear-confining estimates consequently permit $M_4=O(\sqrt L)$, a mean-curve Lipschitz bound $O(\sqrt L)$ and an Euler $L^4$ error constant $C_{E,4}=O(\sqrt L)$. The direct weighted generator estimate displayed above has coefficients independent of $L$ in these coordinates. Combining it with these moments gives conservative bias/coupling constants $O(L)$. For the fixed-rank model and bit-efficient arithmetic, the resulting gate ledger is therefore at most of the form

```math
\widetilde O\!\left(C(T,\kappa,D,q,v,\Sigma_C,c_0,K)
\frac{L^2}{\epsilon}\right),
```

using the direct history implementation. This is a conservative scaling inference from the displayed constants, not a tight gate count or a dimension-independent theorem. It does not include fault-tolerant hardware overhead. Large stiffness, wavevector norms or long-time stability factors remain in $C$.

This accounting also does not define a hard growing-chain family. With fixed wavevectors the Gaussian coil broadens as $L$ grows, and retained density amplitudes can become exponentially small, making an absolute-error target trivial in some parameter ranges. A useful joint $L,\epsilon$ regime must specify nontrivial parameter and observable scaling separately. The stated finite-chain task requires no growing-chain limit.

## 7. Classical methods that can remove the advantage

The primary matched comparator is the direct classical Dynkin/Picard/MLMC construction in [Classical comparison](classical-comparison.md), with cost $\widetilde O(C_{\rm cl}(d,K)\epsilon^{-2})$ under the same analytic assumptions. Stronger smoothness permits the established antithetic particle/time analysis of Szpruch–Tse [ST21], also near quadratic after finite-rank aggregation is used.

For this particular polymer model, the following improvements must be credited:

- **Rouse control variates.** The interaction-free chain is Gaussian and its density and scattering expectations are available analytically. A synchronous reference can reduce variance substantially. At fixed nonzero interaction strength, a fixed reference generally changes a prefactor, not the precision exponent. Improving the reference with accuracy is allowed, with its construction cost counted.
- **QMC and MLQMC.** Brownian bridges, Rouse normal modes, scrambling, weighted constructions and suitable time integrators can reduce stochastic dimension and variance. High nominal state dimension is not evidence that these methods fail.
- **Moment or density closures.** For near-Gaussian states or slowly varying fields, a reduced equation may already meet the requested observable tolerance.
- **Tensor or regression representations.** Local chain interactions, a simple initial state and finite evolution time can make the joint law or backward conditional expectation compressible. The required ranks and approximation error are not bounded merely by finite interaction rank.
- **Temporal spectral reconstruction.** The displayed Fourier/harmonic example is analytic in its coefficients. If suitable uniform temporal analyticity bounds are available, a Chebyshev or graded spectral representation can replace the many temporal nodes used by the finite-regularity theorem. This possibility weakens any claim that the new reconstruction mechanism is necessary for this particular example. It does not provide a cheap high-dimensional stochastic expectation oracle.

No physical crossover is predicted in this repository. The conditional complexity result concerns a class specified by finite regularity and explicit access; it does not assert that every analytic physical specialization is classically difficult.

In particular, at fixed $d$ a sufficiently effective deterministic PDE/spectral solver for this analytic example may eventually beat the proposed precision exponent, while having a prohibitive dimension-dependent prefactor. The example therefore supports a finite high-dimensional benchmark and a restricted cost regime, not an asymptotic separation for a single fixed analytic physical model.

## 8. Illustrative classical experiment

The accompanying research experiment used $L=32$ ($d=96$), three nonorthogonal wavevectors, $T=1$, $\kappa=4$, $D=0.2$, and 32 Euler steps. There were twenty IID runs or independently scrambled Sobol runs at each particle count. Both implementations used the same Rouse-mode initial sampler, Brownian-bridge/Rouse-mode noise transform, drift, aggregation and observables. Both used exact reflection symmetry to set the sine means to zero, leaving three nontrivial cosine means.

| Particles | IID density variance | RQMC density variance | IID/RQMC density variance | IID/RQMC scattering variance |
|---:|---:|---:|---:|---:|
| 64 | $3.065\times10^{-4}$ | $8.433\times10^{-5}$ | 3.63 | 3.71 |
| 128 | $1.103\times10^{-4}$ | $2.012\times10^{-5}$ | 5.48 | 3.54 |
| 256 | $7.262\times10^{-5}$ | $3.600\times10^{-5}$ | 2.02 | 5.90 |
| 512 | $1.692\times10^{-5}$ | $2.544\times10^{-5}$ | 0.665 | 5.27 |
| 1024 | $1.899\times10^{-5}$ | $3.199\times10^{-6}$ | 5.94 | 3.72 |

RQMC usually gave a useful constant-factor improvement. No stable $P^{-2}$ variance regime was established over these five counts. Fits over this short range are approximately $P^{-1.07}$ versus $P^{-0.91}$ for density variance and $P^{-0.98}$ versus $P^{-1.04}$ for scattering variance. These noisy empirical slopes are not asymptotic complexity estimates.

A matched 2048-particle reference started at density amplitude $0.77466$, ending at $0.41420$ with interactions and $0.76457$ with the interaction switched off. Corresponding scattering estimates were $0.38696$ and $0.68002$. The test therefore exercises a material nonzero nonlinear field. It does not certify the exact target expectation.

The experiment is a reproducible implementation comparison at one time mesh, not a proof of hardness, a certification of particle/time bias, a comparison with all classical methods, or a quantum crossover calculation. In particular, Rouse control variates, MLQMC, tensor methods and closure solvers were not tested. The code and raw replicates must accompany any reuse of this table; it should not be presented as evidence of a quantum speedup.
