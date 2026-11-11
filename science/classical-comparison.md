# Classical comparison under matched access

The relevant classical baseline is already nearly quadratic in the inverse error, not a naive cubic-cost particle simulation. Moreover, there are important subclasses where classical methods match the proposed inverse-linear quantum precision exponent. For example, Kloeden–Shardlow give a deterministic $O(\epsilon^{-1}|\log\epsilon|^3)$ method for an appropriately smooth, nondegenerate **one-dimensional** mean-field SDE with a suitably approximated initial distribution [KS17]. No quantum advantage is claimed for that regime.

The contribution studied here is a conditional upper bound for finite-rank, additive-noise mean-field dynamics with explicitly supplied finite-regularity coefficient circuits. It is compared with algorithms that have the same input information and return the same small classical output. It is not a classical lower bound and does not establish that an explicit physical instance defeats every classical method. Full source details and the distinctions between proved and conjectural results appear in [Related work](../references/related-work.md).

## 1. What is matched

Both sides receive the finite programs for the drift, mean observables and terminal payoff, the initial sampler, the diffusion matrix, and the required derivative and stability bounds. The classical side can inspect these programs and exploit their structure. It is not restricted to a black-box sampling oracle if the quantum side receives a succinct circuit. Finite-rank aggregation, known symmetries, Gaussian reference solutions, normal-mode transforms, multilevel coupling, higher-order integration when justified, and classical preprocessing are available to both sides.

The task is an additive-error estimate of $\mathbb E[f(X_T)]$ or a specified mean coordinate at $T$, with a stated success probability. A strong approximation of a random solution path, a full density reconstruction, and a relative-error rare-event estimate are different tasks. Complexity statements for those tasks are not transferred unchanged.

We track serial work, arithmetic or logical gates as appropriate. Parallel hardware, fault-tolerant overhead and implementation constants matter for a practical crossover, but no such crossover is inferred from a precision exponent. The fixed-$T$, fixed-$K$ asymptotics retain explicit dependence on dimension, regularity, stiffness, diffusion norms, initial moments and circuit costs.

## 2. A direct classical analogue under the same analytic assumptions

The classical side can use **the same** global Picard iteration, temporal interpolation, randomized Dynkin details and fine/coarse Euler coupling as the quantum algorithm. Replacing quantum mean estimation by ordinary multilevel Monte Carlo gives a useful matched upper bound without requiring the stronger measure-differentiability assumptions of a particle antithetic theorem.

The derivation below needs only the analytic estimates used by the main algorithm: additive-noise strong order one for a supplied Lipschitz mean curve; a Lipschitz, or appropriately weighted Lipschitz, generator payoff; and the stated moment bounds. In the bounded setting, the sufficient conditions are bounded drift derivatives through order two in the state variable, bounded mean-feature derivatives through order three, and the stated Lipschitz bounds in the mean variable. The linear-growth extension uses the corresponding $L^4$ estimates. No fourth derivative of the mean feature is needed for the Lipschitz temporal reconstruction.

### 2.1 One temporal detail

Fix a Picard input curve $v$. A Faber–Schauder detail supported on an interval of width $H$ is the expectation of the randomized Dynkin payoff

```math
Y_H=\frac H2\,s(\tau)F_v(\tau,X_\tau^v),
```

where $\tau$ is uniform on the support interval and $s$ is $+1$ on its first half and $-1$ on its second. Let $h_\ell\asymp2^{-\ell}$ denote the simulation mesh and let $\Delta_\ell$ be the coupled fine/coarse payoff correction. The main analytic estimates have the form

```math
|\mathbb E Y_{H,h_\ell}-\mathbb E Y_H|\leq B Hh_\ell,
\qquad
V_\ell:=\operatorname{Var}(\Delta_\ell)\leq V H^2h_\ell^2,
\qquad
C_\ell\leq C_{\rm sample}h_\ell^{-1}.
```

The base level has variance $O(H^2)$ and constant-scale path cost. Every sample pays for its full path prefix; a small temporal support does not make that prefix free.

For coefficient accuracy $\delta$, choose the last mesh $h_{\ell_*}\asymp\min\{1,\delta/H\}$, with constants adjusted for the bias budget. Define

```math
S=\sum_{\ell=0}^{\ell_*}\sqrt{V_\ell C_\ell}
\leq C H\sum_{\ell=0}^{\ell_*}h_\ell^{1/2}=O(H).
```

The usual MLMC allocation, for example

```math
M_\ell=\left\lceil2\delta^{-2}S
\sqrt{V_\ell/C_\ell}\right\rceil,
```

makes the sampling variance at most $\delta^2/2$. Its work is bounded by

```math
\begin{aligned}
\sum_\ell M_\ell C_\ell
&\leq2\delta^{-2}S^2+\sum_\ell C_\ell\\
&=O\!\left(C_{\rm cl}\left[1+H^2\delta^{-2}+H\delta^{-1}\right]\right).
\end{aligned}
```

The $1+H\delta^{-1}$ terms pay for at least one sample on each retained level, including the coarse-only case. Median amplification supplies the required failure probability at logarithmic overhead. This is the standard MLMC allocation [G15] applied to the same $H$-dependent bounds as the quantum construction.

### 2.2 All coefficients and Picard iterations

At temporal scale $j$ there are $O(H_j^{-1})$ details of width $H_j\asymp2^{-j}$. Consequently one curve update costs, before confidence and stability factors,

```math
\sum_{j=0}^J\frac1{H_j}
O\!\left(1+H_j^2\delta^{-2}+H_j\delta^{-1}\right)
=O\!\left(2^J+\delta^{-2}\sum_{j=0}^JH_j+(J+1)\delta^{-1}\right).
```

Only one tent per temporal scale contributes at any time, so choosing $\delta\asymp\eta/(J+1)$ yields a curve error $O(\eta)$. The geometric sum of $H_j$ is bounded, and Lipschitz interpolation uses $2^J=O(\eta^{-1})$ at fixed model parameters. The same stable Picard iteration, mean-range projection and scalar terminal stage then give

```math
W_{\rm classical}
=\widetilde O\!\left(C_{\rm cl}(T,d,K)\,\epsilon^{-2}\right).
```

Confidence, representation, arithmetic and finite-bit errors must be budgeted on both sides. This matched construction does not use an interacting particle approximation and has no additional particle-number bias. Sharing paths across several coefficients might improve constants or logarithms; it does not invalidate this upper bound.

### 2.3 What the coarse-level variance does and does not show

For a fixed coarse detail whose variance is nonzero, an IID sample average has variance $V_0/M_0$. Its root-mean-square error therefore requires $M_0=\Omega(\epsilon^{-2})$ within that estimator family. Dynkin integration and finer multilevel corrections do not remove this coarse stochastic uncertainty. This explains why the direct classical transplantation does not automatically inherit the quantum $\epsilon^{-1}$ factor.

That observation is **not** a lower bound for arbitrary classical computation on explicit circuits. A deterministic quadrature rule, an analytic calculation, a sufficiently accurate control variate, a learned conditional expectation, or a low-rank solver can avoid this particular estimator. The following sections deliberately examine those alternatives.

## 3. Established interacting-particle and iterative comparators

Szpruch–Tse [ST21, Theorem 6.2] prove, for constant diffusion and their $\mathcal M_4$ spatial/Lions-derivative assumptions, the antithetic variance bound

```math
\operatorname{Var}(\Delta_{N,h})\leq C(N^{-2}+h^2).
```

Together with their weak-bias analysis and finite-rank per-step cost $O(NC_{\rm step})$, taking $N_\ell\asymp h_\ell^{-1}\asymp2^\ell$ gives $\alpha=1$, $\beta=2$, $\gamma=2$ and hence

```math
O\!\left(C_{\rm step}\epsilon^{-2}\log^2(1/\epsilon)\right).
```

This finite-rank specialization follows from their variance theorem and MLMC allocation. It is not their displayed generic pairwise-interaction cost, which is worse. Their $\mathcal M_4$ definition excludes zeroth derivatives; bounded positive-order derivatives permit the linear harmonic drift in the polymer model. These are nevertheless stronger derivative assumptions than the core $C^2/C^3$ route above.

Bao–Reisinger–Ren–Stockinger [BRRS24] establish strong-order-one and antithetic Milstein results, including delay settings. Their sharper particle/time cost conclusion in Corollary 3.4 uses the multiplicative bound in their equation (3.8). The paper explicitly says that this bound is supported numerically but not proved in general. It must not be treated as an unconditional theorem for the entire present class.

Picard decoupling and multilevel estimation are also established classical ideas. Szpruch–Tan–Tse [STT19] develop conditionally independent iterative particle representations; Hutzenthaler–Kruse–Nguyen [HKN22] give full-history multilevel Picard approximations; Neufeld–Nguyen–Schmocker [NNS26] cover nonconstant diffusion and polynomial dimension/error cost for strong $L^2$ approximation. Their strong random-solution statements are relevant prior art, but are not identical to the small classical expectation output studied here.

The 2026 Bouillon–Samaey preprint [BS26] reports a single-ensemble $O(\epsilon^{-2-\delta})$ result under global Lipschitz assumptions for every fixed $\delta>0$. Only its primary abstract was available for this audit, so its detailed theorem is not used to establish any step of our comparison.

## 4. Can ordinary control variates match the exponent?

They can in favorable instances. Finite regularity alone does not give the required inexpensive control function.

For a fixed mean curve and deterministic initial state, let $u(t,x)=\mathbb E[f(X_T)\mid X_t=x]$. The exact backward solution gives the martingale representation

```math
f(X_T)=u(0,x_0)+\int_0^T\nabla u(t,X_t)\Sigma\,dW_t.
```

Subtracting the stochastic integral with an approximation $\widehat u$ gives an unbiased control variate, under the usual integrability and independent-training conditions. Its ideal continuous-time residual variance is

```math
\mathbb E\int_0^T
\|\Sigma^T(\nabla u-\nabla\widehat u)(t,X_t)\|^2\,dt.
```

Thus sufficiently accurate approximation of the backward gradient can eliminate most sampling work. But constructing that approximation is a high-dimensional conditional-expectation problem. Random initial data add the corresponding initial conditional-expectation term. Discretization, regression, training data and validation costs must be included; the unknown exact backward solution is not a free control variate.

This is a real classical improvement, not a merely hypothetical objection. Belomestny–Häfner–Nagapetyan–Urusov [BHNU18] derive regression-based methods with ordinary-SDE complexity exponents approaching $7/4$ under sufficiently strong approximation and tail assumptions. Their analyzed piecewise-polynomial construction has cost involving $Q^d$ and a favorable regime requiring regression degree $p>(d-2)/2$. These high-order approximation assumptions are not implied by the present fixed finite derivative bounds. The result also concerns a known ordinary SDE; solving the nonlinear mean curve is an additional cost here.

Regression on interacting particles is itself established [BS20]. Stochastic grid bundling approximates conditional expectations with local regression and can improve practical performance [CO19]. These methods should be tried when their learned representations are compact. No dimension-uniform $\epsilon^{-1}$ bound for our full class follows simply from having a regression algorithm; approximation error, sample complexity and nonlinear-solver convergence still need to be established.

### The particularly strong Rouse reference

For the polymer specialization, let $R^0$ be the interaction-free Gaussian Rouse chain driven by the same Brownian path. The identity

```math
\mathbb E[g(R_T)]
=\mathbb E[g(R_T^0)]
+\mathbb E[g(R_T)-g(R_T^0)]
```

supplies a cheap, exactly integrable reference term. At fixed horizon, synchronous stability and bounded Fourier forces give a residual of order the interaction strength, and hence a variance bound proportional to its square. This can be a substantial practical improvement. At a fixed nonzero interaction strength it changes the sampling prefactor rather than guaranteeing a new precision exponent. A hierarchy of increasingly accurate controls could do better, but its construction must be counted.

There is an important bias detail: if the reference path is an Euler-discretized Rouse chain, subtract the **discrete** Gaussian expectation. With $B_h=I-\kappa hA$, its covariance obeys $C_{n+1}=B_hC_nB_h^T+2DhI$. Subtracting a continuous-time Gaussian expectation from this discrete reference without correcting their difference introduces a time-discretization bias. A field-conditioned control must likewise condition on the same learned curve, and particle bias is not removed just by reducing variance.

## 5. Deterministic quadrature and low-dimensional PDEs

The one-dimensional Gauss-quadrature result [KS17, Theorems 5.3–5.4] is a concrete counterregime. Its favorable cost uses the paper's infinitely smooth coefficient and test-function classes, nondegenerate diffusion, exponential integrability, and an initial quadrature that can be prepared to the required accuracy at the stated cost. These assumptions yield a compressed one-dimensional quadrature rule at each step. It is not a general high-dimensional theorem under fixed-order regularity.

Cubature on Wiener space also directly targets mean-field SDE expectations. Chaudru de Raynal–Garcia Trillos [CDRGT15] and Crisan–McMurray [CM19] construct deterministic high-order methods; the latter includes Lagrange interpolation of scalar mean curves. Their smoothness and nondegeneracy/Hörmander analyses are important prior art. The branching and recombination costs of the cubature representation must be charged. A convergence order in the time step is not, by itself, a dimension-robust end-to-end work bound.

For small state dimension, a forward Fokker–Planck or backward Kolmogorov PDE solver may be superior. Sparse grids can exploit mixed smoothness; tensor methods can exploit low rank. Neither full mixed smoothness with dimension-controlled constants nor a rank bound follows from the finite number of state derivatives assumed by the quantum construction. Conversely, high dimension alone does not prove that these structures are absent. Tensor-train backward/regression solvers have demonstrated useful performance on high-dimensional parabolic equations [RSN21] and remain credible adversaries.

Affine-in-state mean-field equations admit especially powerful characteristic-function/Picard techniques [AP21]. Gaussian Rouse dynamics, exact moment closure, or an exact conditional low-dimensional reduction should be used whenever available. The proposed advantage excludes these easy structures rather than hiding them in a generic simulation baseline.

## 6. Global mean curves and temporal spectral shortcuts

Representing the unknown mean curve globally is not new. Agarwal–Amato–Pagliarani–dos Reis [AAPdR26] formulate a finite-dimensional optimization over temporal basis coefficients, use Chebyshev-point interpolation in their experiments, and establish stochastic-gradient convergence results. Their nonconvex convergence-to-stationarity statement does not supply a general epsilon-cost guarantee for locating the desired mean curve. Nevertheless, an effective deterministic or stochastic coefficient solver is a valid classical competitor.

If a particular mean curve admits a dimension-controlled analytic continuation or a certified rapidly convergent temporal spectral approximation, few Chebyshev nodes may suffice. Both classical and quantum algorithms should use that fact. Temporal analyticity may also arise from parabolic smoothing in suitable subclasses; it should not be dismissed solely because one initially states finite-order spatial assumptions.

Such a shortcut reduces the **representation and number of temporal expectation problems**. It does not automatically evaluate each high-dimensional stochastic expectation in inverse-linear classical work. For a known coefficient circuit, a deterministic solver that does both can remove the quantum advantage. Our general finite-regularity statement assumes neither this spectral approximation guarantee nor its failure on every instance. Its reconstruction argument remains valid without needing it.

There is a sharper limitation for the concrete analytic polymer example. At a **fixed** state dimension, a deterministic PDE or spectral method can eventually undercut an inverse-linear precision exponent if sufficiently strong analytic-solution and approximation bounds hold, despite an enormous dimension-dependent prefactor. Analytic coefficients alone do not certify all those bounds, but this is an active competing possibility. Consequently the example is offered as a finite high-dimensional benchmark and a restricted cost regime, not an asymptotic separation for one fixed analytic physical instance. A joint dimension/accuracy claim would require a specified growing model family, including parameter and observable scaling, and a matched analysis of its deterministic solvers.

## 7. QMC and MLQMC

Ben Rached–Haji-Ali–Tempone–Wilkosz [BRTW24] construct low-discrepancy drivers for mean-field particle systems and prove weak convergence for additive-noise SDEs. Enhanced rates are proved in a simplified deterministic mean-field ODE setting. The favorable stochastic weak-error and variance rates used in their Section 5 cost analysis are conjectures supported by numerical results. The paper is written in a scalar-state setting; it does not provide the dimension-controlled weighted mixed-derivative bounds needed to apply those rates uniformly to every high-dimensional instance here.

Their own MLQMC cost analysis gives quadratic work under its conjectures, with independent-shift sampling at coarse levels. That is not a lower bound for all QMC variants. A more effective weighted construction, coordinate ordering, higher-order scheme or coupled control variate can improve on it. When applying their analysis to our model, replace their generic pairwise interaction cost by the actual finite-rank $O(NC_{\rm step})$ cost. Quoting an unmodified pairwise-cost exponent would give the quantum algorithm an unfair comparison.

The finite-order coefficient assumptions do not guarantee favorable effective stochastic dimension or bounded high-order mixed derivatives of the path integrand. Thus no universal classical QMC $\epsilon^{-1}$ theorem has been established by the cited work for the present class. This leaves a conditional route to advantage; it does not prove that QMC cannot achieve that rate on a chosen physical model.

The $d=96$ [polymer experiment](polymer-dynamics.md#8-illustrative-classical-experiment) gives scrambled Sobol sampling both a Brownian bridge and a Rouse-mode ordering. It often improves variance by a factor of three to six, but cannot determine an asymptotic rate, certify total error or rule out stronger QMC designs. It is included to prevent high nominal dimension from being mistaken for a demonstrated obstacle to QMC.

## 8. Scope of the conditional advantage

| Comparator | Established or derived scope | Consequence for the claim |
|---|---|---|
| Same Dynkin/Picard construction with classical MLMC | $\widetilde O(\epsilon^{-2})$ from the same analytic estimates | Main matched upper-bound comparison; no additional particle bias |
| Antithetic particle/time MLMC [ST21] | $O(\epsilon^{-2}\log^2(1/\epsilon))$ after finite-rank cost specialization, under stronger derivative assumptions | Confirms that naive cubic simulation is not the baseline |
| Regression controls [BHNU18] | Subquadratic ordinary-SDE rates under stronger high-dimensional approximation assumptions | Can narrow the gain; training and MV closure are not free |
| One-dimensional Gauss quadrature [KS17] | $O(\epsilon^{-1}\lvert\log\epsilon\rvert^3)$ in the stated smooth, nondegenerate regime | No precision-exponent advantage claimed there |
| QMC, cubature, spectral, tensor or closure methods | Depend on additional regularity, effective dimension, rank or approximation accuracy | May match or beat the quantum algorithm on structured instances |
| Analytic/affine or symmetry-trivial dynamics | Exact or low-dimensional reductions can apply | Use the reduction; no generic sampling advantage claim |

The supported conclusion is an inverse-linear quantum accuracy upper bound, conditional on the specified analytic and reversible-computation hypotheses, compared with the explicit matched classical multilevel construction and established applicable methods. An unconditional separation for explicit coefficient circuits, optimality among classical algorithms, and a physical wall-clock crossover are not established.
