# Thirty milestones: from multilevel sampling to self-consistent dynamics

[Home](../README.md) · [Learning route](README.md) · [Giles bridge](giles-bridge.md) · [Scientific overview](../science/README.md)

This map follows one computational question: how can an explicit classical model whose drift depends on its own unknown expectations produce one useful scalar answer, with the complete cost counted?

The allocation is **12 background/preliminaries milestones and 18 narrative/innovation milestones**. The first phase joins the objects of the problem, stochastic error, classical multilevel allocation, quantum mean estimation and finite computation into one resource framework. The second phase follows the construction from a single frozen feedback curve to a complete algorithm, then examines its memory tradeoff, physical example, classical competitors and evidence. The allocation reflects these dependencies: understanding a fast expectation estimator is only the beginning of understanding a solver that must construct its own sampling dynamics.

The single educational anchor is Michael B. Giles, *Multilevel Monte Carlo methods*, Acta Numerica 24, 259–328 (2015), [DOI 10.1017/S096249291500001X](https://doi.org/10.1017/S096249291500001X). We use the [corrected author PDF][giles], whose header is dated 2018, as identified on the [author's article page](https://people.maths.ox.ac.uk/gilesm/acta/). Its relevant route is §§1.2–1.3, §2.1 and Theorem 1, §3.1, and §§5.1–5.2. Each source pointer below names the closest connection in that route. Where the project goes beyond the review, the canonical project page supplies the additional argument; the anchor is not presented as a source for quantum estimation or nonlinear feedback.

The map gives a common understanding route, not an exhaustive list of every proof step. Each milestone has a central question, an understanding target, earlier dependencies, a role in the story, a conceptual connection and scientific pointers.

## Notation across the route

- $`X_t`$ is the stochastic state, $`m(t)=\mathbb E[g(X_t)]`$ is the vector of $`K`$ mean features, and $`z=\mathbb E[f(X_T)]`$ is the requested scalar.
- $`H`$ is the width of a temporal interpolation support; $`h_\ell`$ is a simulation time step. Their two hierarchies perform different jobs.
- $`\eta`$ is the allowed error in one feedback-map evaluation; $`\delta`$ is one scalar coefficient's tolerance; $`\epsilon`$ is the final output tolerance.
- The repository uses $`\alpha`$ for failure probability. When translating Giles's Theorem 1, write its weak-bias exponent as $`\alpha_{\mathrm{bias}}`$; $`\beta`$ describes variance decay and $`\gamma`$ describes per-sample work growth.
- A correction bound in $`L^2`$ controls its second moment. It does not give a pointwise range bound. The notation $`\widetilde O`$ suppresses the specified logarithms and retains charged model and evaluator factors.

## Phase I: background and preliminaries

### M01. The task begins with a scalar expectation

**Central question.** What must be supplied, and what answer counts as solving the computational task?

**Understanding target.** The input is a finite description of the drift, noise, observables, initial state or permitted initial sampler, and certified bounds. The output is a classical estimate of $`z`$ with additive error $`\epsilon`$ and success probability at least $`1-\alpha`$. A trajectory ensemble, full density and relative-error rare-event estimate have different output requirements.

**Dependencies and story role.** No earlier milestone. This fixes the object whose cost every later step must justify.

**Conceptual connection.** An expectation can be useful even when reconstructing the underlying probability distribution would be expensive.

**Sources.** [Giles §1.3 and §2.1][giles]; [computational task](../science/model-and-results.md#computational-task); [matched output and access](../science/classical-comparison.md#1-what-is-matched).

### M02. Additive noise and finite-time evolution

**Central question.** What does the SDE describe, and what is special about constant diffusion?

**Understanding target.** In $`dX_t=a(X_t,m(t))\,dt+\Sigma\,dW_t`$, drift gives deterministic local motion while Brownian increments give Gaussian fluctuations. Constant $`\Sigma`$ means additive noise. The horizon, dimensions and diffusion norm are supplied parameters; a fixed-horizon estimate does not give a uniform long-time guarantee.

**Dependencies and story role.** M01. The scalar target is now attached to a dynamical model and an obtainable source of randomness.

**Conceptual connection.** Discrete Gaussian increments approximate continuous stochastic evolution; their cost and error must eventually enter the same ledger as expectation estimation.

**Sources.** [Giles §5.1][giles]; [model and norms](../science/model-and-results.md#norms-and-assumptions); [known random initial laws](../science/dynkin-details.md#known-random-initial-laws).

### M03. Finite moment dependence is still a nonlinear law problem

**Central question.** Why do finitely many feedback features not automatically give a closed finite system of ODEs?

**Understanding target.** Each trajectory depends on the deterministic mean curve $`m`$, while that curve is determined by the law of the trajectories. Applying the dynamics to the selected features can create new moments outside the selected list. Finite interaction rank reduces force evaluation costs without necessarily closing the moment hierarchy.

**Dependencies and story role.** M01–M02. This identifies the obstacle beyond simulating an ordinary SDE with known coefficients.

**Conceptual connection.** A small interface between a distribution and its drift can coexist with a high-dimensional distribution.

**Sources.** [Giles §5.1][giles] for the known-coefficient SDE starting point; [self-consistent task](../science/model-and-results.md#computational-task); [explicit nonclosure example](../science/polymer-dynamics.md#4-finite-interaction-rank-does-not-close-the-state-dynamics).

### M04. The generator turns dynamics into observable change

**Central question.** How can the evolution of an expected feature be expressed through local drift and diffusion?

**Understanding target.** For a feature $`g_k`$, the generator observable is

```math
F_k=a\cdot\nabla g_k+\tfrac12\mathrm{tr}(\Sigma\Sigma^{\mathsf T}\nabla^2g_k).
```

Its magnitude and Lipschitz bounds follow from the supplied coefficient and derivative bounds. These derivatives must be computable by charged programs. The terminal observable $`f`$ has a different, weaker regularity contract.

**Dependencies and story role.** M02–M03. This provides the local quantity from which the later temporal estimator will be built.

**Conceptual connection.** Itô's formula links a differential description of a path to an integral description of its expected observables.

**Sources.** [Giles §§5.1–5.2][giles] for the SDE setting; [generator constants](../science/model-and-results.md#constants-used-in-the-proof); [generator identity](../science/dynkin-details.md#the-generator-and-a-temporal-detail).

### M05. Bias, strong error, sampling error and confidence

**Central question.** Which errors concern paths, expectations, random estimates and successful output?

**Understanding target.** Strong error compares coupled random states; bias compares expectations; sampling error concerns an estimator's fluctuation; a failure allowance bounds the probability that a specified tolerance is missed. A small mean-square error and a prescribed high-probability guarantee are related but distinct specifications. Smoothness, norms and constants determine which bounds apply.

**Dependencies and story role.** M01–M04. This establishes the distinctions needed to combine numerical and statistical estimates without confusing their units or guarantees.

**Conceptual connection.** A path estimate can imply a payoff bias and a correction second-moment bound, which then enter different parts of the work calculation.

**Sources.** [Giles §2.1, Theorem 1, and §§5.1–5.2][giles]; [explicit accuracy parameters](../science/model-and-results.md#explicit-accuracy-parameters); [Euler error statement](../science/dynkin-details.md#additive-noise-euler-error-with-an-explicit-constant).

### M06. Coupling makes differences easier to estimate

**Central question.** Why should two approximations share their random inputs?

**Understanding target.** Fine and coarse paths use the same initial draw and compatible increments of one Brownian path. This can make their payoff difference much smaller than either payoff. Independence between the two paths would usually lose that benefit. Coupling must preserve each path's intended marginal law as well as produce correlation.

**Dependencies and story role.** M02 and M05. The route moves from the error of one simulated path to the value of comparing two simulations.

**Conceptual connection.** Shared randomness is the mechanism that turns a numerical convergence estimate into a useful multilevel correction.

**Sources.** [Giles §1.2 and §5.1][giles]; [coupled random-time estimators](../science/dynkin-details.md#coupled-random-time-estimators); [shared partial increments](../science/dynkin-details.md#a-coupled-partial-step).

### M07. Telescoping reorganizes one expectation

**Central question.** How can a sequence of increasingly accurate approximations reduce total sampling work?

**Understanding target.** A finest-level expectation is the coarse expectation plus the sum of adjacent-level expected differences. Each difference needs a coupled pair, while different levels can be estimated separately. The coarse expectation remains part of the computation even when all fine corrections are small.

**Dependencies and story role.** M05–M06. This turns the local coupling idea into a complete estimator hierarchy.

**Conceptual connection.** Telescoping changes how work is distributed without changing the ideal finest-level expectation.

**Sources.** [Giles §§1.2–1.3][giles]; [path-level hierarchy](../science/dynkin-details.md#quantum-multilevel-cost-for-one-coefficient); [classical detail hierarchy](../science/classical-comparison.md#21-one-temporal-detail).

### M08. Classical allocation depends on variance and cost together

**Central question.** Why is the finest discretization error insufficient to determine computational complexity?

**Understanding target.** Classical MLMC balances level variance $`V_\ell`$ against per-sample cost $`C_\ell`$. Its allocation is proportional to $`\sqrt{V_\ell/C_\ell}`$, with a common factor set by the target variance. Rounding sample counts retains a mandatory sample-cost term. Small corrections can still be expensive if evaluating one correction requires a long path. The project's supplied certified constants differ from empirical variance and bias diagnostics.

**Dependencies and story role.** M05–M07. This supplies the classical work calculation that later quantum and classical constructions must both face.

**Conceptual connection.** The important quantity is the cost of obtaining the requested accuracy, including every retained level.

**Sources.** [Giles §1.3, §2.1/Theorem 1 and §3.1][giles]; [matched classical allocation](../science/classical-comparison.md#21-one-temporal-detail).

### M09. Structure can reduce the sampling problem itself

**Central question.** What can a strong classical method exploit before increasing its sample count?

**Understanding target.** A known Gaussian reference, symmetry, control variate or reduced representation can remove part of the uncertainty. Constructing and evaluating a control must be charged. Finite interaction rank permits efficient aggregation on both sides of the comparison. Nominal dimension alone does not establish the absence of a useful reduction.

**Dependencies and story role.** M03, M06 and M08. This prevents the classical baseline from being fixed prematurely to a naive simulator.

**Conceptual connection.** Variance reduction and model structure change the estimator that quantum acceleration is supposed to improve.

**Sources.** [Giles §1.2 and §2.1][giles]; [control variates](../science/classical-comparison.md#4-can-ordinary-control-variates-match-the-exponent); [Rouse reference](../science/classical-comparison.md#the-particularly-strong-rouse-reference).

### M10. Quantum estimation needs coherent sample access

**Central question.** What resource replaces repeated independent averaging in the quantum construction?

**Understanding target.** A deterministic finite sampler acts on uniform random bits and is implemented reversibly, together with its inverse. A supplied second-moment bound $`\mathbb E[Y^2]\leq\sigma^2`$ supports mean-estimation query dependence near $`1+\sigma/e`$ for additive tolerance $`e`$. Binning handles signed outputs and tails; a standard deviation cannot be used as a pointwise normalization bound.

**Dependencies and story role.** M05–M08. This introduces the quantum primitive while keeping its access requirements visible.

**Conceptual connection.** The same small coupled correction can support different classical and quantum estimation costs.

**Sources.** [Giles §1.3 and §2.1][giles] for the classical comparison; [second-moment quantum interface](../science/finite-precision-and-circuits.md#6-an-explicit-second-moment-quantum-mean-estimator).

### M11. A query must become finite computation

**Central question.** What must be counted when an abstract sample becomes a quantum circuit?

**Understanding target.** Gaussian generation, coefficient evaluation, arithmetic, controls, retained history, uncomputation and synthesis all have costs. Smooth functions do not automatically have inexpensive evaluators. Classical circuit-generation work and coherent workspace are separate resources. An exact finite Boolean sample circuit can still require approximate non-Boolean gates in the estimation procedure.

**Dependencies and story role.** M01, M02 and M10. This connects the estimation primitive to the actual programs supplied by the input model.

**Conceptual connection.** Reversibility provides coherent access, while finite precision and memory determine whether that access preserves the promised advantage.

**Sources.** [Giles §2.1 and §3.1][giles] for charging sample work; [bit-model contract](../science/finite-precision-and-circuits.md#1-bit-model-and-the-input-contract); [reversible gates and synthesis](../science/finite-precision-and-circuits.md#7-reversible-gates-synthesis-and-success-probability).

### M12. The complete comparison has one input and one output

**Central question.** What would constitute a justified end-to-end gain for this task?

**Understanding target.** Both methods receive the same explicit programs and may exploit their structure. The comparison includes preparation, construction of the unknown feedback, numerical and statistical accuracy, compilation, memory and scalar readout. A quantum upper bound compared with a specified classical upper bound is not a lower bound against every classical algorithm.

**Dependencies and story role.** M01–M11. This completes the background phase by joining the dynamical model, error budget, estimator and access model.

**Conceptual connection.** The remaining question is whether a self-consistent sampler can be built cheaply enough for quantum mean estimation to matter.

**Sources.** [Giles §2.1/Theorem 1 and §3.1][giles]; [charged model factors](../science/model-and-results.md#a-resource-bound-with-the-model-factors-exposed); [matched access](../science/classical-comparison.md#1-what-is-matched).

## Phase II: narrative and innovation

### M13. Freeze the feedback to obtain an ordinary SDE

**Central question.** How can an unknown nonlinear mean curve be approached using known-coefficient simulations?

**Understanding target.** A supplied curve $`v`$ defines $`X^v`$ and the map $`\Phi(v)(t)=\mathbb E[g(X_t^v)]`$. The desired mean is the fixed point $`m=\Phi(m)`$. Each call uses a known classical curve; uniform bounds over admissible curves make later measured feedback legitimate. The generator bounds make $`\Phi(v)`$ uniformly Lipschitz.

**Dependencies and story role.** M03–M04, M06 and M12. This is the first step from the task to a computable construction.

**Conceptual connection.** Nonlinearity is moved into a classical fixed-point iteration around ordinary stochastic simulations.

**Sources.** [Giles §5.1][giles] for the known-coefficient starting point; [frozen-feedback map](../science/feedback-reconstruction.md#the-frozen-feedback-map).

### M14. A temporal detail is one random-time generator expectation

**Central question.** How can a midpoint interpolation correction be estimated without separately learning three means?

**Understanding target.** Subtracting the generator-integral formulas cancels the common past. For support $`I`$ of width $`H`$, its midpoint detail equals the mean of

```math
Y_{k,I}=\frac H2\,s_I(\tau)F_k(X_\tau^v,v(\tau)),
\qquad \tau\sim\mathrm{Unif}(I).
```

The sign changes at the midpoint. An endpoint coefficient has the corresponding unsigned integral representation.

**Dependencies and story role.** M04 and M13. This isolates the estimator that carries the contribution.

**Conceptual connection.** A deterministic temporal difference becomes a single stochastic integrand with an explicit support-width factor.

**Sources.** [Giles §§1.2–1.3][giles] for estimating corrections; [random-time representation, Lemma 1](../science/dynkin-details.md#the-generator-and-a-temporal-detail).

### M15. Preserve the mixed support and mesh factor

**Central question.** Why do both the bias and coupled correction scale with $`Hh_\ell`$?

**Understanding target.** Additive-noise Euler has fixed-time $`L^2`$ strong order one under the stated regularity of the frozen curve. Sharing the random terminal time and Brownian increments, including the partial final increment, transfers that estimate through the Lipschitz generator. Multiplication by $`H`$ gives correction $`L^2`$ size $`O(Hh_\ell)`$ and variance at most $`O(H^2h_\ell^2)`$.

**Dependencies and story role.** M05–M06 and M13–M14. This proves the estimate on which the subsequent cost cancellation depends.

**Conceptual connection.** Temporal localization reduces the payoff scale even though every sample still simulates its entire prefix from zero.

**Sources.** [Giles §§5.1–5.2][giles]; [Euler estimate, Lemma 2](../science/dynkin-details.md#additive-noise-euler-error-with-an-explicit-constant); [mixed estimate, Lemma 3](../science/dynkin-details.md#coupled-random-time-estimators).

### M16. Choose a different path depth for each support

**Central question.** How does the mixed estimate change the work needed for one coefficient?

**Understanding target.** A coefficient tolerance $`\delta`$ permits a final mesh of order $`\delta/H`$ when refinement is needed. Combining the $`Hh_\ell`$ second-moment scale with a full-path sample cost of order $`h_\ell^{-1}`$ gives near $`1+H/\delta`$ quantum work per coefficient. The one-sample term and the level-zero case remain included.

**Dependencies and story role.** M07–M08, M10 and M15. This converts the local estimate into its first computational consequence.

**Conceptual connection.** A short support permits a coarser simulation; imposing one global finest mesh would discard that saving.

**Sources.** [Giles §2.1/Theorem 1 and §3.1][giles] for bias/cost allocation; [one-coefficient quantum cost](../science/dynkin-details.md#quantum-multilevel-cost-for-one-coefficient).

### M17. Learn a curve through dyadic details

**Central question.** How can many noisy coefficients represent the whole mean curve with controlled uniform error?

**Understanding target.** A piecewise-linear dyadic interpolant is an initial value, an endpoint increment and a hierarchy of midpoint hats. At any time, only one hat per temporal scale is nonzero. Consequently coefficient errors accumulate with the number of scales, rather than the total number of coefficients. Lipschitz interpolation needs order $`\eta^{-1}`$ stored nodes at fixed model parameters.

**Dependencies and story role.** M14 and M16. The route moves from one coefficient to one approximate feedback-map evaluation.

**Conceptual connection.** Local support makes a large classical representation compatible with a small uniform error.

**Sources.** [Giles §1.3][giles] for hierarchical correction bookkeeping; [dyadic mean-curve representation](../science/feedback-reconstruction.md#a-dyadic-representation-of-the-whole-mean-curve).

### M18. Preserve the regularity needed by the next iteration

**Central question.** Why is closeness in value insufficient for a usable reconstructed feedback curve?

**Understanding target.** Euler's error estimate also needs a uniform time-Lipschitz bound. Hat derivatives amplify coefficient errors by inverse support width. The chosen coefficient tolerance and temporal mesh bound that amplification. Projection of grid nodes into the permitted mean box preserves both the relevant error bound and a uniform slope bound.

**Dependencies and story role.** M04 and M17. This closes the invariant that allows the estimator to be reused on its own output.

**Conceptual connection.** Numerical reconstruction must preserve the analytic assumptions of the following computation, not only approximate the present values.

**Sources.** [Giles §2.1][giles] for checking estimator assumptions; [reconstruction and slope bound, Lemma 5](../science/feedback-reconstruction.md#reconstruction-error-and-a-uniform-slope-bound).

### M19. Global Picard iteration has factorial stability

**Central question.** How can inexact feedback iteration converge without assuming a short-time contraction?

**Understanding target.** Synchronous state coupling gives a Volterra integral bound for differences between frozen-feedback maps. Repeated integration yields a factorially decaying initial-error term. The accumulated map error is bounded by $`\eta e^{(L_x+L_gL_m)T}`$. Fixed finite time permits logarithmically many sufficient iterations, while the displayed stability factors remain charged.

**Dependencies and story role.** M13 and M18. This upgrades a reliable single map evaluation to a solver for the unknown mean curve.

**Conceptual connection.** The causal time-integral structure supplies more information than a single global norm-contraction factor.

**Sources.** [Giles §2.1][giles] for separating approximation errors; [Volterra and inexact Picard stability, Lemma 4](../science/feedback-reconstruction.md#stability-without-a-short-time-contraction-assumption).

### M20. Adaptive classical feedback retains a global confidence guarantee

**Central question.** Why do earlier random estimates not invalidate the next estimator's assumptions?

**Understanding target.** Conditional on previous measured outputs, the next curve is a fixed classical input. Range projection and the slope check ensure that every invoked sampler is admissible. Uniform conditional bounds support a union bound without independence between calls. Coefficients and terminal estimation each receive total ideal failure allowance $`\alpha/4`$; synthesis receives a separate $`\alpha/2`$ output-law allowance.

**Dependencies and story role.** M05, M10 and M18–M19. This attaches the promised success probability to the full nonlinear construction.

**Conceptual connection.** Classical adaptivity avoids coherent nesting of expectation estimators while retaining an auditable probability budget.

**Sources.** [Giles §2.1 and §3.1][giles] for accuracy allocation; [complete algorithm](../science/feedback-reconstruction.md#the-complete-algorithm); [conditional confidence argument](../science/feedback-reconstruction.md#why-adaptivity-does-not-multiply-failure-probabilities).

### M21. Read out the requested observable

**Central question.** Once the mean curve is accurate, what remains before the scalar answer is available?

**Understanding target.** Stability transfers curve error into terminal expectation error through the explicit factor $`D_f`$. A final coupled hierarchy estimates $`\mathbb E[f(X_T^{v_R})]`$ and consumes the remaining output-error budget. A bounded Lipschitz $`f`$ suffices: no derivative or generator evaluator for $`f`$ is required. A retained mean coordinate can instead be read from its reconstructed endpoint.

**Dependencies and story role.** M05–M06, M15 and M19–M20. This closes the output contract established in M01.

**Conceptual connection.** Learning the feedback and evaluating the final observable are separate computations with linked tolerances.

**Sources.** [Giles §§5.1–5.2][giles]; [terminal readout and cost](../science/feedback-reconstruction.md#terminal-readout-and-its-resource-cost); [density and scattering outputs](../science/polymer-dynamics.md#3-useful-output).

### M22. Sum the entire direct algorithm

**Central question.** Why does reconstructing all coefficients not restore an extra inverse power of accuracy?

**Understanding target.** At scale $`j`$, the $`2^j`$ supports have width $`H_j=T2^{-j}`$, so their total width is $`T`$. Summing $`1+H/\delta`$ includes both support-weighted sampling work and the unavoidable number of coefficients. The scale count, Picard count and confidence overhead add logarithms. Terminal readout has the same direct near $`\epsilon^{-1}`$ precision exponent.

**Dependencies and story role.** M16–M17 and M19–M21. This integrates the analytic construction into its complete direct work bound.

**Conceptual connection.** The support-width saving is useful only because it survives summation over the entire unknown curve.

**Sources.** [Giles §1.3 and §2.1/Theorem 1][giles]; [whole-curve cost sum](../science/feedback-reconstruction.md#summing-the-cost-of-a-whole-feedback-curve); [exposed resource factors](../science/model-and-results.md#a-resource-bound-with-the-model-factors-exposed).

### M23. Finite-bit paths need not telescope exactly

**Central question.** What replaces exact telescoping after endpoint, Gaussian and arithmetic discretization?

**Understanding target.** Quantized fine increments need not sum to the separately quantized coarse law. Each finite correction is instead coupled to its ideal correction with a summable $`L^2`$ error budget. Aligned endpoint bins preserve sign and update masks. Explicit Gaussian tails, finite ranges, update errors and output rounding yield logarithmic word lengths for fixed admissible models.

**Dependencies and story role.** M11, M15–M16 and M20–M22. This turns the real-arithmetic hierarchy into the finite-program object used by the theorem.

**Conceptual connection.** Approximate agreement with a common ideal hierarchy is sufficient; false exact identities between discretized distributions are unnecessary.

**Sources.** [Giles §2.1 and §5.1][giles] for the ideal hierarchy; [aligned endpoints and stability](../science/finite-precision-and-circuits.md#4-aligned-endpoint-bits-and-numerical-stability); [nonmatching marginals and budgets](../science/finite-precision-and-circuits.md#5-nonmatching-discrete-marginals-and-the-full-coefficient-error-budget).

### M24. Compress randomness while preserving both moments

**Central question.** What must pseudorandom replacement preserve for quantum multilevel estimation to remain inexpensive?

**Understanding target.** The finite classical sampler is a bounded-space, fixed-clock computation. HashPRG approximates its final output distribution. An all-seed payoff bound converts that distributional error into separate bounds on the mean and second moment. A mean-only guarantee would not preserve the cost estimate. The quantum algorithm evaluates the resulting deterministic function over uniform short seeds.

**Dependencies and story role.** M10–M11 and M23. This addresses the long random tape in the direct coherent implementation.

**Conceptual connection.** Classical distributional pseudorandomness can support a quantum sampler without a claim that the generator fools quantum distinguishers.

**Sources.** [Giles §2.1/Theorem 1][giles] for the separate roles of bias and variance; [explicit HashPRG](../science/memory-and-recomputation.md#2-an-explicit-block-addressable-classical-prg); [both-moment transfer](../science/memory-and-recomputation.md#23-first-and-second-moments-must-both-be-preserved).

### M25. Recompute paths to reduce coherent memory

**Central question.** How can an irreversible sample program be run coherently without retaining its entire history?

**Understanding target.** A clean reversible bridge computes a step into a target register and uncomputes temporary data, even when the original step is not invertible. Ternary checkpointing implements $`3^r`$ steps with $`5^r`$ leaf calls and logarithmically many live checkpoints. The resulting compact gate exponent is $`\log_3 5`$, approximately $`1.465`$, with fixed-model polylogarithmic coherent workspace.

**Dependencies and story role.** M11 and M24. This gives a concrete time-memory tradeoff rather than treating reversible workspace as free.

**Conceptual connection.** Regenerating the same pseudorandom blocks makes deterministic replay possible; extra work replaces stored coherent history.

**Sources.** [Giles §2.1][giles] for re-evaluating sample cost; [clean bridge and ternary schedule](../science/memory-and-recomputation.md#3-explicit-reversible-checkpoint-circuits); [compact multilevel sum](../science/memory-and-recomputation.md#5-the-multilevel-and-temporal-cost-sums).

### M26. Compile feedback access within the same budget

**Central question.** Does a small quantum sampler hide expensive access to a large classical curve?

**Understanding target.** Deterministic path-grid feedback values can be precomputed by sorted merges and embedded as constants. Processing path levels in the outer loop permits sequential coefficient and table scans. Checkpoint tape travel is of the same order as the leaf count. This charges classical compilation and avoids quantum address lookup. The classical curve still occupies near inverse-linear storage.

**Dependencies and story role.** M17, M22–M23 and M25. This closes the access-model and compilation parts of the compact implementation.

**Conceptual connection.** A fixed circuit schedule can use large classical data without placing that data in coherent random-access memory.

**Sources.** [Giles §3.1][giles] for implementation-aware allocation; [uniform compilation and storage](../science/memory-and-recomputation.md#4-uniform-compilation-and-classical-feedback-storage); [complete compact theorem](../science/memory-and-recomputation.md#6-confidence-quantum-gates-and-the-complete-resource-theorem).

### M27. Instantiate the task with a finite polymer chain

**Central question.** What makes the abstract input and output concrete without granting free physical preparation?

**Understanding target.** A specified finite Rouse chain with a finite Fourier mean field has an explicit Gaussian coil sampler, computable initial features and $`O(LK)`$ force/generator arithmetic for $`L`$ beads. Density relaxation and intrachain scattering are distinct useful outputs. The linear harmonic drift requires the separate fourth-moment extension. Finite rank, finite bead count, stiffness and model approximation remain explicit.

**Dependencies and story role.** M02–M04, M21 and M23–M26. This connects the algorithmic construction to a physically interpretable finite input model.

**Conceptual connection.** An application supports useful computation only when its preparation, observables and regularity satisfy the same contract as the theorem.

**Sources.** [Giles §§5.1–5.2][giles]; [polymer input and outputs](../science/polymer-dynamics.md); [linear-confining estimates](../science/dynkin-details.md#linear-confining-drift); [concrete finite arithmetic](../science/finite-precision-and-circuits.md#10-concrete-arithmetic-for-the-finite-bead-fourier-model).

### M28. Challenge the result with matched classical methods

**Central question.** Which comparison is established, and which classical structures can eliminate the gain?

**Understanding target.** The identical Dynkin/Picard construction with ordinary MLMC has a near $`\epsilon^{-2}`$ upper bound under the same analytic assumptions. Finite-rank antithetic methods, Rouse controls, regression, QMC, quadrature, PDEs, closures and tensor methods are additional competitors. One-dimensional smooth quadrature and fixed-dimensional analytic spectral regimes can match or beat the precision exponent when their additional approximation assumptions hold. Coarse-level IID variance is not a general classical lower bound.

**Dependencies and story role.** M08–M09, M12, M22 and M26–M27. This gives the algorithmic bound its justified comparative meaning.

**Conceptual connection.** A useful conditional advantage requires both a fast quantum construction and a credible accounting of the classical alternatives available from the same input.

**Sources.** [Giles §1.2 and §2.1/Theorem 1][giles]; [same-assumption classical derivation](../science/classical-comparison.md#2-a-direct-classical-analogue-under-the-same-analytic-assumptions); [conditional comparison and counterregimes](../science/classical-comparison.md#8-scope-of-the-conditional-advantage).

### M29. Interpret validation at its actual scope

**Central question.** What do the numerical and exact finite checks establish?

**Understanding target.** Coupling slopes probe the $`Hh`$ mechanism; the classically soluble example probes full feedback, slope control and terminal bookkeeping; the polymer experiment measures fixed-mesh randomization variance; exact finite resource checks probe indexing, masks, budgets, hashes and cleanup. None alone proves the general theorem, certifies all model promises, executes the quantum algorithm or establishes a hardware crossover.

**Dependencies and story role.** M15, M21, M23–M28. This connects the theoretical account to reproducible evidence without replacing analysis by fitted rates.

**Conceptual connection.** Different checks answer different questions, and their limitations are part of interpreting a computational result.

**Sources.** [Giles §3.1 and §§5.1–5.2][giles]; [verification scope and recorded checks](../verification/README.md); [classical experiment interpretation](../science/polymer-dynamics.md#8-illustrative-classical-experiment).

### M30. State the complete conditional contribution

**Central question.** What survives when the construction, resource accounting, classical challenge and evidence are considered together?

**Understanding target.** The contribution combines support-scaled random-time details, coefficient-dependent path depth and stable reconstruction of unknown feedback, with finite-bit and resource closure. At fixed admissible models it yields direct near $`\epsilon^{-1}`$ work, or compact near $`\epsilon^{-\log_3 5}`$ work with polylogarithmic coherent workspace. The claim retains evaluator, dimension, stability and classical-memory costs, and does not establish unrestricted classical separation or physical crossover.

**Dependencies and story role.** M12 and M22–M29. This completes the narrative by joining the useful scalar task, mechanism, implementation choices and limits.

**Conceptual connection.** Established multilevel, estimation and reversible-computation tools become a consequential construction only when their interaction solves the unknown-feedback problem at the claimed total cost.

**Sources.** [Giles §2.1/Theorem 1][giles] for the classical complexity framework; [main theorem and scope](../science/model-and-results.md#main-theorem); [compact resource theorem](../science/memory-and-recomputation.md#6-confidence-quantum-gates-and-the-complete-resource-theorem); [inherited and derived components](../references/related-work.md#what-is-inherited-and-what-is-derived-here).

[giles]: https://people.maths.ox.ac.uk/gilesm/files/acta15.pdf
