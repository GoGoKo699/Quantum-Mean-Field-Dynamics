# Quantum Mean-Field Dynamics

## Question and result

Can a quantum computer estimate a useful observable of a nonlinear mean-field diffusion when its feedback curve is itself unknown?

The input here is an ordinary classical description of the coefficients, a known initial state or a specified efficiently sampled initial law, and one requested observable. The task is to estimate

```math
dX_t=a(X_t,m(t))\,dt+\Sigma\,dW_t,\qquad
m(t)=\mathbb E[g(X_t)],\qquad
z=\mathbb E[f(X_T)]
```

to additive error $\epsilon$, with a specified success probability. The diffusion matrix is constant and the feedback uses a finite vector of moments. The algorithm constructs that vector's time dependence before returning one classical number.

For fixed model parameters, bounded drift with two bounded spatial derivatives, bounded mean observables with three bounded derivatives, Lipschitz moment dependence, and the stated finite-bit evaluation contract, the construction has logical gate cost $\widetilde O(\epsilon^{-1})$. A separate corollary covers a linear confining drift plus a bounded nonlinear force. The complete assumptions, explicit stability factors and proof are in [Model and results](science/model-and-results.md).

| Implementation or comparator | Precision dependence of work | Relevant memory and scope |
|---|---:|---|
| Quantum, retained path history | $\widetilde O(\epsilon^{-1})$ | Up to $\widetilde O(\epsilon^{-1})$ coherent workspace at fixed model size |
| Quantum, ternary reversible recomputation | $\widetilde O(\epsilon^{-\log_3 5})$ | Coherent workspace polynomial in logarithmic inverse precision at fixed model parameters; $\log_3 5<3/2$ |
| Classical MLMC using the same feedback construction | $\widetilde O(\epsilon^{-2})$ | Same input information, coupling and temporal reconstruction; no particle-number bias |

All three methods also store the classical feedback curve, with $\widetilde O(\epsilon^{-1})$ entries under the basic regularity assumptions. Input reading, initial sampling, derivative evaluation, classical curve reconstruction, circuit generation, numerical precision and scalar output are charged. The tilde suppresses logarithmic factors. Dimension, noise, regularity, stiffness and finite-time stability factors remain part of the cost; these are not hardware runtime estimates.

The quantum bounds are algorithmic upper bounds. Their comparison with classical sampling gives a conditional route to advantage in high-dimensional instances where inexpensive quadrature, closure, effective control variates or low-rank representations do not resolve the expectation. There is no unconditional classical lower bound for the supplied coefficient programs. In a suitable smooth one-dimensional regime, deterministic Gauss quadrature already has essentially inverse-linear cost. [Classical comparison](science/classical-comparison.md) treats that counterregime and the strongest applicable alternatives directly.

## What makes the construction work

A time detail of the mean curve, supported on an interval of width $H$, can be expressed as the expectation of one generator value at a uniformly sampled time. Coupling fine and coarse Euler paths with the same Brownian increments gives

```math
\text{detail bias}=O(Hh),\qquad
\|\text{fine/coarse correction}\|_{L^2}=O(Hh).
```

The factor $H$ permits each short temporal support to use a coarser simulation grid at the same coefficient accuracy. Quantum multilevel estimation then costs $\widetilde O(1+H/\delta)$ for a coefficient of accuracy $\delta$, including the complete path prefix. There are $O(H^{-1})$ supports at that temporal scale, so the factors cancel when constructing the whole curve. A global Picard argument controls the nonlinear feedback error without requiring a short-time contraction.

Dynkin's formula, dyadic interpolation, quantum mean estimation, quantum MLMC and reversible simulation are established ingredients. The contribution is their complete use in reconstructing unknown finite-moment feedback, with the mixed temporal/path estimate, support-dependent simulation meshes, uniform regularity of measured iterates, and explicit end-to-end resource accounting. The [related-work account](references/related-work.md) identifies the inherited ingredients and the closest quantum and classical antecedents.

## Reading routes

| Goal | Start here |
|---|---|
| Understand the computational idea | [Scientific overview](science/README.md) |
| Build the prerequisites and follow the whole argument | [Learning route](learning/README.md), [Giles bridge](learning/giles-bridge.md) and [30 milestones](learning/milestones.md) |
| Check the exact guarantee and all model constants | [Model and results](science/model-and-results.md) |
| Check the new sampling estimate | [Random-time Dynkin details](science/dynkin-details.md) |
| Check nonlinear convergence, confidence and scalar output | [Feedback reconstruction](science/feedback-reconstruction.md) |
| Audit preparation, numerical arithmetic and gates | [Finite precision and circuits](science/finite-precision-and-circuits.md) |
| Audit coherent workspace and recomputation | [Memory and recomputation](science/memory-and-recomputation.md) |
| Examine a concrete finite-dimensional task | [Polymer dynamics](science/polymer-dynamics.md) |
| Challenge the advantage claim | [Classical comparison](science/classical-comparison.md) and [related work](references/related-work.md) |
| Reproduce the numerical checks | [Verification](verification/README.md) |

The learning route uses Giles's *Multilevel Monte Carlo methods* as its single educational anchor, in the corrected author copy of the 2015 review. It starts from the classical baseline, works through a complete random-time detail example, and connects the estimator to unknown feedback and paid quantum resources. The 30 milestones comprise 12 background/preliminaries targets and 18 narrative/innovation targets.

## A concrete task and its limits

The application is a finite Gaussian bead-spring chain with a finite Fourier mean-field interaction. Its parameters, Gaussian initial sampler, initial Fourier moments, forces and requested density/scattering observables are explicit. One force-and-generator evaluation takes $O(LK)$ elementary arithmetic for $L$ beads and $K$ mean features, with the bit cost charged separately. The physical model is supplied as the input; approximating an atomistic system or a bead-continuum limit is a separate problem.

The numerical materials check the predicted $Hh$ coupling law, the full feedback reconstruction against a classically soluble example, and classical IID versus randomized-QMC sampling for a 96-coordinate chain. They do not simulate a quantum computer or establish a practical crossover. The QMC comparison includes Brownian-bridge and normal-mode ordering and finds useful classical variance reductions. It does not certify particle bias, time bias or classical optimality.

## Reproduction

With Python and the pinned dependencies in [requirements.txt](requirements.txt):

```bash
python -m pip install -r requirements.txt
python verification/check_random_time_coupling.py
python verification/check_closed_mean.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verification/classical_polymer_stress.py
python verification/check_resource_construction.py
```

The scripts write their results beside the source. [Verification](verification/README.md) records the settings, expected numerical behavior, mathematical purpose and limits of each check. Randomized experiments use fixed seeds; runtimes depend on the machine.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).
