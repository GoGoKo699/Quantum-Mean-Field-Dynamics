# Learning quantum mean-field dynamics

[Home](../README.md) · [Scientific overview](../science/README.md) · [Giles bridge](giles-bridge.md) · [Thirty milestones](milestones.md)

The computational task starts with an explicit stochastic model and ends with one classical observable. The difficult step is that the model's drift depends on a curve of expectations that has not yet been computed. A useful quantum estimator must pay for constructing that curve, simulating its trajectories and returning the requested answer.

This route develops that complete calculation. It begins with a strong classical multilevel method, identifies the support-width factor in a new stochastic estimator, and follows the factor through feedback reconstruction, finite computation, memory and the classical comparison.

## One source, with a focused reading route

The educational anchor is **Michael B. Giles, Multilevel Monte Carlo methods**, *Acta Numerica* **24**, 259–328 (2015), [DOI 10.1017/S096249291500001X](https://doi.org/10.1017/S096249291500001X). Use the [corrected 70-page author PDF](https://people.maths.ox.ac.uk/gilesm/files/acta15.pdf), linked from the [author's article page](https://people.maths.ox.ac.uk/gilesm/acta/). Its header says 2018; the journal citation remains 2015.

The relevant sections are **1.2–1.3** for coupling and telescoping, **2.1 and Theorem 1** for bias/variance/cost, **3.1** for allocation and refinement, and **5.1–5.2** for SDE coupling and numerical convergence. The [source bridge](giles-bridge.md) translates the notation, explains what those sections supply and develops the additional project-specific steps. Reading the entire review is not a prerequisite for following the bridge.

Quantum mean estimation, pseudorandom generators and reversible computation have their own primary research citations in the scientific account. Those citations establish proof dependencies and attribution; Giles remains the single educational anchor.

## The connected route

The [30-milestone map](milestones.md) has two phases. **M01–M12** join the scalar task, stochastic model, generator, coupling, classical allocation, quantum access and finite computation into a common comparison framework. **M13–M30** follow the estimator from one frozen feedback curve to the complete nonlinear algorithm, then examine the physical example, classical threats and evidence. Every milestone states its understanding target, earlier dependencies, conceptual role and source pointers.

| Part of the story | Milestones | Main scientific destination |
|---|---|---|
| Define the task, dynamics and generator | M01–M04 | [Model and results](../science/model-and-results.md) |
| Understand error, coupling and classical allocation | M05–M09 | [Classical comparison](../science/classical-comparison.md) |
| Connect quantum access to finite computation | M10–M12 | [Finite precision and circuits](../science/finite-precision-and-circuits.md) |
| Build the random-time detail and retain its support width | M13–M16 | [Dynkin details](../science/dynkin-details.md) |
| Reconstruct the mean, preserve regularity and close the algorithm | M17–M22 | [Feedback reconstruction](../science/feedback-reconstruction.md) |
| Pay for finite bits, randomness, coherent memory and compilation | M23–M26 | [Memory and recomputation](../science/memory-and-recomputation.md) |
| Evaluate the application, competitors and supported conclusion | M27–M30 | [Polymer dynamics](../science/polymer-dynamics.md), [classical comparison](../science/classical-comparison.md) and [verification](../verification/README.md) |

## A small complete example

The bridge computes a midpoint detail of the curve $\phi(t)=\mathbb E[\cos W_t]=e^{-t/2}$. One signed generator value at a uniformly sampled time has exactly that detail as its expectation. The example makes the sign, support multiplier, random inputs and scalar output explicit. Its answer is analytically known, so it illustrates the mechanism without serving as an advantage example. [Work through the calculation](giles-bridge.md#5-a-complete-random-time-detail-example).

The nonlinear problem then adds the work absent from this toy calculation: a full simulated path prefix, an unknown mean curve, uniform control over measured iterates, finite-bit sampling and confidence across all calls.

## What the final synthesis must preserve

At fixed admissible model parameters, the direct-history quantum construction has near inverse-linear logical work in the output accuracy. The separate compact coherent-memory construction has exponent $\log_3 5$. Classical curve storage and input-program costs remain charged in both cases.

The same-assumption classical multilevel construction has a near quadratic upper bound. Stronger classical structure can improve that comparison, including a smooth one-dimensional quadrature regime with inverse-linear work. A conditional algorithmic comparison does not establish a universal classical lower bound, a hard explicit physical family or a hardware crossover. These qualifications belong to the result itself and remain part of the final milestones.

Readers checking the proof can go directly to [the exact theorem and resource ledger](../science/model-and-results.md). Readers checking empirical support can go to [the reproducible checks and their limits](../verification/README.md).
