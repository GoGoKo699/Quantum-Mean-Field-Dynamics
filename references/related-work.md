# Related work and source audit

This bibliography records the quantum, classical and physical sources used in the scientific account. The audit date is **11 November 2026**. A published citation is distinguished from the manuscript version actually inspected; a result reported only in an abstract is not used as a theorem premise. The structured counterpart is [source-audit.json](source-audit.json).

## What is inherited and what is derived here

Dynkin's identity, Picard iteration, multilevel Monte Carlo, mean-field particle approximation, temporal interpolation and quantum mean estimation are established tools. Classical work already reconstructs mean curves globally and combines Picard iteration with multilevel sampling. Quantum parametric integration already treats function-valued outputs by multiscale interpolation.

The proposed contribution is the specific feedback construction described in the scientific notes: a randomized generator representation of each temporal detail, retention of its support-width factor through the coupled simulation bounds, coefficient-specific simulation depths, and a stable reconstruction of the unknown mean curve. Merely applying amplitude estimation to a particle simulator, or merely using a temporal basis, is not presented as new.

The same-assumption classical $\widetilde O(\epsilon^{-2})$ result in [Classical comparison](../science/classical-comparison.md#2-a-direct-classical-analogue-under-the-same-analytic-assumptions) is a derivation in this repository using ordinary MLMC allocation. It is an upper bound, not a citation to an optimality theorem. The finite-rank specialization of the Szpruch–Tse particle/time result is also identified as a cost inference. Neither establishes a lower bound for all classical algorithms with access to explicit coefficient programs.

## 1. Multilevel and interacting-particle methods

### G15: Multilevel Monte Carlo allocation

Michael B. Giles. **Multilevel Monte Carlo methods.** *Acta Numerica* **24**, 259–328 (2015). [DOI: 10.1017/S096249291500001X](https://doi.org/10.1017/S096249291500001X). [Author's corrected full text](https://people.maths.ox.ac.uk/gilesm/files/acta15.pdf), linked with its correction notice from the [article page](https://people.maths.ox.ac.uk/gilesm/acta/). The inspected 70-page corrected PDF has a 2018 header; the journal article remains the 2015 publication.

**Checked:** Sections 1.2–1.3 and Theorem 1 in Section 2.1. Used for the standard variance/cost sample allocation. Our comparison displays the allocation and sums the actual coefficient costs, rather than treating a time-step convergence rate as a work bound.

### ST21: Antithetic particle and time discretization

Łukasz Szpruch and Alvin Tse. **Antithetic multilevel sampling method for nonlinear functionals of measure.** *Annals of Applied Probability* **31**(3), 1100–1139 (2021). [DOI: 10.1214/20-AAP1614](https://doi.org/10.1214/20-AAP1614). [Inspected author manuscript, arXiv:1903.07063v3](https://www.pure.ed.ac.uk/ws/portalfiles/portal/227385780/1903.07063.pdf).

**Checked:** Theorem 6.2, weak-bias analysis in Appendix B, and Definition A.3. Constant diffusion and $\mathcal M_4$ regularity give the stated antithetic variance estimate. Zeroth derivatives are excluded from the definition; a linear harmonic drift is therefore permitted. These conditions are stronger than the main finite-regularity theorem here. The finite-rank work specialization is our inference.

### BRRS24: Milstein and antithetic methods, including delay

Jianhai Bao, Christoph Reisinger, Panpan Ren and Wolfgang Stockinger. **Milstein schemes and antithetic multilevel Monte Carlo sampling for delay McKean–Vlasov equations and interacting particle systems.** *IMA Journal of Numerical Analysis* **44**(4), 2437–2479 (2024; first online 2023). [DOI: 10.1093/imanum/drad064](https://doi.org/10.1093/imanum/drad064). [Inspected manuscript, arXiv:2005.01165v3](https://arxiv.org/pdf/2005.01165).

**Checked:** Section 3, equation (3.8), and Corollary 3.4. The sharper particle/time complexity conclusion depends on a multiplicative variance assumption that the paper does not prove generally. Its proven temporal/strong estimates and its conditional cost result must be distinguished.

### STT19: Iterative multilevel particle approximation

Łukasz Szpruch, Shuren Tan and Alvin Tse. **Iterative multilevel particle approximation for McKean–Vlasov SDEs.** *Annals of Applied Probability* **29**(4), 2230–2265 (2019). [DOI: 10.1214/18-AAP1452](https://doi.org/10.1214/18-AAP1452). [Author manuscript, arXiv:1706.00907](https://www.pure.ed.ac.uk/ws/files/179707667/1706.00907.pdf).

Establishes prior art for Picard decoupling, conditional independence and multilevel mean-field approximation. Its generic interaction-cost theorem is not quoted unchanged as the baseline for a finite-rank model.

### HKN22: Full-history multilevel Picard

Martin Hutzenthaler, Thomas Kruse and Tuan Anh Nguyen. **Multilevel Picard approximations for McKean–Vlasov stochastic differential equations.** *Journal of Mathematical Analysis and Applications* **507**(1), 125761 (2022). [DOI: 10.1016/j.jmaa.2021.125761](https://doi.org/10.1016/j.jmaa.2021.125761). [arXiv:2103.10870](https://arxiv.org/abs/2103.10870).

Relevant for polynomial-dimension multilevel Picard approximation. A strong approximation of a random solution is a different output from the scalar classical expectation considered here.

### NNS26: Multilevel Picard with nonconstant diffusion

Ariel Neufeld, Tuan Anh Nguyen and Philipp Schmocker. **Multilevel Picard approximations for McKean–Vlasov stochastic differential equations with nonconstant diffusion.** *IMA Journal of Numerical Analysis*, online first, article drag078 (2026). [DOI: 10.1093/imanum/drag078](https://doi.org/10.1093/imanum/drag078). [Inspected manuscript, arXiv:2502.03205v2](https://arxiv.org/html/2502.03205v2).

The global-Lipschitz strong-$L^2$ analysis and polynomial dimension/error dependence are important comparators. They do not directly certify an identical work exponent for our requested expectation output.

### BS26: Recent single-ensemble result, abstract-only audit

Arne Bouillon and Giovanni Samaey. **Convergence of a single-ensemble multilevel scheme for McKean–Vlasov SDEs.** Preprint, [arXiv:2609.10000](https://arxiv.org/abs/2609.10000) (2026).

The primary abstract reports $O(\epsilon^{-2-\delta})$ cost for fixed $\delta>0$ under global Lipschitz assumptions. **Only the abstract was available in this audit.** No detailed theorem from this preprint is used to prove the comparison.

## 2. Strong classical alternatives

### KS17: A deterministic inverse-linear counterregime

Peter Kloeden and Tony Shardlow. **Gauss-quadrature method for one-dimensional mean-field SDEs.** *SIAM Journal on Scientific Computing* **39**(6), A2784–A2807 (2017). [DOI: 10.1137/16M1095688](https://doi.org/10.1137/16M1095688). [Inspected manuscript, arXiv:1608.06741v3](https://arxiv.org/pdf/1608.06741).

**Checked:** Assumptions 1.1, 5.1–5.2, Theorems 5.3–5.4, and Corollary 5.9. The favorable one-dimensional regime uses the specified infinite-smoothness classes, nondegeneracy, initial integrability and controlled initial quadrature. It gives $O(\epsilon^{-1}|\log\epsilon|^3)$ work. No quantum precision-exponent advantage is claimed in that regime.

### BRTW24: QMC and MLQMC for mean-field particles

Nadhir Ben Rached, Abdul-Lateef Haji-Ali, Raúl Tempone and Leon Wilkosz. **Forward propagation of low discrepancy through McKean–Vlasov dynamics: from QMC to MLQMC.** Preprint, [arXiv:2409.09821v1](https://arxiv.org/abs/2409.09821) (2024). [Inspected full manuscript](https://arxiv.org/pdf/2409.09821).

**Checked:** Sections 3–5. Additive-noise weak convergence is proved; enhanced rates are proved for a simplified deterministic mean-field ODE. Favorable stochastic rates used for Section 5 complexity are conjectures. The paper's own MLQMC analysis gives quadratic cost under those conjectures, not a universal classical lower bound. Its scalar-state presentation does not supply a general high-dimensional weighted-QMC guarantee.

### BHNU18: Regression control variates

Denis Belomestny, Stefan Häfner, Tigran Nagapetyan and Mikhail Urusov. **Variance reduction for discretised diffusions via regression.** *Journal of Mathematical Analysis and Applications* **458**(1), 393–418 (2018). [DOI: 10.1016/j.jmaa.2017.09.002](https://doi.org/10.1016/j.jmaa.2017.09.002). [Inspected manuscript, arXiv:1510.03141v4](https://arxiv.org/pdf/1510.03141).

**Checked:** Approximation assumptions and Theorems 6.1–6.2. The subquadratic ordinary-SDE result pays for training and regression, with dimension-dependent approximation costs. Its favorable high-order regime is not implied by the fixed finite derivatives assumed here. Mean-curve solution would be an additional step.

### AAPdR26: A global mean-curve optimization approach

Ankush Agarwal, Andrea Amato, Stefano Pagliarani and Gonçalo dos Reis. **Numerical approximation of McKean–Vlasov SDEs via stochastic gradient descent.** *Advances in Applied Probability*, online first, 1–33 (2026). [DOI: 10.1017/apr.2026.10057](https://doi.org/10.1017/apr.2026.10057). [Inspected manuscript, arXiv:2310.13579v3](https://arxiv.org/html/2310.13579v3).

**Checked:** Mean-curve parameterization, Theorem 3.2 and numerical interpolation choices. Global temporal basis representations, including the Chebyshev experimental choice, are prior art. The stated nonconvex optimization result does not by itself certify global epsilon-cost for the desired fixed point.

### CDRGT15: Mean-field cubature

Pierre-Emmanuel Chaudru de Raynal and Camilo A. Garcia Trillos. **A cubature based algorithm to solve decoupled McKean–Vlasov forward–backward stochastic differential equations.** *Stochastic Processes and their Applications* **125**, 2206–2255 (2015). [DOI: 10.1016/j.spa.2014.11.018](https://doi.org/10.1016/j.spa.2014.11.018). [Author manuscript](https://discovery.ucl.ac.uk/id/eprint/1464371/1/de%20Raynal_Cubature_MKVFBSDE_v3.pdf).

Deterministic cubature is an established alternative. Time-step accuracy must be combined with the actual branching/recombination cost before inferring end-to-end work.

### CM19: Cubature with smooth scalar interaction

Dan Crisan and Eamon McMurray. **Cubature on Wiener space for McKean–Vlasov SDEs with smooth scalar interaction.** *Annals of Applied Probability* **29**(1), 130–177 (2019). [DOI: 10.1214/18-AAP1407](https://doi.org/10.1214/18-AAP1407). [Primary journal issue](https://www.imstat.org/publications/aap/aap_29_1/aap_29_1.pdf).

The smoothness and nondegeneracy/Hörmander analysis, including Lagrange interpolation of mean curves, is relevant prior art. A general dimension-controlled work bound does not follow merely from its order of convergence.

### AP21: Structured Fourier/Picard solution

Ankush Agarwal and Stefano Pagliarani. **A Fourier-based Picard-iteration approach for a class of McKean–Vlasov SDEs with Lévy jumps.** *Stochastics* **93**(4), 592–624 (2021). [DOI: 10.1080/17442508.2020.1771337](https://doi.org/10.1080/17442508.2020.1771337). [Author manuscript, arXiv:1812.05026](https://eprints.gla.ac.uk/215023/7/215023.pdf).

Its affine-in-state structure allows characteristic-function methods. Such an exact reduction should be used when available; the nonlinear state force in the polymer example is not affine.

### BS20: Regression on interacting particles

Denis Belomestny and John Schoenmakers. **Optimal stopping of McKean–Vlasov diffusions via regression on particle systems.** *SIAM Journal on Control and Optimization* **58**(1), 529–550 (2020). [DOI: 10.1137/18M1195590](https://doi.org/10.1137/18M1195590). [arXiv:1806.09483](https://arxiv.org/abs/1806.09483).

Used as evidence that particle-regression methods are established. Its optimal-stopping task is not assigned an unchanged complexity for the expectation problem here.

### CO19: Stochastic grid bundling

Ki Wai Chau and Cornelis W. Oosterlee. **Stochastic grid bundling method for backward stochastic differential equations.** *International Journal of Computer Mathematics* **96**(11), 2272–2301 (2019). [DOI: 10.1080/00207160.2019.1658868](https://doi.org/10.1080/00207160.2019.1658868). [Author manuscript](https://ir.cwi.nl/pub/28852/28852.pdf).

Local regression and conditional-expectation approximation are credible practical alternatives. No general mean-field inverse-linear work theorem is attributed to this paper.

### RSN21: Tensor-train backward/PDE representations

Lorenz Richter, Leon Sallandt and Nikolas Nüsken. **Solving high-dimensional parabolic PDEs using the tensor train format.** *Proceedings of the 38th International Conference on Machine Learning*, PMLR **139**, 8998–9009 (2021). [Primary proceedings page and paper](https://proceedings.mlr.press/v139/richter21a.html).

Demonstrates a serious representation-based classical alternative. Useful empirical ranks do not establish uniform rank bounds for every coefficient program considered here.

## 3. Polymer application and its classical reductions

### FO14: Dynamic polymer mean-field theory

Glenn H. Fredrickson and Henri Orland. **Dynamics of polymers: a mean-field theory.** *Journal of Chemical Physics* **140**(8), 084902 (2014). [DOI: 10.1063/1.4865911](https://doi.org/10.1063/1.4865911). [Inspected manuscript, arXiv:1312.7518](https://arxiv.org/pdf/1312.7518).

Provides the physical dynamic single-chain mean-field setting and discussion of additional closures. The supplied finite Fourier potential and Gaussian quench in this repository are explicit specializations, not an asserted model reduction theorem from this paper.

### ES13: A rigorous polymer mean-field limit

Weinan E and Hao Shen. **Mean field limit of a dynamical model of polymer systems.** *Science China Mathematics* **56**(12), 2591–2598 (2013). [DOI: 10.1007/s11425-013-4713-y](https://doi.org/10.1007/s11425-013-4713-y). [Primary author-institution record](https://collaborate.princeton.edu/en/publications/mean-field-limit-of-a-dynamical-model-of-polymer-systems/).

The primary abstract describes a Brownian elastic-string mean-field limit under regularity assumptions. **Abstract/metadata verification only.** It is cited for the model family, not as a proof of our finite-chain specialization or algorithm.

### DM06: An established classical SCMF implementation

Kostas Ch. Daoulas and Marcus Müller. **Single chain in mean field simulations: quasi-instantaneous field approximation and quantitative comparison with Monte Carlo simulations.** *Journal of Chemical Physics* **125**(18), 184904 (2006). [DOI: 10.1063/1.2364506](https://doi.org/10.1063/1.2364506). [Primary article abstract indexed in PubMed](https://pubmed.ncbi.nlm.nih.gov/17115792/).

The abstract describes field-based single-chain simulations and quantitative Monte Carlo comparisons. This supports treating SCMF as a strong classical method, not claiming exactness for arbitrary nonequilibrium reductions.

### G26: A recent Gaussian-closure alternative

Andriy Goychuk. **Dynamical mean-field theory of active, self-interacting heteropolymers.** *Physical Review E* **114**, 045406 (2026). [DOI: 10.1103/g3mp-mqmt](https://doi.org/10.1103/g3mp-mqmt). [Primary article](https://journals.aps.org/pre/abstract/10.1103/g3mp-mqmt). Earlier preprint: [arXiv:2604.02085](https://arxiv.org/abs/2604.02085), under a different title.

The primary abstract supplies a current Gaussian-closure comparator. This is an approximation for related polymer dynamics, not an exact solver for the full joint finite-chain equation here.

## 4. Function-valued quantum integration antecedent

### W04: Quantum parametric integration

Carsten Wiegand. **Quantum complexity of parametric integration.** *Journal of Complexity* **20**(1), 75–96 (2004). [DOI: 10.1016/j.jco.2003.10.001](https://doi.org/10.1016/j.jco.2003.10.001). [Inspected manuscript, arXiv:quant-ph/0305103v2](https://arxiv.org/pdf/quant-ph/0305103).

**Checked:** Theorem 2.9 and the multilevel interpolation construction. Function-valued quantum integration and multiscale interpolation are established. This source studies joint finite-dimensional smoothness in a query model; it does not itself supply the stochastic-path, feedback, preparation or reversible-circuit accounting of the current construction.

### KU26: Quantum integration using Faber–Schauder decompositions

Bernd Käßemodel and Tino Ullrich. **Quantum integration in Besov spaces with mixed smoothness.** [arXiv:2610.17204v1](https://arxiv.org/html/2610.17204v1), 12 October 2026.

**Checked:** Sections 4–5. The tensor-product Faber–Schauder decomposition and multilevel quantum summation are explicit antecedents. This work concerns integration of periodic functions with mixed smoothness under function-query access. It does not furnish the unknown mean-field curve, the stochastic-path mixed bound or the paid finite-bit path sampler here. The basis and its combination with quantum summation are not claimed as new.

## 5. Quantum estimation and mean-field simulation

### MY26: The closest quantum mean-field algorithm

Koichi Miyamoto. **Quantum algorithm for solving McKean-Vlasov stochastic differential equations.** *EPJ Quantum Technology* **13**, article 3 (2026), published online 10 December 2025. [DOI and full text](https://doi.org/10.1140/epjqt/s40507-025-00453-z). [Preprint, arXiv:2507.10926](https://arxiv.org/abs/2507.10926).

**Checked:** Task definition and Section 3's sequential-time construction and complexity. This is direct prior art for finite-moment feedback and a terminal expectation. Its stated step-query bound is $O(\epsilon^{-1-2/p})$ for weak order $p\in(1,2]$, including repeated preparation of the evolving distribution. The present construction obtains a different bound for the narrower additive-noise class using global feedback reconstruction and charged sample circuits. It does not improve every model in Miyamoto's more general diffusion setting. His comparison to a cubic particle baseline is not used as our strongest classical benchmark.

### ALLMSW21: Quantum multilevel Monte Carlo

Dong An, Noah Linden, Jin-Peng Liu, Ashley Montanaro, Changpeng Shao and Jiasu Wang. **Quantum-accelerated multilevel Monte Carlo methods for stochastic differential equations in mathematical finance.** *Quantum* **5**, 481 (2021). [Primary published paper](https://quantum-journal.org/papers/q-2021-06-24-481/). [Inspected preprint, arXiv:2012.06283v2](https://arxiv.org/pdf/2012.06283v2).

**Checked:** Theorem 2 and the SDE specialization in Theorem 3. Quantum MLMC and the inverse-linear accuracy behavior of suitable known-coefficient SDE expectations are inherited. The feedback curve is unknown in the present task; supplying it without additional accuracy powers is the extra construction to be proved. Our level-by-level sum also includes finite-bit perturbations and the mandatory minimum number of calls.

### M15: Mean estimation from a second-moment bound

Ashley Montanaro. **Quantum speedup of Monte Carlo methods.** *Proceedings of the Royal Society A* **471**, 20150301 (2015). [DOI](https://doi.org/10.1098/rspa.2015.0301). [Inspected author preprint, arXiv:1504.06987v3](https://arxiv.org/pdf/1504.06987v3), revised 11 July 2017.

**Checked:** Section 2.2, Algorithm 2 and Lemma 4; Section 2.3 distinguishes the variance-based centered version. The bounded-second-moment route is sufficient here. The explicit signed dyadic-bin implementation in the circuit chapter pays a harmless additional logarithmic factor and retains a comparator-based finite-bit encoding.

### BHMT02: Amplitude estimation

Gilles Brassard, Peter Høyer, Michele Mosca and Alain Tapp. **Quantum Amplitude Amplification and Estimation.** *Contemporary Mathematics* **305**, 53–74 (2002). [DOI](https://doi.org/10.1090/conm/305/05215). [Inspected primary preprint, arXiv:quant-ph/0005055v1](https://arxiv.org/pdf/quant-ph/0005055), 2000.

**Checked:** Theorem 12. Its probability-dependent error and constant success probability underlie the explicit bin estimator. The sample unitary, its inverse, marking, reflections, confidence amplification and gate synthesis are all counted in this repository.

## 6. Pseudorandomness and reversible implementation

### N92: Space-bounded pseudorandomness

Noam Nisan. **Pseudorandom Generators for Space-Bounded Computation.** *Combinatorica* **12**, 449–461 (1992). [DOI](https://doi.org/10.1007/BF01305237). [Author-institution record](https://cris.huji.ac.il/en/publications/pseudorandom-generators-for-space-bounded-computation-13/).

The foundational space-bounded PRG principle is inherited. The explicit inequality instantiated here is taken from the checked HashPRG analysis below, rather than attributed to an unspecified unit-cost PRG oracle.

### KPTW24: Explicit block-addressable HashPRG

Praneeth Kacham, Rasmus Pagh, Mikkel Thorup and David P. Woodruff. **Pseudorandom Hashing for Space-bounded Computation with Applications in Streaming.** [Inspected preprint, arXiv:2304.06853v2](https://arxiv.org/html/2304.06853v2), 4 January 2024.

**Checked:** Sections 3.1–3.2 and equation (1). The memory chapter specializes the construction to binary branching, supplies a concrete pairwise-independent Toeplitz hash circuit, and substitutes parameters that preserve both payoff moments. Only a classical finite-state output distribution is approximated. No quantum indistinguishability or cryptographic hardness is assumed.

### B89: Reversible checkpointing

Charles H. Bennett. **Time/Space Trade-Offs for Reversible Computation.** *SIAM Journal on Computing* **18**, 766–776 (1989). [DOI](https://doi.org/10.1137/0218053). [Author-institution record](https://research.ibm.com/publications/timespace-trade-offs-for-reversible-computation).

The checkpointing principle is inherited. The memory chapter gives its own clean XOR-bridge proof for fixed-arity circuit segments, including the exact recurrence and workspace. This avoids importing a Turing-machine time bound without specifying memory access and the actual leaf circuits.

### DN06: Finite-set gate synthesis

Christopher M. Dawson and Michael A. Nielsen. **The Solovay–Kitaev algorithm.** *Quantum Information and Computation* **6**(1), 81–95 (2006). [Primary preprint, arXiv:quant-ph/0505030](https://arxiv.org/abs/quant-ph/0505030).

The constructive synthesis result supplies polynomial-logarithmic gate and classical compilation overhead in inverse operator accuracy. The repository allocates the cumulative gate error across the actual predetermined calls; it does not assume exact arbitrary rotations.

## 7. The initial algorithmic lead

[**Algorithms and compilation**](https://github.com/GoGoKo699/Quantum-Off-the-Wall/blob/c2fd6333a0d2f19fd109a10b9baaba6b5c138923/science/algorithms-and-compilation.md), *Quantum Off the Wall*, pinned commit `c2fd6333a0d2f19fd109a10b9baaba6b5c138923`, supplied the initial leads toward quantum numerical algorithms and end-to-end compilation costs. It is a research map, not a proof of the construction in this repository. The primary papers above carry the specific imported results.

## Audit boundaries

The most consequential numerical rates are tied above to checked theorem locations or marked as derived here. Abstract-only items supply context and do not carry a proof obligation. None of the cited numerical experiments is a lower bound on classical computation. The illustrative polymer calculation is likewise reported with its actual mesh, particle counts and methods, without a certified total-error or physical-crossover claim.
