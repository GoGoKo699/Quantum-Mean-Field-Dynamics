# Memory, pseudorandomness, and reversible recomputation

[Home](../README.md) · [Model and results](model-and-results.md) · [Finite precision and circuits](finite-precision-and-circuits.md)

The direct implementation has near-linear precision dependence in logical gates, but it can retain a long random tape and path history. This chapter gives a constructive alternative. It combines a classical space-bounded pseudorandom generator with reversible checkpointing. The pseudorandom generator is used only to approximate the distribution of a classical finite-bit payoff. A quantum algorithm then estimates the mean of that deterministic payoff over a uniform short seed.

For the explicit ternary checkpoint schedule, put

```math
a_3=\log_3 5=1.464973\ldots\lt \tfrac32.
```

At fixed admissible model parameters, the alternative has

```math
G_{\rm compact}=\widetilde O(\epsilon^{-a_3}),
\qquad
Q_{\rm compact}=\mathrm{poly}\!\left(
\vert P\vert ,d,K,q,p,\log(1/\epsilon),\log(1/\alpha)
\right),
```

where $`p`$ is the charged working word length, including integer ranges and numerical condition factors. At fixed model parameters, $`p`$ grows only logarithmically in inverse accuracy and confidence. Thus $`\widetilde O(\epsilon^{-3/2})`$ is a simpler sufficient gate bound. The result counts the Boolean programs, random-bit generation, reversible simulation, amplitude estimation, gate synthesis, classical compilation and output specified below. It retains the separately charged classical moment curve. All polynomial-time input-program requirements are the explicit ones in [Finite precision and circuits](finite-precision-and-circuits.md#12-supplied-programs-rather-than-free-function-oracles). Large physical constants can enlarge $`p`$; no polynomial bound in their binary description length alone is inferred.

The PRG ingredient is inherited from Nisan and the explicit HashPRG analysis of Kacham, Pagh, Thorup and Woodruff. The checkpoint construction is a circuit-level version of Bennett's reversible simulation. The application-specific work is the compatibility of these ingredients with the random-time multilevel corrections and their first- and second-moment budgets. Neither inherited ingredient is claimed as new.

## 1. The finite sampler as a sequence of small configurations

Fix the previous classical feedback curve, a temporal coefficient, and a coupled path level. Section 3 of [Finite precision and circuits](finite-precision-and-circuits.md#3-the-exact-finite-grid-correction-program) is a fixed-clock streaming program. It stores two state vectors, a coarse increment accumulator, the random endpoint, the two frozen moment arguments, the eventual payoff, and a bounded amount of local arithmetic state. Initial sampling and payoff evaluation are included. Let

```math
N=O(n_l+1),\qquad
S=S(P,d,K,q,p,\log n_l)
```

bound its number of macrosteps and its number of bits retained between macrosteps. A macrostep has a Boolean circuit of size $`C_{\rm step}`$ and local scratch $`W_{\rm step}`$. These are charged quantities. For the explicit coefficient families considered here, they are polynomial in model size and precision. A safe bound is $`W_{\rm step}\leq O(C_{\rm step})`$.

Each random macrostep consumes one fresh block containing enough bits for its prescribed normal draws or endpoint draw. An unused suffix of a block is ignored. We may pad to $`N_R\leq O(N)`$ random blocks; the initialization can use a fixed number of additional blocks, or a larger first block. The entire sequence is padded to its worst-case clock. Inactive Euler steps are computed with zero masks and still incur their allotted work.

Two different resources must be distinguished. A fresh random block is an **input symbol** of the finite-state machine used to prove the PRG guarantee. It is not part of that machine's retained state between symbols. In the circuit implementation, regenerating and evaluating this symbol does use real gates and scratch qubits; those are included in $`C_{\rm step}`$ and $`W_{\rm step}`$. In particular the choice of PRG block length below has no circular requirement that it exceed three times its own length.

For the distributional proof, append a clock of $`\lceil\log_2(N_R+1)\rceil`$ bits to the retained configuration. This converts the known sequence of transition functions into one finite-state transition function over an alphabet of random blocks. If $`S`$ already includes that clock, do not count it again. A valid state-width bound is

```math
w\geq S+\lceil\log_2(N_R+1)\rceil+2.
```

The final payoff bits are part of this state. The transition function is implemented by the supplied arithmetic programs and the fixed feedback schedule. It is **never** constructed by enumerating its $`2^w`$ states or its full transition table. Its existence as a finite function is enough for the PRG theorem; its explicit step circuit supplies the implementation cost.

## 2. An explicit block-addressable classical PRG

### 2.1 Binary HashPRG with pairwise independent hashes

Let $`k=\lceil\log_2\max\lbrace 2,N_R\rbrace \rceil`$ and choose a block length $`b`$ below. The construction outputs $`2^k`$ blocks of $`b`$ bits. Let $`x\in\lbrace 0,1\rbrace ^b`$ be uniform. For every depth $`i=0,\ldots,k-1`$, choose two independent pairwise independent hash functions

```math
h_i^{(0)},h_i^{(1)}:\lbrace 0,1\rbrace ^b\longrightarrow\lbrace 0,1\rbrace ^b.
```

All choices are included in the uniformly random seed. Define $`G_0(x)=x`$ and recursively concatenate

```math
G_{i+1}(x)=
G_i\!\left(h_i^{(0)}(x)\right)
\;\Vert\;
G_i\!\left(h_i^{(1)}(x)\right),
```

where each $`G_i`$ uses the already chosen lower-level hashes. If $`j`$ has binary digits $`j_{k-1}\cdots j_0`$, its block is

```math
G_k(x)[j]
=h_0^{(j_0)}\!\left(h_1^{(j_1)}\!\left(
\cdots h_{k-1}^{(j_{k-1})}(x)\cdots\right)\right).
```

The block address is the known classical macrostep index. A requested block requires $`k`$ hash evaluations; generating all earlier blocks is unnecessary. This is the binary branching case of [HashPRG, Section 3.1](https://arxiv.org/html/2304.06853v2#S3.SS1).

Here is a concrete hash family with no field-construction subroutine. For $`z\in\lbrace 0,1\rbrace ^b`$, use

```math
h_{a,c}(z)=T_a z\oplus c,
```

where $`T_a`$ is a $`b\times b`$ binary Toeplitz matrix, its $`2b-1`$ diagonals are independent uniform bits, and $`c`$ is an independent uniform $`b`$-bit vector. The seed length of one hash is $`3b-1`$. It is pairwise independent: for nonzero $`z`$, $`T_a z`$ is uniform. To see this, take the largest index of a nonzero coordinate of $`z`$. The corresponding $`b`$ diagonal variables form a triangular linear system with unit diagonal in the $`b`$ output coordinates, after the other diagonals are fixed. Adding independent $`c`$ makes $`h(u)`$ uniform independently of $`h(u)\oplus h(v)=T_a(u\oplus v)`$ for distinct $`u,v`$.

Schoolbook Toeplitz multiplication uses $`O(b^2)`$ binary AND/XOR operations. A block evaluation uses $`O(kb^2)`$ operations. Compute-copy-uncompute gives a reversible block evaluator with $`O(kb^2)`$ gates and a conservative $`O(kb^2)`$ scratch bound. The complete seed has exactly

```math
r_{\rm seed}=b+2k(3b-1)
```

bits. Preparing its uniform superposition uses that many Hadamard gates.

### 2.2 An explicit total-variation guarantee

Let $`Q`$ be any finite-state machine with at most $`2^w`$ states over the alphabet $`\lbrace 0,1\rbrace ^b`$. Write $`Q(D)`$ for its matrix of final-state transition probabilities under input distribution $`D`$. The matrix norm here is the maximum row $`\ell^1`$ norm. Equation (1) of [HashPRG, Section 3.2](https://arxiv.org/html/2304.06853v2#S3.SS2) states, in the binary case, that for any $`\theta\gt 0`$,

```math
\Pr_h\!\left\lbrace 
\left\Vert Q(G_k(*,h))-Q(U_b^{\otimes 2^k})\right\Vert \gt \theta
\right\rbrace 
\leq
\frac{(2^k-1)^2k}{\theta^2}\,2^{3w-b}.
```

The probability is over the hash functions; the star denotes a uniform initial block $`x`$. The hash family may be any pairwise independent family, including the explicit Toeplitz family above.

For a desired $`0\lt \chi\leq1/4`$, choose the integer block length at least the number of random bits needed in one macrostep and at least

```math
b\geq
3w+2k+
\left\lceil\log_2(16(k+1))\right\rceil+
3\left\lceil\log_2(1/\chi)\right\rceil+2.
```

Set $`\theta=\chi/2`$. Substitution bounds the bad-hash probability by $`\chi/16`$. On the complementary event, final-state total variation is at most $`\theta/2=\chi/4`$. Averaging over all hash seeds therefore gives

```math
\mathrm{TV}\!\left(
\mathrm{state}(G_k(U_{r_{\rm seed}})),
\mathrm{state}(U_b^{\otimes 2^k})
\right)
\leq\chi/4+\chi/16\lt \chi.
```

There is no search for a good fixed hash family: all hash bits and $`x`$ remain uniformly random. Nor is an exponential transition matrix materialized. The finite-state proof supplies a bound valid for every permitted transition circuit of this width. The actual block evaluator above remains polynomial in $`w,k,\log(1/\chi)`$.

### 2.3 First and second moments must both be preserved

Let $`\widehat D_l`$ be the finite-bit correction from the preceding chapter. Suppose its proved bounds are

```math
\vert \widehat D_l\vert \leq B_l,
\qquad
\Vert \widehat D_l\Vert _2\leq\sigma_l.
```

Choose $`B_l`$ by a deterministic all-seed range certificate, not by a high-probability bound. In the bounded-generator core one can take $`B_l=2H(B_F+u_F)+u_{\rm out}`$. The affine linear-growth extension instead uses its explicit finite-state range, and $`B_l`$ may grow polynomially with $`n_l`$. That growth is allowed here.

Write $`D_l^{\rm PRG}`$ for the deterministic finite program evaluated on HashPRG output. Total variation contracts under taking the final payoff. Because the payoff has range $`[-B_l,B_l]`$ and its square has range $`[0,B_l^2]`$,

```math
\left\vert \mathbb E D_l^{\rm PRG}-\mathbb E\widehat D_l\right\vert 
\leq2B_l\chi_l,
\qquad
\left\vert \mathbb E(D_l^{\rm PRG})^2-
\mathbb E\widehat D_l^{\,2}\right\vert 
\leq B_l^2\chi_l.
```

For the coefficient budget of the preceding chapter, set

```math
\rho=\frac{\delta}{16(L+1)},
\qquad
0\lt \chi_l\leq
\min\left\lbrace 
\frac14,\frac{\rho}{2B_l},\frac{\sigma_l^2}{B_l^2}
\right\rbrace .
```

Choose a dyadic $`\chi_l`$ below this minimum. A certified $`B_l=0`$ or $`\sigma_l=0`$ permits the identically zero correction to be omitted. Otherwise the construction gives

```math
\left\vert \mathbb E D_l^{\rm PRG}-\mathbb E\widehat D_l\right\vert \leq\rho,
\qquad
\Vert D_l^{\rm PRG}\Vert _2\leq\sqrt2\,\sigma_l.
```

Thus the summed additional coefficient bias is at most $`\delta/16`$. The multilevel second-moment estimate is retained with a constant factor. In particular one may use

```math
\sigma_l\leq S_kHh_l+\rho
```

for corrections and the analogous bounded base-level norm. A purely additive PRG mean guarantee would be insufficient: it could increase the correction variance and destroy the multilevel cost bound.

For fixed model parameters all tolerances and all finite ranges in this construction are polynomial in $`\epsilon`$ or $`1/\epsilon`$, up to logarithms. Consequently $`\log(1/\chi_l)`$ and $`b`$ are polynomial in model size and logarithmic accuracy parameters. An all-seed state magnitude of order $`\sqrt{n_l\log n_l}`$ in the linear-growth extension does not alter the precision exponent; it changes $`\log B_l`$ and the explicit constants in the seed budget.

### 2.4 What the guarantee says about quantum use

The PRG guarantee concerns the original **classical, one-pass, finite-bit sample program**. After replacing its random tape, its output is a deterministic function of a short uniform seed. A reversible circuit computes exactly that function. Quantum mean estimation is applied to its actual seed distribution, with the moment bounds just proved.

No claim that HashPRG fools quantum distinguishers is made or needed. The later coherent replay of a seed does not modify the function whose classical distribution was analyzed. Likewise, storing the classical feedback curve outside the sampler does not turn the PRG claim into a quantum-memory oracle: for each fixed prior curve, its fixed transition functions are included in the finite-state theorem. Conditional guarantees therefore remain valid through the adaptive classical iteration.

## 3. Explicit reversible checkpoint circuits

### 3.1 A clean bridge for one irreversible step

After the PRG substitution, each macrostep is a total deterministic Boolean function

```math
f_i:\lbrace 0,1\rbrace ^S\times\lbrace 0,1\rbrace ^{r_{\rm seed}}
\longrightarrow\lbrace 0,1\rbrace ^S.
```

Its index $`i`$ is known in the classical circuit schedule. The seed is read-only. Arithmetic is a finite Boolean map, with overflow behavior defined on every bit string; the analytic range certificate proves that valid trajectories do not overflow. This makes $`f_i`$ a total function even on unreachable configurations.

Compute the Boolean circuit for $`f_i`$ into fresh temporary wires, XOR its output into a target configuration register, and reverse the temporary computation. The result is the exact unitary

```math
B_i:\ \vert x\rangle\vert y\rangle\vert s\rangle
\longmapsto \vert x\rangle\vert y\oplus f_i(x;s)\rangle\vert s\rangle,
```

with clean temporary wires. It costs $`O(C_{\rm step})`$ gates and $`O(C_{\rm step})`$ temporary qubits, where $`C_{\rm step}`$ now includes $`O(kb^2)`$ PRG block regeneration. Include at least $`S`$ copying gates in this bound. The operation is its own inverse. It does not require $`f_i`$ to be invertible.

### 3.2 The ternary schedule

Pad the number of macrosteps to $`3^r`$, with $`3^r\lt 3\max\lbrace 1,N\rbrace `$, using identity transitions. For a segment of length $`3^r`$, let $`F`$ be the composition of its irreversible steps. We construct a clean bridge

```math
B_F:\ \vert x\rangle\vert y\rangle
\longmapsto \vert x\rangle\vert y\oplus F(x)\rangle.
```

For one step this is the construction above. For a longer segment write $`F=F_3\circ F_2\circ F_1`$ for its three equal consecutive subsegments. Starting with two zero checkpoint registers $`a,b`$, perform

1. $`B_{F_1}(x,a)`$;
2. $`B_{F_2}(a,b)`$;
3. $`B_{F_3}(b,y)`$;
4. $`B_{F_2}^{-1}(a,b)`$;
5. $`B_{F_1}^{-1}(x,a)`$.

After the third operation, $`y`$ has been XORed with $`F_3(F_2(F_1(x)))`$. The last two operations erase $`b`$ and $`a`$ in that order without changing $`x,y`$. Induction proves the claimed clean bridge for every input basis state and therefore for arbitrary superpositions. A checkpoint is a coherently retained computational encoding, copied with CNOTs; the construction does not copy an arbitrary unknown quantum state.

There are exactly $`5^r`$ leaf bridge calls. At most two intermediate $`S`$-bit registers per recursion level are live. Consequently

```math
G_{\mathrm{bridge}}
\leq O(C_{\rm step})\,5^r
\leq O(C_{\rm step})(3N)^{\log_3 5},
```

and

```math
Q_{\mathrm{bridge}}
\leq (2r+2)S+r_{\rm seed}+O(C_{\rm step}).
```

The displayed gate and workspace bounds include the regeneration of the appropriate pseudorandom block every time a leaf step is replayed. No long random tape or path history is retained.

Prepare the seed uniformly, apply the full bridge from a fixed zero initial configuration into a zero target configuration, copy the payoff bits to a designated output register, and apply the inverse bridge. This produces a clean sampler retaining only the seed and payoff, with the same asymptotic bounds. Initial-state and endpoint generation are ordinary macrosteps of this bridge. The sampler's inverse is the reversed gate list.

### 3.3 A general fixed tradeoff

The same proof with a fixed branching factor $`r_0\geq2`$ uses $`r_0-1`$ forward checkpoint segments, one final segment, and $`r_0-1`$ reverse segments. Its recurrence and workspace are

```math
\begin{aligned}
C(r_0 n)&=(2r_0-1)C(n),\\
a_{r_0}&=\log_{r_0}(2r_0-1),\\
Q&=O((r_0-1)S\log_{r_0}N)+r_{\rm seed}+O(C_{\rm step}).
\end{aligned}
```

For any **specified fixed** $`\kappa\gt 0`$, choose an integer $`r_0\geq\max\lbrace 2,2^{1/\kappa}\rbrace `$. Then $`a_{r_0}\lt 1+\kappa`$. This proves the general $`\widetilde O(\epsilon^{-1-\kappa})`$ corollary directly in the circuit model. Constants grow with $`r_0`$, so taking $`\kappa`$ to zero as a function of accuracy is not part of the claim. The explicit ternary choice is the concrete resource statement used here.

This construction implements Bennett's checkpointing principle without transferring a Turing-machine time bound to an unspecified quantum memory model. Each leaf is an actual finite Boolean circuit, and each replayed segment has an explicit deterministic gate schedule.

## 4. Uniform compilation and classical feedback storage

The checkpoint recursion decides all segment boundaries and leaf indices classically. A depth-first schedule generator uses a classical recursion stack of depth $`O(\log N)`$ and emits the leaf bridge circuits in the stated order. Inverse segments are emitted in the reverse order with inverse gates. Schedule generation is linear in the emitted gate list, apart from the explicitly charged polynomial-bit coefficient evaluations and synthesis.

At a known leaf index $`i`$, the required feedback constant is $`v(i h_l)`$, or the known coarse-grid counterpart. With ordinary classical addressable memory, binary search among the $`M`$ stored curve nodes and finite interpolation use $`O(\log M)`$ comparisons and polynomial-bit arithmetic. Each addressed read transfers its actual word bits. This is classical memory used to emit a fixed circuit constant, not an operation on a superposition of addresses.

There is also a purely sequential compilation route, so the exponent does not depend on that classical access convention. At the beginning of one Picard iteration, compute the feedback tables for every path level $`0\leq l\leq L_{\max}`$ that can occur. For each level, merge the sorted list of $`M`$ curve nodes with the sorted grid, computing each required linear interpolation and rounding it inside $`D`$. Use the largest working precision needed by any call in the iteration. The total work is

```math
O\!\left(K\mathrm{poly}(p)
\left[M(L_{\max}+1)+n_{\max}\right]\right),
\qquad n_{\max}=2^{L_{\max}},
```

because $`\sum_l2^l\lt 2n_{\max}`$. The storage is $`O(K\mathrm{poly}(p)(M+n_{\max}))`$ bits. Both fit $`\widetilde O(\epsilon^{-1})`$ at fixed model parameters. Put the fine and coarse feedback constants for each macrostep in one fixed-size row of a classical tape, with harmless padding for initialization, finalization and identity steps. One table is reused by all temporal supports requiring that path level.

To avoid repeated travel between distinct grid tables, execute path levels in the outer loop within one Picard map. At level $`l`$, scan the list of temporal coefficients and features, estimate that level's mean for every coefficient whose predetermined cutoff satisfies $`L_{k,H}\geq l`$, and add it to a stored finite coefficient accumulator. Inactive entries are skipped during the scan. The cutoffs and confidence allocations are known before any estimate is made, and the level means enter only through a classical sum, so this reordering leaves the algorithm and error proof unchanged. Accumulator scans cost $`O(KM(L_{\max}+1)\mathrm{poly}(p))`$. The feedback-grid tapes are visited once in level order, with total inter-table travel $`O(n_{\max})`$. A coefficient accumulator is rounded only within its already allocated final or per-level arithmetic tolerance.

For an arity-$`r_0`$ bridge of padded length $`r_0^r`$, the sequence of leaf indices begins and ends at the first index of its interval. Its total adjacent-index travel $`D_r`$ satisfies

```math
D_0=0,\qquad
D_r=(2r_0-1)D_{r-1}+2(r_0-1)r_0^{r-1},
\qquad
D_r=2\left[(2r_0-1)^r-r_0^r\right].
```

The recurrence follows by adding the travel inside the $`2r_0-1`$ child bridges and the $`2(r_0-1)`$ transitions of one child-interval length between their first indices. Therefore a two-way classical tape delivers all feedback rows in $`O((2r_0-1)^r)`$ row moves, the same order as the leaf count. The recursive compiler emits each leaf circuit immediately. Repeated checkpoint visits incur only charged tape movement and polynomial-bit work; no transition table or random-access oracle is required. The endpoint-dependent choice of the final frozen feedback value remains the masked sequential accumulator inside the sample program.

Neither the quantum state nor the seed controls a memory address into the feedback table. The only quantum conditionals are Boolean gates on the explicitly represented masks and state bits. No RAM access is translated into a free quantum operation. The prior classical curve, coefficient descriptions and compiled circuit descriptions do occupy classical storage.

For the Lipschitz temporal reconstruction there are $`M=O(1/\eta)`$ stored nodes per moment, with model-dependent constants and finite word lengths. Classical reconstruction and explicit output of this curve cost at least their ordinary $`O(KM)`$ size. The curve and optional path-grid tables may be maintained in $`\widetilde O(K/\epsilon)`$ classical bits at fixed model size; a fully materialized quantum gate list is larger and has size bounded by its charged gate count. Alternatively a classical controller can stream the recursively generated schedule. The small-space theorem concerns **coherent quantum workspace**, not the elimination of these classical resources.

## 5. The multilevel and temporal cost sums

Let $`a=a_3\gt 1`$ for the ternary version, or $`a=a_{r_0}\gt 1`$ for another fixed choice. With $`h_l=2^{-l}`$ on the normalized horizon, a clean coupled sample and its inverse cost

```math
C_{\rm compact}(p,w,k,\chi_l)\,(h_l^{-1}+1)^a,
```

where $`C_{\rm compact}`$ is polynomial in the actual step-evaluation, Gaussian-generation and PRG parameters. For fixed admissible model parameters this contributes only powers of logarithms in accuracy and confidence. Its coherent workspace is polynomial in those same parameters and $`\log(1/h_l)`$.

The PRG-corrected norm bound is at most a constant times $`S Hh_l+\rho`$. Use the equal per-level mean accuracy $`e_l=\delta/[4(L+1)]`$ and the explicit mean-estimation interface from [Finite precision and circuits](finite-precision-and-circuits.md#6-an-explicit-second-moment-quantum-mean-estimator). Up to its displayed logarithms, the work for one coefficient is bounded by

```math
C_{\rm compact}
\sum_{l=0}^L h_l^{-a}
\left[1+\frac{(S Hh_l+\rho)(L+1)}{\delta}\right].
```

The base level is included by enlarging the supplied constant $`S`$. Since $`\rho(L+1)/\delta=1/16`$, its contribution is absorbed in the one-sample term. The two geometric sums are bounded by constants times

```math
h_L^{-a}
\quad\text{and}\quad
\frac{SH(L+1)}{\delta}\,h_L^{1-a}.
```

The analytic bias choice gives $`h_L^{-1}=O(1+HB_k/\delta)`$. This yields the sufficient fixed-model bound

```math
W(H,\delta)
=\widetilde O\!\left(1+(H/\delta)^a\right).
```

This estimate remains valid when the analytic bias vanishes and $`L=0`$: the variance-aware mean work $`O(1+H/\delta)`$ is then bounded by the right-hand side because $`a\gt 1`$. The mandatory cost of at least one finest-level sample has been retained throughout. Optimizing the statistical allocation is unnecessary for this exponent.

At temporal scale $`H=2^{-j}`$ there are $`O(1/H)`$ coefficients. For one reconstructed component their total cost is therefore at most

```math
\widetilde O\!\left(H^{-1}+H^{a-1}\delta^{-a}\right).
```

Summing scales down to $`H_{\min}=\Theta(\eta)`$ gives an output-size term $`O(\eta^{-1})`$ and the convergent geometric factor

```math
\sum_{j\geq0}2^{-j(a-1)}
=\frac{1}{1-2^{-(a-1)}}.
```

The coefficient accuracy is $`\delta=\Theta(\eta/[\sqrt K(J+1)])`$ with $`J=O(\log(1/\eta))`$. At fixed model parameters, $`\eta`$ is proportional to $`\epsilon`$ up to the explicit stability factors in [Model and results](model-and-results.md), and the number of Picard iterations is logarithmic in inverse accuracy. It follows that all mean-curve calls cost $`\widetilde O(\epsilon^{-a})`$ gates. The terminal expectation is one ordinary coupled hierarchy, with the same bound. Classical node reconstruction, initial-mean computation, numerical parameter selection and compilation fit the displayed cost under the stated polynomial-bit input contract.

## 6. Confidence, quantum gates, and the complete resource theorem

Finite arithmetic and PRG replacement are deterministic expectation biases already budgeted inside $`\delta`$; they are not additional failure events. For ideal quantum statistical estimates, allocate total failure probability at most $`\alpha/4`$ to all coefficient calls and at most $`\alpha/4`$ to the terminal call. Conditional bounds remain valid for an adaptively selected prior curve because every sampler and PRG bound above is uniform over the admissible fixed curves.

The actual quantum sample circuits consist of exact Boolean reversible gates, uniform-seed Hadamards and the amplitude-estimation control circuitry. Comparator-based amplitude encoding is exact on finite dyadic outputs. Approximate Fourier-transform rotations and any chosen finite-basis synthesis use the global hybrid-error budget of [Finite precision and circuits](finite-precision-and-circuits.md#72-the-only-approximation-in-the-quantum-gates). There are a priori finite worst-case call and gate counts. Choosing each synthesized rotation's operator error below $`\alpha/[16\max\lbrace 1,N_{\rm rot}\rbrace ]`$ suffices for an additional final-output total-variation allowance below $`\alpha/2`$. The constructive synthesis work and gate overhead are polynomial in logarithmic accuracy and confidence. Thus the actual output satisfies the claimed error event with probability at least $`1-\alpha`$.

**Theorem: compact coherent implementation.** Assume the analytic conditions of [Model and results](model-and-results.md), the explicit uniform polynomial-bit coefficient/initial-mean contract of [Finite precision and circuits](finite-precision-and-circuits.md), and its finite sampler and deterministic range bounds. Include polynomial-bit classical circuit generation in this contract, as in the preceding chapter's $`C_{\rm build}`$ and $`C_{\rm init}`$ bounds. For each specified fixed checkpoint arity $`r_0\geq2`$, there is a uniformly generated family of finite-basis quantum circuits, with adaptive classical Picard control, that returns the stated scalar output with additive error $`\epsilon`$ and success probability at least $`1-\alpha`$. At fixed model parameters, its logical gate and classical circuit-generation counts are

```math
\widetilde O(\epsilon^{-a_{r_0}}),
\qquad a_{r_0}=\log_{r_0}(2r_0-1),
```

and its coherent workspace is polynomial in the charged working word length, model evaluator sizes and logarithmic accuracy and confidence parameters. In particular it is polylogarithmic in inverse accuracy for each fixed admissible model. The finite classical curve and coefficient descriptions are separately stored and charged. All arithmetic and state-range factors remain explicit input-program costs; a model with superpolynomial evaluator complexity does not inherit the simplified precision statement.

For $`r_0=3`$, $`a_{r_0}\lt 3/2`$. Consequently even the conservative serial upper bound on logical qubit-time, $`Q_{\rm compact}G_{\rm compact}`$, has the same accuracy exponent $`a_3`$, up to logarithms and fixed-model factors. This addresses a coherent-memory objection to the direct implementation. It is not a physical-qubit crossover estimate, and it does not bound a different accounting rule that multiplies the persistent **classical** memory by the full runtime.

| Implementation | Logical gates, fixed model | Coherent workspace | Separately retained classical curve |
|---|---:|---:|---:|
| Direct history and true random tape | $`\widetilde O(\epsilon^{-1})`$ | Can be $`\widetilde O(\epsilon^{-1})`$ | $`\widetilde O(K/\epsilon)`$ bits |
| Binary HashPRG and ternary checkpoints | $`\widetilde O(\epsilon^{-\log_3 5})`$ | Fixed-model polylogarithmic accuracy dependence; evaluator workspace charged | $`\widetilde O(K/\epsilon)`$ bits |
| Simpler bound for the same ternary circuit | $`\widetilde O(\epsilon^{-3/2})`$ | Same | Same |

The direct option remains available when gates matter more than coherent workspace. No PRG assumption or hardness conjecture is needed for either implementation.

## 7. Verification and attribution

[The resource-construction checker](../verification/check_resource_construction.py) tests finite instances of the explicit constructions: the clean checkpoint bridge for noninvertible step maps, its exact leaf count and cleaned registers, block-addressed versus recursively generated HashPRG output, Toeplitz pairwise independence, the numerical parameter inequalities, and aligned endpoint masks. It does not simulate quantum amplitude estimation or validate a hardware crossover. The general guarantees are the proofs above and in the finite-precision chapter.

Primary sources for inherited results are:

- Noam Nisan, *Pseudorandom Generators for Space-Bounded Computation*, Combinatorica 12 (1992), 449–461. [Author-institution record](https://cris.huji.ac.il/en/publications/pseudorandom-generators-for-space-bounded-computation-13/).
- Praneeth Kacham, Rasmus Pagh, Mikkel Thorup and David P. Woodruff, *Pseudorandom Hashing for Space-bounded Computation with Applications in Streaming*, [arXiv:2304.06853v2](https://arxiv.org/html/2304.06853v2). Section 3.1 gives the construction and block formula; Section 3.2, Equation (1), gives the explicit error inequality used here. The binary specialization and concrete budget substitution are written above.
- Charles H. Bennett, *Time/Space Trade-Offs for Reversible Computation*, SIAM Journal on Computing 18 (1989), 766–776. [DOI](https://doi.org/10.1137/0218053) and [author-institution record](https://research.ibm.com/publications/timespace-trade-offs-for-reversible-computation). The fixed-arity checkpoint principle is inherited; the finite circuit recurrence, leaf schedule and small-space bound needed here have been proved explicitly.

The mean-estimation, amplitude-estimation and gate-synthesis sources are listed in [Finite precision and circuits](finite-precision-and-circuits.md#11-scope-and-provenance).
