# Finite precision, sample circuits, and quantum mean estimation

This chapter makes the sample-access assumption in [Model and results](model-and-results.md) constructive. It proves that finite random registers, arithmetic, confidence amplification, gate synthesis, classical compilation, and output rounding preserve the stated accuracy exponents. The construction is a uniform family of Boolean and quantum circuits under the input contract below. It does not assume an exact real-arithmetic machine, a Gaussian amplitude-loading oracle, or a quantum random-access memory.

The analytic inputs are the Dynkin identity and the coupled Euler estimates proved in [Dynkin details](dynkin-details.md). The mean-estimation primitive is inherited from Brassard–Høyer–Mosca–Tapp and Montanaro; an explicit implementation with a slightly larger logarithmic factor is given here. The optional smaller-workspace implementation is proved separately in [Memory and recomputation](memory-and-recomputation.md).

## 1. Bit model and the input contract

### 1.1 Time normalization and constants

Rescale a supplied horizon $T$ to $[0,1]$. The transformed drift and diffusion are $\bar a=T a$ and $\bar\Sigma=\sqrt T\Sigma$. Their evaluation and the computation of $\sqrt T$ are charged. All constants below refer to the transformed system. This is a change of units, not a uniform-in-$T$ assertion.

Write $q$ for the number of independent Brownian coordinates. The bounded core has

```math
|a(x,m)|\leq B_a,\qquad
|a(x,m)-a(y,n)|\leq L_x|x-y|+L_m|m-n|.
```

For one scalar mean feature, write its generator as $F(x,m)$ and assume the already established bounds

```math
|F(x,m)|\leq B_F,\qquad
|F(x,m)-F(y,n)|\leq L_{Fx}|x-y|+L_{Fm}|m-n|.
```

The state and moment norms are Euclidean unless another norm is explicitly substituted with its induced constants. Let $\nu=\|\bar\Sigma\|_F$ and $\nu_o=\|\bar\Sigma\|_{\rm op}$. The affine-drift extension used for polymers is treated in Section 9.

### 1.2 Supplied programs, rather than free function oracles

The ordinary classical input contains a finite binary description $P$ of the coefficients, initial law and observable, their numerical parameters, and certified bounds used in the analysis. It also specifies a uniform evaluator with the following contract.

Given a precision request $p$, a state range bound $R\geq1$, and binary inputs representing exact dyadic $x,m$ in the stated domain, a classical circuit generator emits fixed-clock Boolean circuits for $a$, the required $F$, and the terminal observable. The input also includes a certified evaluator for the initial moment vector: $g(x_0)$ in the deterministic case and $\mathbb E g(X_0)$ in the permitted random-initial-law extension. The output errors, measured in the stated scalar or vector norm, are at most $2^{-p}$. Let

```math
C_{\rm eval}(P,p,\lceil\log_2 R\rceil),\quad
W_{\rm eval}(P,p,\lceil\log_2 R\rceil),\quad
C_{\rm gen}(P,p,\lceil\log_2 R\rceil)
```

bound their Boolean operation count, working bits, and classical circuit-generation work. These are charged functions, not unit-cost queries. The precision-exponent conclusions apply when they are polynomial in the model description, $d,K,q,p$, and $\log R$. A coefficient implementation requiring more work retains that cost in the theorem.

The circuit basis may be AND, XOR and NOT, with fan-out represented by copying bits. Integer addition, multiplication, comparison and division can be implemented with schoolbook polynomial-bit algorithms. A bounded loop is unrolled or padded to its worst-case length. A data-dependent table lookup is compiled into actual Boolean gates and included in $C_{\rm eval}$; it is not granted a RAM or QRAM query cost.

The contract for $F$ includes derivative computation and the Hessian contraction. Smoothness alone does not imply cheap evaluation. The Fourier polymer formulas in [Polymer dynamics](polymer-dynamics.md) meet this contract explicitly; Section 10 gives the arithmetic construction.

The feedback curve is classical. Its dyadic coefficients or finite rational node values are stored explicitly. Every required $v(t_i)$ is computed classically and hardcoded in the circuit for deterministic grid index $i$. Rounding is chosen inside the prescribed moment box $D$: round each component and clip it to a nearest representable value in its interval. The extra grid-spacing error is included in $u_v$. This preserves the domain of the certified drift bounds even at a box boundary. The complete feedback table is not copied into a quantum register.

### 1.3 Initial states

The basic case is a specified dyadic $x_0$. An initial approximation to another computable $x_0$ is included in the numerical error budget. The supported Gaussian extension is $X_0=\mu+L_0Z_0$, with an ordinary explicit matrix or structured factor $L_0$. Its matrix arithmetic, storage and Gaussian generation are charged. The fine and coarse paths in a correction use the same initial draw.

If $Z_0$ has $q_0$ coordinates, use independent calls of Lemma 1 with coordinate error $\eta_0$, and approximations satisfying $|\widehat\mu-\mu|\leq u_\mu$, $\|\widehat L_0-L_0\|_{\rm op}\leq u_L$. With final transformation arithmetic error at most $u_{\rm init}$, an explicit initial coupling bound is

```math
\|\widehat X_0-X_0\|_{L^4}
\leq
u_\mu+\|L_0\|_{\rm op}\sqrt{q_0}\,\eta_0
u_L\sqrt{q_0}(3^{1/4}+\eta_0)+u_{\rm init}.
```

Its all-seed range is at most

```math
R_{X_0}
=|\widehat\mu|+
3\|\widehat L_0\|_{\rm op}\sqrt{q_0}\sqrt{s_0+1}
+u_{\rm init},
```

where $s_0$ is the normal cutoff parameter. These formulas determine the initial sampler's requested bits, rather than assuming exact Gaussian state preparation.

More generally, a supplied finite-bit initial sampler must provide a proved $L^2$ coupling error, or an $L^4$ coupling error for Section 9, at its charged time and space cost. This is a stated additional input property, not a claim about arbitrary sampleable laws. The polymer Gaussian sampler is constructed from the elementary Gaussian routine below. Its known initial feature means are evaluated classically with certified rounding error.

For every supported initial law, compute the initial moment vector once to the outer algorithm's Euclidean error $\eta/4$ and charge that evaluation. An accurate path sampler alone does not grant an inexpensive initial mean.

## 2. A finite-bit Gaussian generator with a coupling proof

The algorithm needs independent finite-bit classical random variables in superposition. It does not need a quantum state whose amplitudes are supplied Gaussian probabilities.

### Lemma 1: Gaussian generation

For each dyadic $0<\eta_Z<1$, there is an explicit deterministic fixed-clock Boolean circuit $\mathsf{Normal}_{\eta_Z}$ taking $2m$ uniform bits and returning a dyadic $\widehat Z$ such that a coupling to $Z\sim N(0,1)$ satisfies

```math
\|\widehat Z-Z\|_{L^4}\leq\eta_Z,
\qquad |\widehat Z|\leq 3\sqrt{s+1},
\qquad m,s=O(1+\log(1/\eta_Z)).
```

Its gate count, workspace and classical generation time are polynomial in $1+\log(1/\eta_Z)$. In particular the same coupling bound holds in $L^2$.

**Construction and proof.** Let $U,V$ be independent continuous uniform variables. The Box–Muller identity [BM58] gives a standard normal

```math
Z=\sqrt{-2\log U}\cos(2\pi V).
```

Choose an integer $s\geq1$ such that

```math
4\,2^{-s}(s^2+2s+2)\leq(\eta_Z/3)^4,
```

and put $a=2^{-s}$. Such an $s$ is found by incrementing an integer and checking a rational inequality; its size is $O(\log(1/\eta_Z))$. Set $U_a=\max(U,a)$ and $Z_a=\sqrt{-2\log U_a}\cos(2\pi V)$. Since the two outputs differ only for $U<a$,

```math
\mathbb E|Z-Z_a|^4
\leq\int_0^a(-2\log u)^2\,du
=4a\left[(\log(1/a))^2+2\log(1/a)+2\right]
\leq(\eta_Z/3)^4.
```

The last inequality uses $\log2<1$. This is the explicit normal-tail contribution.

From two independent $m$-bit integers $J_1,J_2$, form midpoint uniforms

```math
U_m=(J_1+1/2)2^{-m},\qquad V_m=(J_2+1/2)2^{-m},
```

and put $\widetilde U=\max(U_m,a)$. Couple $J_1=\lfloor2^mU\rfloor$ and $J_2=\lfloor2^mV\rfloor$. Clamping is nonexpansive, so $|\widetilde U-U_a|\leq2^{-m}$. For $u,v\geq a$,

```math
|\sqrt{-2\log u}-\sqrt{-2\log v}|
\leq\sqrt{2|u-v|/a}.
```

Consequently the error before numerical evaluation is bounded pointwise by

```math
2^{(1+s-m)/2}+8\sqrt{2s}\,2^{-m}.
```

Choose $m$ so each term is at most $\eta_Z/6$. For example it suffices to impose

```math
m\geq s+1+\left\lceil2\log_2(6/\eta_Z)\right\rceil,
\qquad
m\geq\left\lceil\log_2(96(s+1)/\eta_Z)\right\rceil.
```

Finally evaluate the log, square root, cosine and multiplication with total absolute output error at most $\eta_Z/3$. The triangle inequality in $L^4$ proves the claim. The output bound follows from $-2\log\widetilde U\leq2s\log2$ and the chosen arithmetic accuracy.

There is a direct polynomial-bit elementary-function implementation. One may use exact rational partial sums and outward rational error bounds, followed by one final dyadic rounding; their intermediate numerator and denominator lengths remain polynomial in the requested precision. In detail:

- For $\log x$, write $x=2^{-k}y$, $1\leq y<2$, and use $\log y=2\sum_{j\geq0}z^{2j+1}/(2j+1)$, where $z=(y-1)/(y+1)$ has magnitude at most $1/3$. The remainder after terms $j=0,\ldots,N-1$ is at most $2|z|^{2N+1}/[(2N+1)(1-z^2)]$. The same series with $z=1/3$ computes $\log2$. Since $k\leq s+1$, a number of terms linear in the requested bits and $\log(s+1)$ suffices.
- Compute square roots by dyadic bisection, using exact integer comparisons of squares. The number of iterations is linear in the output bits and the integer-part length.
- Compute $\pi=16\arctan(1/5)-4\arctan(1/239)$ by the alternating arctangent series, whose error is at most its first omitted term. Then compute cosine by its Taylor series on $[0,7]$; the remainder after degree $2N$ is at most $7^{2N+2}/(2N+2)!$. Increment $N$ until this rational upper bound is below the requested tolerance. Linear-in-precision term counts suffice.
- Schoolbook arithmetic with the additional guard bits makes each of these fixed-clock programs polynomial in $m,s$ and the output precision.

For completeness, square-root propagation does not assume a bounded derivative at zero. To obtain absolute root error $\gamma$, enclose its nonnegative argument to error at most $\gamma^2/16$, replace any negative numerical endpoint by zero, and bisect to error at most $\gamma/2$. The inequality $|\sqrt x-\sqrt y|\leq\sqrt{|x-y|}$ controls the rest. In the Box–Muller product, the radial factor is at most $\sqrt{2s}$, so cosine error $O(\eta_Z/\sqrt{s+1})$ and radial error $O(\eta_Z)$ suffice. These requests add only a constant factor or $O(\log(s+1))$ guard bits. All loop limits are obtained by checking the displayed rational bounds before the circuit is emitted.

No rejection loop, inverse-CDF oracle or probability table is involved. Reversible simulation of these finite Boolean programs is exact. Independent calls using disjoint uniform bits give an independent $q$-vector of approximate normals. Minkowski's inequality gives vector $L^p$ error at most $\sqrt q\,\eta_Z$ for $p=2,4$.

## 3. The exact finite-grid correction program

Let $I=[u,u+H]$ be a dyadic temporal support, $H=2^{-j}$. Let $h_l=2^{-l}$ and $n_l=2^l$. The signed detail multiplier is $H/2$; the endpoint coefficient uses multiplier $H$. We bound both by $H$.

For an ideal continuous endpoint $\tau$ and exact normals, define $Y_l$ as in [Dynkin details](dynkin-details.md): Euler is stopped at $\tau$ with a partial final step, and the generator uses $v(h_l\lfloor\tau/h_l\rfloor)$. Let $D_0=Y_0$ and $D_l=Y_l-Y_{l-1}$ for $l\geq1$. The ideal coupled variables telescope in expectation.

The finite program runs every fine-grid step, including the inactive suffix. At fine index $i$, it computes

```math
\Delta_i=\min\{h_l,\max(\widehat\tau-i h_l,0)\},
\qquad \Delta W_i=\sqrt{\Delta_i}\,\widehat Z_i.
```

It advances the fine state using the hardcoded value $v(i h_l)$. It accumulates two consecutive fine increments and, at the end of that pair, advances the coarse state with its hardcoded left-endpoint feedback and elapsed duration $\min\{2h_l,\max(\widehat\tau-2kh_l,0)\}$. A missing last fine increment is zero. The same initial state and the same fine increments are used for both paths.

The two frozen feedback values needed at the final generator evaluation can be accumulated by a one-hot conditional addition over the deterministic grid. Exactly one interval contains the endpoint. This costs $O(Kp n_l)$ bit operations and uses $O(Kp)$ working bits; it is not a quantum address lookup.

All comparisons, masks and loop indices are explicit. The program uses $O(n_l)$ normal blocks and coefficient evaluations. It stores only two current states, the coarse increment accumulator, the endpoint, the frozen feedback values, and local arithmetic workspace. This small classical streaming space will matter for the PRG construction.

## 4. Aligned endpoint bits and numerical stability

### 4.1 The endpoint causes no hidden discontinuity

Couple $\tau=u+HU$ to

```math
\widehat\tau=u+H\left(\lfloor2^{b_\tau}U\rfloor+1/2\right)2^{-b_\tau},
\qquad r=H2^{-b_\tau}.
```

Require $b_\tau\geq1$ and $j+b_\tau\geq L$, where $L$ is the largest path level used for this coefficient. Then $r$ divides every path-grid step. The sign switch at the midpoint of $I$ and every Euler-grid boundary are boundaries of endpoint bins. Except on a null event, $\tau$ and $\widehat\tau$ therefore have the same sign, the same active-step masks, and the same frozen feedback indices at every level.

Only the duration of the final partial step changes. With the same normal for that step, $|\sqrt a-\sqrt b|\leq\sqrt{|a-b|}$ shows, in the bounded core,

```math
\|X^{h_l}_{\tau}-X^{h_l}_{\widehat\tau}\|_2
\leq B_a r+\nu\sqrt r.
```

The coarse path obeys the same estimate: its partial increment is a sum of the shared fine increments, only one of which changes. The frozen mean argument is unchanged. This proves the endpoint approximation bound directly, without assuming that a discontinuous mask is Lipschitz.

### 4.2 A conservative perturbation bound

Let $n_*=2^L$. Suppose every elementary finite-grid update, after combining its arithmetic operations, has an absolute vector rounding error at most $u_{\rm upd}$. Let $u_a,u_F,u_v$ bound drift evaluation, generator evaluation and hardcoded-feedback errors. Let $e_0$ bound the initial-state coupling error. Then discrete Gronwall and the triangle inequality give, for either path at any used level,

```math
E_x\leq B_a r+\nu\sqrt r
+e^{L_x}\left[e_0+\nu_o\sqrt{q n_*}\,\eta_Z
                +u_a+L_m u_v+n_*u_{\rm upd}\right].
```

Indeed $\sum_i\sqrt{\Delta_i}\leq\sqrt{n_l}$, the sum of the elapsed time steps is at most one, and $n_l\leq n_*$. This intentionally conservative estimate does not rely on cancellation between rounding errors. It applies separately to the fine and coarse states generated from the same approximate normal blocks.

For a correction, including the base level with an unnecessary factor of two, the resulting payoff error satisfies

```math
\|\widehat D_l-D_l\|_2
\leq 2H\left(L_{Fx}E_x+L_{Fm}u_v+u_F\right)+u_{\rm out},
```

where $u_{\rm out}$ accounts for the final multiplication, subtraction and output rounding. These inequalities are an explicit finite-bit error certificate. A compiler chooses dyadic tolerances satisfying them; it need not estimate their errors experimentally.

For example define a supplied upper constant $C_{\rm fp}\geq1$ larger than every coefficient multiplying $\sqrt r$, $e_0$, $\sqrt{n_*}\eta_Z$, $u_a$, $u_v$, $n_*u_{\rm upd}$, $u_F$ and $u_{\rm out}/H$ in the last two displays; use $r\leq\sqrt r$ for $r\leq1$. For a requested correction perturbation $\rho$, set $z=\min\{1/4,\rho/(16HC_{\rm fp})\}$ and enforce

```math
r\leq z^2,\qquad e_0\leq z,\qquad
\eta_Z\leq z/\sqrt{n_*},\qquad
u_a,u_v,u_F\leq z,\qquad
u_{\rm upd}\leq z/n_*,\qquad u_{\rm out}\leq Hz.
```

Increasing $C_{\rm fp}$ once if necessary makes $\|\widehat D_l-D_l\|_2\leq\rho$. The primitive arithmetic is evaluated with enough guard bits that its composite update error obeys $u_{\rm upd}$. One can obtain this by interval error propagation through the explicit arithmetic circuit; all operations here have polynomially bounded condition factors on the finite range. Square roots at zero use the bisection absolute-error bound rather than a derivative bound.

### 4.3 Range and word lengths

Finite normals satisfy the deterministic bound in Lemma 1. Under the bounded drift, all seed strings obey the finite-path range estimate

```math
|\widehat X_i|
\leq R_X:=|x_0|+B_a+2+
3\nu_o\sqrt{q n_*}\sqrt{s+1},
```

after allowing the allocated accumulated numerical error. For a Gaussian initial law, add its analogous deterministic finite-normal range bound. A register with integer part longer than $\log_2(2R_X)$ cannot overflow. Since $R_X$ is polynomial in $n_*$ and the model parameters, its integer-part length is logarithmic in $n_*$ for a fixed model.

With the preceding tolerance choices and $n_*=O(1+HB_k/\delta)$, or fixed-model $O(1+H/\delta)$, a sufficient common word length is

```math
p=O\!\left(
1+\log n_*+\log(2+H/\rho)+\log(1+C_{\rm fp})+
\log(1+d+q+K+|P|)
\right),
```

plus the explicitly supplied parameter bit lengths and any model-dependent integer range. In particular $p=O(\log(1/\epsilon))$ for a fixed admissible model. This is proved using the elementary Gaussian construction and the above update bounds, rather than postulating an accuracy-$\epsilon$ sampler at unit cost.

## 5. Nonmatching discrete marginals and the full coefficient error budget

Quantized fine Gaussian increments generally do not sum to the same distribution as a separately quantized coarse normal. Thus the implemented corrections do **not** have to telescope exactly. The comparison is to the ideal coupled hierarchy.

For a coefficient tolerance $\delta$, choose $L$ so that the analytic Euler bias is at most $\delta/4$, and choose

```math
\rho=\frac{\delta}{16(L+1)}.
```

Use Section 4 to ensure $\|\widehat D_l-D_l\|_2\leq\rho$ at every level. Then

```math
\left|\sum_{l=0}^L\mathbb E\widehat D_l-\mathbb E Y_L\right|
\leq(L+1)\rho=\delta/16.
```

If the ideal second-moment bounds are $\|D_l\|_2\leq S Hh_l$ for $l\geq1$ and $\|D_0\|_2\leq B_FH$, then

```math
\|\widehat D_l\|_2\leq S Hh_l+\rho.
```

These are directly usable bounds, even when the analytic bias constant is zero. Alternatively choose $L$ using a positive supplied majorant of all analytic constants; then $\rho=O(Hh_L)$ and the bounds retain the form $O(Hh_l)$ with an enlarged constant. No small pointwise range proportional to $Hh_l$ is inferred.

Estimate each finite-bit mean to additive error $e_l=\delta/[4(L+1)]$. The total statistical allowance is $\delta/4$. Allow $\delta/16$ for rounding the returned coefficient. The optional PRG adds at most another $\delta/16$, as proved in [Memory and recomputation](memory-and-recomputation.md). Together with the ideal Euler bias and finite-bit bias, these contributions total at most $11\delta/16<\delta$. The spare margin covers a rounded supplied constant or an initial-mean error if the outer algorithm elects to allocate it here.

This bound is conditional on the previously computed feedback curve satisfying its stated regularity invariant. The guarantees and operation counts are uniform over all such curves. A conditional failure bound at every call, followed by a union bound, therefore applies to adaptive classical Picard iterations. The circuit always runs its predetermined clock and has its predetermined register range; an earlier estimation failure does not cause an unbounded run.

## 6. An explicit second-moment quantum mean estimator

### Lemma 2: Mean estimation interface

Suppose a finite deterministic sampler $A$ produces a signed dyadic variable $X$ from uniform random bits, $\mathbb E X^2\leq\sigma^2$, and $\sigma$ is a supplied dyadic power-of-two upper bound. With a coherent reversible implementation of $A$ and its inverse, its mean can be estimated to additive error $e>0$ and failure probability at most $\alpha$ with

```math
O\!\left[
\left(1+\frac{\sigma}{e}\right)
\log^2\!\left(2+\frac{\sigma}{e}\right)
\log\!\left(\frac{2+\log(2+\sigma/e)}{\alpha}\right)
\right]
```

calls to the sample circuit and its inverse. The extra arithmetic and quantum registers have size polynomial in the input word length and $\log(2+\sigma/e)$. The bound uses a second moment, so no centering procedure is necessary.

**Proof and circuit.** If $e\geq\sigma$, output zero, using $|\mathbb E X|\leq\sigma$. Otherwise set $a=\sigma/e>1$ and $J=\lceil\log_2(8a)\rceil$. Discard the part with $|X|>2^J\sigma$. Its mean contribution has magnitude at most

```math
\frac{\mathbb E X^2}{2^J\sigma}\leq e/8.
```

Separate signs and divide magnitudes into bin zero $|X|\leq\sigma$ and bins $2^{j-1}\sigma<|X|\leq2^j\sigma$, $1\leq j\leq J$. For each sign and bin define the nonnegative dyadic variable

```math
f_j^{\pm}(X)=\frac{|X|}{2^j\sigma}\,
\mathbf 1\{X\text{ has the chosen sign and lies in bin }j\}.
```

For $j=0$, use denominator $\sigma$. Every $f_j^\pm$ lies in $[0,1]$. Its expectation $p_j^\pm$ obeys $p_0^\pm\leq1$ and $p_j^\pm\leq2^{1-2j}$ for $j\geq1$, by the second-moment bound.

Amplitude encoding is exact. If the computed dyadic $f$ is $k/2^b$, append a uniform $b$-bit integer $R$ using Hadamards and mark the predicate $R<k$. Store $k$ with one additional bit so $f=1$ marks every possible $R$. Its success probability is exactly $f$. All required operations are Boolean comparisons, bit shifts, Hadamards and reversible arithmetic; an arbitrary controlled rotation is not assumed. Since $\sigma$ is a power of two, bin normalization is a bit shift with a finite exact denominator.

Apply the amplitude-estimation algorithm of BHMT, Theorem 12 [BHMT02], with a power-of-two query parameter

```math
M\geq128\pi a(J+1).
```

For each probability, one run obeys

```math
|\widetilde p-p|\leq
\frac{2\pi\sqrt{p(1-p)}}{M}+\frac{\pi^2}{M^2}
```

with probability at least $8/\pi^2>3/4$. Multiplying by $2^j\sigma$ and summing over both signs and all bins bounds the first terms by $4\pi\sqrt2\sigma(J+1)/M$ and the second terms by $4\pi^2\sigma2^J/M^2$. With the stated $M$, these are less than $e/8$ in total. Adding the discarded tail still leaves error less than $e/4$.

Repeat each probability estimate an odd number at least $8\log(2(J+1)/\alpha)$ of times and take its median. Hoeffding's inequality and a union bound make all bin estimates valid with probability at least $1-\alpha$. The amplitude-estimation measurement is an integer; its reported probability $\sin^2(\pi y/M)$ is computed classically with the elementary arithmetic above. Accuracy $e/[16(J+1)2^J\sigma]$ in each such probability makes the weighted contribution of all these evaluations at most $e/8$. Round the final signed sum to error at most another $e/8$. This proves the claimed accuracy with room in the constants. The total calls are $O(M(J+1)\log((J+1)/\alpha))$, which is the displayed bound.

This is a transparent dyadic-bin version of Montanaro's second-moment construction, Algorithm 2 and Lemma 4 [M15]. The slightly larger logarithmic factor is acceptable here. The general quantum MLMC result is Theorem 2 of An et al. [ALLMSW21]; the path-level cost calculation below does not need to assume exact finite-bit telescoping.

## 7. Reversible gates, synthesis, and success probability

### 7.1 Direct-history sample circuit

A finite Boolean circuit can be evaluated reversibly by writing each gate result into a clean bit, copying the desired output, and reversing the computation. AND uses a Toffoli gate; XOR and NOT use CNOT and $X$. This gives an exact reversible circuit with gate count linear in the Boolean circuit size. Local arithmetic can be uncomputed after each path step while retaining the old states and random bits required for reversal.

For the streaming sample program with $n_l$ steps, write

```math
C_s(p)=\operatorname{poly}(d,q,K,p,|P|)
+O\bigl(C_{\rm eval}(P,p,\log R_X)\bigr)
```

for a fully charged per-step Boolean bound, including the Gaussian routine, coefficient and generator computations, controlled partial steps, and arithmetic guards. Then a sample and its inverse cost $O(C_s(p)(n_l+1))$ logical gates. A conservative workspace bound is $O(C_s(p)(n_l+1))$ bits; with local arithmetic uncomputed, one can instead keep $O((d+q+K)p(n_l+1)+W_{\rm eval})$ bits. The small-space alternative is not silently included in this direct bound.

Amplitude estimation uses reflections about the sample-prepared state and about the marked predicate. Reflection about the all-zero computational state is implementable with linearly many elementary reversible gates and ancillas. Controlled Grover iterations can be written using unconditional sample preparation/unpreparation and controlled reflections, so no unavailable black-box controlled sampler is required. Toffoli gates and these Boolean controls have exact Clifford+$T$ decompositions.

### 7.2 The only approximation in the quantum gates

The sampling arithmetic and the comparator amplitude encoding are exact finite Boolean unitaries. The Fourier transforms used by the selected amplitude-estimation algorithm contain rotations. Compile these fixed-dimensional gates into a finite inverse-closed universal gate set, for example Clifford+$T$, using a constructive Solovay–Kitaev algorithm [DN06]. This adds a polynomial in the logarithm of the inverse synthesis tolerance to both gate count and classical synthesis work.

Let $N_{\rm rot}$ be a deterministic upper bound on all rotation occurrences over every coefficient, level, repetition and Picard iteration. Let the sum of the ideal estimation failure allowances be at most $\alpha_{\rm tot}/2$. Synthesize each rotation to operator-norm error at most

```math
\epsilon_{\rm gate}=\frac{\alpha_{\rm tot}}{16\max(1,N_{\rm rot})}.
```

The telescoping product inequality for unitaries bounds the total state error by the sum of the gate errors; the induced measurement total variation error is at most twice that sum. The same hybrid argument applies to adaptive classical control, conditional on the prior outcomes. Thus synthesis changes the final failure probability by less than $\alpha_{\rm tot}/2$. It does not require an unaccounted additive numerical error in the coefficient.

The worst-case count is polynomial in $1/\epsilon$ for a fixed model, so $\log(1/\epsilon_{\rm gate})=O(\log(1/\epsilon)+\log(1/\alpha_{\rm tot}))$. Gate synthesis preserves the precision exponent. Approximate inverse circuits are obtained by exactly reversing the synthesized gate list, avoiding an independent inconsistency between a unitary and its inverse.

## 8. End-to-end resource theorem in the bounded core

### Theorem 1: Finite-bit direct-history implementation

Assume the analytic conditions and feedback reconstruction theorem of the preceding chapters, together with the evaluator contract in Section 1. For any requested additive error $0<\epsilon<1$ and failure probability $0<\alpha_{\rm tot}<1/2$, the algorithm has a uniform finite-circuit implementation over a fixed universal gate set. Its classical output satisfies the analytic error guarantee with probability at least $1-\alpha_{\rm tot}$, using

```math
G_{\rm direct}=
\widetilde O\!\left(
\frac{C_{\rm model}\,
C_s(p)}{\epsilon}
\right).
```

Here $p$ is the word length in Section 4 with the logarithmic confidence and synthesis precision included. The factor $C_{\rm model}$ includes the explicit stability, derivative, diffusion and norm factors, the component count and the Picard/reconstruction allocation. It is not asserted to be polynomial in every physical parameter. The tilde hides only logarithms and the polynomial-bit elementary arithmetic and synthesis factors already described.

Let $C_{\rm init}(p)$ be the charged certified initial-mean evaluation cost, and let $C_{\rm build}(p)$ majorize $C_s(p)$, the circuit-generation cost $C_{\rm gen}(P,p,\log R_X)$, and the per-step constant-generation arithmetic. A sufficient classical input, preparation, compilation and output-work bound is

```math
W_{\rm classical}
\leq\widetilde O\!\left(
|P|+C_{\rm init}(p)+
\frac{C_{\rm model}C_{\rm build}(p)}{\epsilon}
\right).
```

The simplified near-linear end-to-end precision exponent applies when all these input-program costs are polynomial in the charged word length at fixed model parameters. The full model input and classical feedback/output storage are separately charged; neither a slow initial-mean evaluator nor slow circuit generation is hidden inside a unit-cost query.

**Proof.** For one support of width $H$ and tolerance $\delta$, the analytic cutoff gives $n_*=O(1+H/\delta)$ with its constants displayed in the analytic theorem. Finite-bit errors and the allocated statistical estimates satisfy Section 5. Lemma 2, with the finite-bit second-moment bound, gives a per-level cost bounded up to logarithms by

```math
C_s(p)\,h_l^{-1}
\left[1+\frac{(SHh_l+\rho)(L+1)}{\delta}\right].
```

The geometric sum of mandatory sample costs is $O(C_s(p)n_*)$. The $SHh_l$ terms sum to $O(C_s(p)SH(L+1)^2/\delta)$, and the $\rho$ terms are absorbed in the mandatory costs because $\rho(L+1)/\delta=1/16$. Hence one coefficient costs $\widetilde O(C_s(p)(1+H/\delta))$.

There are $O(1/H)$ coefficients at a temporal scale and $O(\log(1/\epsilon))$ scales. With the allocated coefficient accuracy $\delta=\Theta(\epsilon/\operatorname{polylog}(1/\epsilon))$, their costs sum to the claimed exponent. The Lipschitz-only representation has $O(K/\epsilon)$ classical coefficients, so reading and writing it is within the bound. The terminal scalar stage has the same path-level construction, with $f$ replacing the generator payoff.

With ordinary classical addressable memory, generating the deterministic feedback constants costs $O(K\operatorname{poly}(p)\log M)$ per grid point if the previous curve has $M$ nodes, counting the transferred word bits. Alternatively precompute all path-grid tables once per Picard iteration by sorted merges in $O(K\operatorname{poly}(p)[M(L_{\max}+1)+n_{\max}])$ work and stream the relevant table. [Memory and recomputation](memory-and-recomputation.md#4-uniform-compilation-and-classical-feedback-storage) proves that even the checkpoint replay can read these tables sequentially within its gate bound. Across path levels of one coefficient the number of grid points is $O(n_*)$. Thus classical schedule generation, Boolean circuit emission, and gate synthesis also have the same precision exponent. Reading the initial model description and emitting the requested output are included. This completes the construction.

## 9. Linear-confining drift and unbounded generator values

The polymer application has $a(x,m)=-Ax+c(x,m)$ with $c$ bounded and globally Lipschitz. Its generator satisfies a weighted, rather than global, Lipschitz bound. Assume the preceding analytic chapter supplies uniform fourth moments and an $L^4$ Euler strong-order-one estimate, and that

```math
|F(x,m)-F(y,n)|
\leq C_F(1+|x|+|y|)|x-y|+L_{Fm}|m-n|.
```

Use Lemma 1 in $L^4$ and allocate an $L^4$ initial-state error. Let $M_4^\star$ denote the supplied uniform fourth-moment norm of the **ideal Euler states**, as in the preceding analytic chapter; it may be larger than the exact diffusion's bound. With $a_0=\|A\|$ and $B_0=\|c\|_\infty$, the explicit version of Section 4's bound is

```math
E_{x,4}\leq(a_0M_4^\star+B_0)r+3^{1/4}\nu\sqrt r
+e^{L_x}\left[e_{0,4}+\nu_o\sqrt{q n_*}\,\eta_Z
                +u_a+L_m u_v+n_*u_{\rm upd}\right].
```

The factor $3^{1/4}$ follows from the fourth moment of a Gaussian vector. The same normal-coupling and discrete Gronwall argument controls the remaining terms. The implemented state has fourth-moment norm at most $M_4^\star+E_{x,4}$. Hölder's inequality yields

```math
\|F(\widehat X,\widehat m)-F(X,m)\|_2
\leq C_F(1+2M_4^\star+E_{x,4})E_{x,4}+L_{Fm}u_v.
```

Choosing the state tolerances smaller by this supplied factor gives exactly the correction perturbation requirement in Section 5. No fourth spatial derivative of $g$ is introduced by this finite-bit step.

For all seed strings, a conservative deterministic state bound is

```math
R_X\leq e^{\|A\|}
\left[R_{X_0}+\|c\|_\infty+2+
3\nu_o\sqrt{q n_*}\sqrt{s+1}\right].
```

Consequently the finite-bit generator range can grow like $H\sqrt{n_*}\operatorname{polylog}(n_*)$, with additional model factors. It is not generally polylogarithmic. The integer registers still need only $O(\log n_*)$ bits for a fixed model. This larger range is explicitly used in the PRG error budget in the next chapter. No clipping of the physical state is assumed.

The evaluator contract must hold up to this range, and its dependence on $\log R_X$ is charged. Stiffness, matrix norm and initial covariance factors are part of $C_{\rm model}$. A continuum refinement that increases them requires a separate joint-parameter bound.

## 10. Concrete arithmetic for the finite-bead Fourier model

For $L$ beads and $K=2R_f$ mean features, store the ordinary parameter list $(\kappa,D,T,q_a,v_a,c_0,\Sigma_C)$. A square root or explicit factor of the fixed-dimensional center covariance can be supplied or computed to certified precision. The initial sampler in [Polymer dynamics](polymer-dynamics.md) is a linear transformation of $O(L)$ independent normals and uses $O(L)$ arithmetic. Its $L^4$ coupling error is bounded by the norm of that explicit linear map times the normal-vector error. This determines the normal tolerance; it is not a separate state-preparation assumption.

The open-chain spring matrix is applied by nearest-neighbor differences. At every bead and mode, compute a three-dimensional dot product, its sine and cosine, and the finite Fourier force. The explicit generator formulas use the same drift and trigonometric values. A fixed finite-model evaluator therefore uses $O(LK)$ elementary arithmetic operations, each with polynomial bit cost, rather than a dense Hessian routine.

For potentially large phase arguments $q_a\cdot r_j$, first reduce modulo $2\pi$. Computing $\pi$ to $p+O(\log(2+|q_a|R_X))$ bits makes the reduced argument accurate to $O(2^{-p})$. Integer division and the bounded-argument sine/cosine series are polynomial in that bit length. The drift and generator sums receive additional $O(\log(LK))$ guard bits and guards for the explicit coefficient norms. Gaussian characteristic-function formulas for the initial means use the same trigonometric arithmetic and exponential evaluation. Exponentials on a nonpositive argument can be truncated to zero once a certified bound is below the target error; on the remaining bounded range, range reduction and a Taylor series give polynomial-bit work.

Thus this application satisfies the input contract with explicit polynomial-bit $C_{\rm eval}$, $W_{\rm eval}$ and $C_{\rm gen}$. The arithmetic count is not a physical fault-tolerant gate estimate, and it does not establish a crossover against QMC or other structured classical methods.

## 11. Scope and provenance

The new mean-field argument uses a randomized-time generator payoff and its mixed temporal/path estimate. This chapter supplies a constructive precision and resource closure under explicit program contracts. Gaussian transformations, reversible Boolean evaluation, amplitude estimation, second-moment binning, quantum MLMC and gate synthesis are established ingredients.

- **[BHMT02]** G. Brassard, P. Høyer, M. Mosca and A. Tapp, *Quantum Amplitude Amplification and Estimation*. Theorem 12 is the amplitude-error and success-probability bound used in Lemma 2. [Primary preprint](https://arxiv.org/abs/quant-ph/0005055).
- **[BM58]** G. E. P. Box and M. E. Muller, *A Note on the Generation of Random Normal Deviates*, The Annals of Mathematical Statistics 29 (1958), 610–611. [Publisher DOI](https://doi.org/10.1214/aoms/1177706645). The finite cutoff, quantization and coupling budget used here are proved in Lemma 1.
- **[M15]** A. Montanaro, *Quantum speedup of Monte Carlo methods*. Algorithm 2 and Lemma 4 treat bounded second moments; Algorithm 3 and Theorem 5 treat a variance bound by centering. [Primary preprint, version 3](https://arxiv.org/abs/1504.06987v3); [published paper](https://doi.org/10.1098/rspa.2015.0301).
- **[ALLMSW21]** D. An, N. Linden, J.-P. Liu, A. Montanaro, C. Shao and J. Wang, *Quantum-accelerated multilevel Monte Carlo methods for stochastic differential equations in mathematical finance*. Theorem 2 is the general quantum MLMC theorem; Theorem 3 is the SDE specialization. [Primary paper](https://quantum-journal.org/papers/q-2021-06-24-481/); [preprint, version 2](https://arxiv.org/abs/2012.06283v2).
- **[DN06]** C. M. Dawson and M. A. Nielsen, *The Solovay–Kitaev algorithm*. Theorem 1 and the constructive analysis give polynomial-logarithmic circuit length and classical synthesis time in inverse operator-norm accuracy. [Primary preprint](https://arxiv.org/abs/quant-ph/0505030).

The repository's budget checker verifies the numerical inequalities and counting recurrences used here. It does not constitute a quantum execution or a verification of all input-model promises.
