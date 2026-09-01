# Theory notes — what the measured structure is, and what it is not

**Scope.** These notes turn the batch-1 measurements into statements one can try to prove. Each item
is tagged [PROVED] / [VERIFIED …] / [CONJ] / [OPEN] exactly as in `paper/claims_ledger.md`; a proof
sketch is given where one exists, and where the intended description was *refuted* that is recorded
as prominently as a success. Scripts: `experiments/E03_syzygy_and_degree_fall/fall_generators.py`
(`--sbox`, `--kernel`, `--orbits`, `--all`), data:
`results/E03_theory_fall_generators.csv`, `results/E03_theory_orbit_lemma.csv`,
`results/E06_theory_fall_generators.csv`.

---

## 0. The tool: trace coordinates and Frobenius orbits

Write a field polynomial $F \in \mathbb{F}_{2^n}[x,y]$ in the *field-valued Boolean* form
$F=\sum_m c_m\,m$, where $m$ runs over squarefree monomials in the $2n$ bits and $c_m\in
\mathbb{F}_{2^n}$ (`src/aim3/fieldpoly.py`; a cross-check against `legacy/aim3_macaulay.py` is in
`tests/test_core.py`). The $\mathbb{F}_2$-span of $F$'s $n$ Boolean coordinates is then

$$V(F)=\{\operatorname{Tr}(\gamma F):\gamma\in\mathbb{F}_{2^n}\},\qquad \dim_{\mathbb{F}_2}V(F)\le n .$$

Two facts do all the work:

1. **Boolean degree = exponent weight.** $x^a$ is the product of the $\mathbb{F}_2$-linear forms
   $x^{2^k}$, $k\in\operatorname{supp}(a)$, so as a Boolean function it has degree $\operatorname{wt}(a)$.
   In the field-valued representation this is automatic, because monomial keys are *sets* and
   $m\cup m=m$ implements $x_j^2=x_j$.
2. **Frobenius collapse.** $\operatorname{Tr}(\gamma z^{2^s})=\operatorname{Tr}(\gamma^{2^{-s}}z)$, so
   $V(x^{a}y^{b})=V(x^{2^s a}y^{2^s b})$ (exponents mod $2^n-1$). Two trace terms can cancel **only**
   when their exponent pairs lie in one Frobenius orbit.

**Why every Macaulay row is of this shape.** A Boolean bit $x_c$ is $\operatorname{Tr}(\delta_c x)$
for a dual-basis element $\delta_c$, so a squarefree bit-monomial of degree $\le 2$ is a
Boolean-degree-$\le2$ *function* of $(x,y)$, and those are exactly the $\mathbb{F}_{2^n}$-span of the
monomials $x^ay^b$ with $\operatorname{wt}(a)+\operatorname{wt}(b)\le2$ (uniqueness of the reduced
polynomial representation of a function $\mathbb{F}_{2^n}\to\mathbb{F}_{2^n}$). Since
$\mu\cdot\operatorname{Tr}(z)=\operatorname{Tr}(\mu z)$ for $\mu\in\{0,1\}$, every row of
$\mathrm{Mac}_4(\{P_1,P_2,P_3\})$ is $\operatorname{Tr}(\gamma\,M\,P_k)$ for such an $M$. **The 18
families $\{1,x^{2^i},y^{2^j},x^{2^i}y^{2^j},x^{2^i+2^j},y^{2^i+2^j}\}\times\{P_1,P_2,P_3\}$
therefore span the whole degree-4 row space** — the question is never "is it in the module", but
"which families does a fall need".

---

## 1. The $3n$ structural syzygies at $D=3$ [PROVED]

$P_2=xP_1$ and $P_3=yP_1$ as field polynomials. Hence each coordinate of $P_2$ is an
$\mathbb{F}_2$-combination of the rows $x_j\cdot(\text{coordinate of }P_1)$, giving $n$ dependencies
among the $\mathrm{Mac}_3$ rows; the same for $P_3$; and $y\cdot P_2=xyP_1=x\cdot P_3$ gives $n$ more.
That is $3n$ of the measured $8n$. The remaining $5n$ are [VERIFIED n=11..16, every admissible e]
(36 runs, `results/E03_sbox_degree_fall.csv`) and **[OPEN]** as a theorem; the mechanism is expected
to be Frobenius collapse, but see §2 for why the analogous $D=4$ question turned out to have no small
structured answer.

---

## 2. The $D=4$ degree falls: the intended description is **refuted**

### 2.1 Where a degree-4 part can cancel

For $i,j\ge1$,

$$x^{2^i}y^{2^j}P_1=\underbrace{x^{2^i+1}y^{2^j+1}}_{\text{Boolean degree }4}
  +\underbrace{x^{2^i+2^e}y^{2^j}}_{\le3}+\underbrace{x^{2^i}y^{2^j}}_{2},$$

and $\operatorname{Tr}(\gamma x^{2^i+1}y^{2^j+1})$ is **homogeneous** of degree 4 (every Frobenius
image of the exponent pair has the same weight). So a fall inside this family is exactly a solution
of $\sum_{i,j}\operatorname{Tr}(\gamma_{ij}x^{2^i+1}y^{2^j+1})=0$, i.e. a pure Frobenius coincidence.
Generators with $i=0$ or $j=0$, and $x^{2^k}P_2$, $y^{2^k}P_3$, have no degree-4 part at all
($\operatorname{wt}(2^0+1)=\operatorname{wt}(2)=1$, and $x^{2^k}P_2=x^{2^k+2}y+x^{2^k+2^e+1}+x^{2^k+1}$
has degree $\le3$), so they contribute their whole $n$ dimensions to the kernel for free.

### 2.2 The orbit lemma [PROVED, checked n=5..20]

> **Lemma.** For $1\le i,j\le n-1$, the pairs $(2^i+1,2^j+1)$ and $(2^{i'}+1,2^{j'}+1)$ lie in one
> Frobenius orbit iff $(i',j')=(i,j)$, or $i=j$ and $(i',j')=(n-i,n-i)$.
>
> *Proof.* $2^s(2^i+1)$ reduces mod $2^n-1$ to the element with support $\{i+s,\,s\}\pmod n$, which
> has the form $\{i',0\}$ only if $s\equiv0$ (returning $(i,j)$) or $i+s\equiv0$, i.e. $s=n-i$. In the
> second case the other coordinate $2^s(2^j+1)$ has support $\{j-i,\,n-i\}$, which contains $0$ only
> if $j\equiv i$. $\square$

So there are exactly $c=\lfloor(n-1)/2\rfloor$ coincidences, one per diagonal pair $\{i,n-i\}$.
`fall_generators.py --orbits` checks the lemma and the counts for $n=5..20$
(`results/E03_theory_orbit_lemma.csv`).

### 2.3 Closed form for the module, and its exact match with the measurement

With $n^2+2n$ generators and $n$ trace coefficients each:

$$\text{top-degree rank}=n\bigl((n-1)^2-c\bigr),\qquad
  \dim\ker=n\bigl(n^2+2n-(n-1)^2+c\bigr),\qquad c=\lfloor(n-1)/2\rfloor .$$

At $n=11$: predicted rank $1045$, kernel $528$ — and `fall_generators.py --sbox 11 3` measures
exactly $1045$ and $528$. The accounting is complete: $473=11\cdot43$ kernel dimensions come from the
generators with no degree-4 part and $55=11\cdot5$ from the diagonal Frobenius coincidences.

### 2.4 …and none of it is new [VERIFIED (n,e) = (11,3), (13,3)]

**Every one of those $528$ (resp. $741$ at $n=13$) kernel elements already lies in
$\mathrm{Mac}_3$.** The module
$\{x^{2^i}y^{2^j}P_1\}\cup\{x^{2^k}P_2\}\cup\{y^{2^k}P_3\}$ — the description proposed for
Theorem T1.2 — explains **$0$ of the $2n$** new cubics, at all four measured pairs
$(n,e)=(11,3),(13,3),(11,4),(13,5)$. Furthermore, **no single one of the 18 multiplier families
produces a new cubic on its own** (all 18 measured, all $0$, at $(11,3)$ and $(11,4)$).

This is a negative result, and it is the main theoretical output of this batch: the first-fall cubics
are *not* generated by a one- or two-family field identity of the syzygy kind. They are genuine
cross-family cancellations.

### 2.5 A structural obstruction to the "two field cubics" picture

$\mathrm{falls}_4=2n$ is a multiple of $n$, which invites the guess that the quotient is
$V(g_1)\oplus V(g_2)$ for two field cubics — $n$ coordinates each. **That cannot be the general
shape**: the classes $\{-3,(-3)^{-1}\}$ and $\{n/2-1\}$ give $\mathrm{falls}_4=2n+1$
(`results/E03_summary.md`), and $2n+1$ is not a multiple of $n$ for $n>1$. Any correct description
must therefore produce a lone extra element in those classes — a $\dim$-$1$ phenomenon on top of two
$n$-blocks, not a third block.

### 2.6 Leave-one-out: exactly two families are indispensable

Since the 18 families together span the whole row space, the meaningful question is which of them
cannot be dropped. Removing one family at a time and recomputing the falls modulo $\mathrm{Mac}_3$:

| dropped family | falls left (of $2n=22$) |
|---|---|
| any of the other 16 | 22 |
| $x^{2^i}y^{2^j}\cdot P_2$ | **11 $=n$** |
| $y^{2^i+2^j}\cdot P_3$ | **11 $=n$** |

[VERIFIED $(n,e)=(11,3)$ **and** $(11,4)$ — all 18 families measured in each, and the
two indispensable ones are **the same pair** for both exponents]

So $\mathrm{falls}_4=2n$ splits as $n+n$: one block needs the multipliers $x^{2^i}y^{2^j}$ on $P_2$,
the other needs $y^{2^i+2^j}$ on $P_3$. Note both are **weight-2** multipliers on $P_2/P_3$, whereas
the refuted ansatz used the **weight-1** ones $x^{2^k}P_2$, $y^{2^k}P_3$.

### 2.7 Why exactly those two — the degree-4 count

Expand the three candidate shapes as multiples of $P_1$ and count the degree-4 field monomials:

| generator | $=$ | degree-4 monomials | Boolean degree |
|---|---|---|---|
| $x^{2^k}P_2 = x^{2^k+1}P_1$ | $x^{2^k+2}y+x^{2^k+1+2^e}+x^{2^k+1}$ | **0** | $\le3$ |
| $y^{2^k}P_3$ | (symmetrically) | **0** | $\le3$ |
| $x^{2^i}y^{2^j}P_1$ | $x^{2^i+1}y^{2^j+1}+x^{2^i+2^e}y^{2^j}+x^{2^i}y^{2^j}$ | **1** | 4 |
| $x^{2^i}y^{2^j}P_2=x^{2^i+1}y^{2^j}P_1$ | $x^{2^i+2}y^{2^j+1}+\mathbf{x^{2^i+1+2^e}y^{2^j}}+x^{2^i+1}y^{2^j}$ | **2** | 4 |
| $y^{2^i+2^j}P_3=y^{2^i+2^j+1}P_1$ | $xy^{2^i+2^j+2}+\mathbf{x^{2^e}y^{2^i+2^j+1}}+y^{2^i+2^j+1}$ | **2** | 4 |

(verified symbolically with `src/aim3/fieldpoly.py`). This explains all three observations at once:

* $x^{2^k}P_2$ and $y^{2^k}P_3$ have **no** degree-4 part at all — they are already inside
  $\mathrm{Mac}_3$ and can contribute nothing. That is one half of why the ansatz fails.
* $x^{2^i}y^{2^j}P_1$ has a **single** degree-4 monomial, so its only cancellations are the
  $\lfloor(n-1)/2\rfloor$ diagonal Frobenius coincidences of §2.2 — and those, measured, all land in
  $\mathrm{Mac}_3$ too. That is the other half.
* The two indispensable families have **two** degree-4 monomials, and in each case the second one
  **carries $2^e$** (bold above). Cancellation can therefore happen between the two monomials of one
  generator, and the condition involves the exponent. **This is the mechanism by which
  $\mathrm{falls}_4$ depends on the multiplicative class of $e$ modulo $n$** — a single-monomial
  family could not see $e$ at all.

[CONJ] A proof of Theorem T1.2 should start here: solve
$\operatorname{Tr}(\gamma\,x^{2^i+2}y^{2^j+1})+\operatorname{Tr}(\gamma'\,x^{2^i+1+2^e}y^{2^j})=0$
for the exponent coincidences, and the same for the $P_3$ family; the class table of
`results/E03_summary.md` should drop out, including the sign asymmetry of $3e\equiv\pm1$ and the lone
extra element of the $2n+1$ classes (which §2.5 says cannot be an $n$-block).

### 2.8 An explicit basis, and its kernel vectors

`fall_generators.py --basis n e` produces the deliverable directly: every row of $\mathrm{Mac}_4$ is
labelled by its (equation, multiplier) pair; the rows whose multiplier has weight $\le1$ are already
of degree $\le3$ and span $\mathrm{Mac}_3$ exactly, so they are eliminated **first**; anything that
then acquires a pivot in the low region is new modulo $\mathrm{Mac}_3$, and its accumulated
combination *is* the kernel vector. The count comes out at exactly $2n$ every time, and the basis
splits cleanly in two ($n$ + $n$), by whether the kernel vector uses any weight-1 multiplier at all:

| $(n,e)$ | $2n$ | block A (uses wt-1 multipliers) | block B (uses none) | $P_3$ rows in A | $P_3$ rows in B | kernel-vector size |
|---|---:|---:|---:|---:|---:|---:|
| (11,3) | 22 | 11 | 11 | 46–76 | 634–699 | 1645–2677 |
| (11,4) | 22 | 11 | 11 | 55–73 | 643–701 | 1656–2668 |
| (13,3) | 26 | 13 | 13 | 53–97 | 1234–1283 | 2774–4734 |
| (13,5) | 26 | 13 | 13 | 83–107 | 1182–1286 | 2807–4729 |

(`results/E03_theory_fall_basis.csv`, one row per basis element.) The $n+n$ split is the same one the
leave-one-out finds, seen from the other side: block B is the $P_3$-heavy half.

**Two caveats, and one more refutation.** The basis is echelon-derived, so an individual element is
not canonical — only the counts and the coarse A/B split are. And the kernel vectors are **dense**:
1645–4734 of the 8382 (resp. 13728) labelled rows take part in each one. There is therefore no sparse
identity of the syzygy kind ($3n$ of the $D=3$ syzygies each use $O(1)$ rows) waiting to be found;
that is independent confirmation of §2.4, from the kernel side rather than the module side.

---

## 3. The kernel model's $6n$ new quadratics at $D=3$ — fully explained [VERIFIED n=7,9]

The kernel model (E06) has, per S-box,
$$B_i=a_ix_i^{2^{e_i}}+a_i^{2^{e_i}}x_i+a_i,\qquad K_i=t_ia_i+a_i^{2^{e_i}},$$
all bilinear, $6n$ equations in $6n$ variables. The $D=3$ pass adds exactly $6n=\#\mathrm{eqs}$ new
quadratics (§1 of the batch-2 wording correction). They are

$$\boxed{\;a_i\cdot K_i=t_ia_i^{2}+a_i^{2^{e_i}+1},\qquad t_i\cdot K_i=t_i^{2}a_i+t_ia_i^{2^{e_i}}\;}
\qquad i=0,1,2\ \ (t_2=pt+ct),$$

six field polynomials, $n$ coordinates each. Measured: they lie in the $D=3$ row space, they are all
new modulo the equations, and the unexplained remainder is **0** at both $n=7$ and $n=9$.

**Why they do not raise the degree.** $K_i$ is quadratic and $a_i$, $t_i$ are linear, so a generic
product would be cubic. It is not, because $a_i^2$ and $t_i^2$ are $\mathbb{F}_2$-linear (Frobenius)
and $\operatorname{wt}(2^{e}+1)=2$: every monomial of $a_iK_i$ and $t_iK_i$ has exponent weight sum
$2$. **No cancellation is involved at all** — unlike the $D=4$ falls of §2, this is a pure
degree-of-the-exponent phenomenon.

**Two consequences worth stating.**

1. It is the **$K$ equations alone** that produce all $6n$; the $B$ equations contribute none. The
   guess in the batch brief --- "$t_i\cdot(B\text{-equation})$ against $a_i\cdot(K\text{-equation})$"
   --- is half right: both generators are multiples of $K_i$.
2. It explains the control. A **random** bilinear system of the same shape has no Frobenius structure,
   so multiplying an equation by a variable really does give a cubic --- which is exactly why E06
   measures **0** new elements for the random bilinear control at $D=3$ and $D=4$, against the kernel
   model's $6n$. The gap is not "the kernel model is more degenerate"; it is "$a^2$ is linear".

**What this is not.** The S-box relations $P_1^{(i)}$ are *implied* by the model through the identity
$B_i+x_i\,K_i=a_i\,P_1^{(i)}$ (a one-line check: the $a_i^{2^e}x_i$ terms cancel), but $P_1^{(i)}$
itself is **not** in the $D=3$ row space, and the $P_1^{(i)}$ explain **none** of the $6n$ new
quadratics (measured). Only the degree-3 multiple $a_iP_1^{(i)}$ is there.

---

## 4. Open problems, with the search already done

1. **[OPEN] The $5n$ non-trivial syzygies at $D=3$.** §1. The $D=4$ refutation (§2.4) suggests the
   same caution: try the leave-one-out analysis before assuming a two-family identity exists.
2. **[OPEN] The sign asymmetry $3e\equiv-1$ vs $3e\equiv+1$.** The first gives $2n+1$ falls, the
   second the generic $2n$ ($(11,4),(13,9),(14,5),(16,11)$ all satisfy $3e\equiv+1$ and all give
   $2n$). Every congruence family $k\cdot e\equiv\pm1,\pm2$ for $k\le6$ — for $e$ and for $e^{-1}$ —
   was tabulated in `results/E03_summary.md`; none of the others separates the two classes. A proof
   presumably has to exhibit the extra element, which §2.4 says is not a single-family object.
3. **[OPEN] Why $e=2$ leaves $2n$ cubics of $I(V_e)$ unreachable at $D=4$, while $e=2^{-1}$ leaves
   none.** $\mathrm{falls}_4$ is symmetric under $e\leftrightarrow e^{-1}$; `outside_ideal` is not
   (measured at $n=11,13,15$). The obvious candidate — that the transposition
   $(x,y)\mapsto(y,x)$ relates the two — is not available, because $P_1$ is not symmetric in $x$ and
   $y$. **Search range for a counterexample:** all admissible $e$ with $n\le16$ have been measured
   and only $e=2$ shows the defect; the next test is $n=17..20$ (server), where $e=2$ is admissible
   for odd $n$ only.
4. **[CONJ] The fall-count class rule** (Table in `results/E03_summary.md`): a fit to 36 runs at
   $n\le16$, untested above. Cheapest falsification: `sweep.sh 17 18 19 20` then `summarise.py`,
   which reports mismatches automatically.
5. **[CONJ] The reused-key coupling $C(q,2)\cdot n$** (`results/E09_theory_reused_key.csv`): verified
   in the non-saturated runs at $q=2,3,4$; the two most saturated runs deviate, which is expected but
   not proved to be the cause. Next: $q=2..4$ at $n=13$ with $\ell=2$, where saturation is lower.
