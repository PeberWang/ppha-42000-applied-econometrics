# Lecture 4 — Linear Space Theory

*PPHA 42000 Applied Econometrics I, Steven N. Durlauf, Fall 2026*

## Why Linear Spaces

Linear statistical models are ubiquitous in theoretical and empirical work, and the mathematics behind them is the theory of linear spaces (sometimes called vector spaces). The statistical models we study live in particular types of linear spaces, so the abstract structure developed here — following Royden (1988) and Simmons (1963), with Ash (1976) ch. 3 as a supplement — is the foundation for everything that follows in linear regression theory.

---

## Basic Definitions

**Definition 4.1 (Linear space).** Let $\Gamma$ denote a non-empty set. For each pair of elements $\gamma_i$ and $\gamma_j$ contained in $\Gamma$, assume there exists an operation $+$, called addition, which combines pairs of elements to produce an element $\gamma_i + \gamma_j$ also contained in $\Gamma$. Suppose the following hold:

$$
\text{i. } \gamma_i + \gamma_j = \gamma_j + \gamma_i \qquad \text{ii. } (\gamma_i + \gamma_j) + \gamma_k = \gamma_i + (\gamma_j + \gamma_k)
$$

$$
\text{iii. } \Gamma \text{ contains a unique element } 0 \text{ such that } 0 + \gamma = \gamma \ \forall \gamma \in \Gamma \qquad \text{iv. } \text{for each } \gamma \in \Gamma \ \exists \ (-\gamma) \text{ such that } \gamma + (-\gamma) = 0
$$

Further, for any scalar $c$ from the complex numbers $\mathbb{C}$, assume $c$ may be combined with any $\gamma \in \Gamma$ to produce an element $c\gamma \in \Gamma$ (scalar multiplication), with:

$$
\text{v. } c(\gamma_i + \gamma_j) = c\gamma_i + c\gamma_j \qquad \text{vi. } (c + d)\gamma = c\gamma + d\gamma \qquad \text{vii. } (cd)\gamma = c(d\gamma) \qquad \text{viii. } 1\gamma = \gamma
$$

Then the set $\Gamma$ together with the addition and multiplication operators defines a **linear space**.

Two comments. First, one can replace $\mathbb{C}$ with the real numbers $\mathbb{R}$ and define a linear space analogously; we use real and complex weights in different contexts. Second, the operators are defined by how they act on elements of $\Gamma$, so the definition is quite general — we will work with linear spaces of real-valued vectors, of functions, and of random variables.

In essence, a linear space must be **closed** under addition of any pair of its elements and under scalar multiplication. The definition does **not** require completeness: a space is complete when the limit of every (Cauchy) sequence of elements in the space is also an element of the space. Completeness matters for optimization — we need limits to exist inside the space to guarantee that minimizers are attainable.

---

## Norms, Convergence, and Banach Spaces

To discuss completeness rigorously, the space must be equipped with a notion of distance.

**Definition 4.2 (Norm).** For a linear space $\Gamma$, a norm $\|\cdot\|$ is a mapping from $\Gamma$ to $\mathbb{R}$ such that

$$
\text{i. } \|\gamma\| \geq 0; \ \|\gamma\| = 0 \Rightarrow \gamma = 0 \qquad \text{ii. } \|c\gamma\| = |c|\,\|\gamma\| \ \forall c \in \mathbb{C} \qquad \text{iii. } \|\gamma_i + \gamma_j\| \leq \|\gamma_i\| + \|\gamma_j\|
$$

The first condition makes the norm interpretable as a length; the second requires linearity in scale; the third is the **triangle inequality** — a hint that intuitions from 3-dimensional Euclidean geometry carry over to general spaces.

**Definition 4.3 (Normed linear space).** A normed linear space is a linear space equipped with a norm.

With a norm, one can discuss convergence of sequences of elements.

**Definition 4.4 (Cauchy sequence).** The sequence $\{\gamma_i\}$ is a Cauchy sequence if, for any $\varepsilon > 0$, there exists an $I$ (which may depend on $\varepsilon$) such that $\|\gamma_i - \gamma_j\| < \varepsilon$ for all $i, j > I$.

**Definition 4.5 (Limit of a Cauchy sequence).** The element $\gamma$ is the limit of the Cauchy sequence $\{\gamma_i\}$ if, for any $\varepsilon > 0$, there exists an $I$ such that $\|\gamma - \gamma_i\| < \varepsilon$ for all $i > I$.

**Definition 4.6 (Completeness).** A linear space $\Gamma$ is complete if the limit of every Cauchy sequence in the space is an element of the space.

It is not automatic that the limit point of a Cauchy sequence lies in the normed linear space. Complete normed linear spaces form an important class:

**Definition 4.7 (Banach space).** A Banach space is a complete normed linear space.

---

## Inner Products and Hilbert Spaces

The norm measures the magnitude of an element $\gamma_i$ and the distance $\|\gamma_i - \gamma_j\|$ between two elements. To characterize the *relationship* between pairs of elements (angles, orthogonality), we need an inner product. For a complex number $c$, $\bar{c}$ denotes its complex conjugate.

**Definition 4.8 (Inner product).** An inner product $\langle \cdot, \cdot \rangle$ is a mapping from $\Gamma \times \Gamma$ to $\mathbb{C}$ such that

$$
\text{i. } \langle a\gamma_i + b\gamma_j, \gamma_k \rangle = a\langle \gamma_i, \gamma_k \rangle + b\langle \gamma_j, \gamma_k \rangle \qquad \text{ii. } \langle \gamma_i, \gamma_j \rangle = \overline{\langle \gamma_j, \gamma_i \rangle}
$$

$$
\text{iii. } \langle \gamma_i, \gamma_i \rangle \geq 0 \qquad \text{iv. } \langle \gamma_i, \gamma_i \rangle = 0 \Rightarrow \gamma_i = 0
$$

**Definition 4.9 (Hilbert space).** A Hilbert space is a Banach space with an inner product such that the inner product defines the norm:

$$
\langle \gamma_i, \gamma_i \rangle = \|\gamma_i\|^2
$$

Hilbert spaces provide the foundation for the theory of linear models (and of time series): they are complete inner-product spaces, so both optimization (via completeness) and orthogonality arguments (via the inner product) are available.

---

## Orthogonality and Projections

**Definition 4.10 (Orthogonality).** Two elements $\gamma_i$ and $\gamma_j$ are orthogonal, written $\gamma_i \perp \gamma_j$, if

$$
\langle \gamma_i, \gamma_j \rangle = 0
$$

Orthogonality allows one to construct decompositions of elements of Hilbert spaces, and hence of the spaces themselves. Take any two elements $\gamma_i, \gamma_j$ and consider the element

$$
\gamma_i - \frac{\langle \gamma_i, \gamma_j \rangle}{\langle \gamma_j, \gamma_j \rangle}\, \gamma_j
$$

This element belongs to the space (linearity), and it is orthogonal to $\gamma_j$, since

$$
\left\langle \gamma_i - \frac{\langle \gamma_i, \gamma_j \rangle}{\langle \gamma_j, \gamma_j \rangle}\gamma_j,\ \gamma_j \right\rangle = \langle \gamma_i, \gamma_j \rangle - \frac{\langle \gamma_i, \gamma_j \rangle}{\langle \gamma_j, \gamma_j \rangle}\langle \gamma_j, \gamma_j \rangle = 0
$$

This construction — subtracting off the component of $\gamma_i$ along $\gamma_j$ — is the atomic form of every projection argument in econometrics.

**Theorem 4.1 (Decomposability of a Hilbert space into orthogonal subspaces).** Suppose $\Gamma$ and $\Gamma_1$ are Hilbert spaces with $\Gamma_1 \subseteq \Gamma$. Then there exists a unique Hilbert space $\Gamma_2 \subseteq \Gamma$ such that

$$
\text{i. } \Gamma = \Gamma_1 \oplus \Gamma_2 \qquad \text{ii. } \Gamma_1 \perp \Gamma_2
$$

The operator $\oplus$ is the **direct sum**: it produces the space whose elements are all linear combinations of elements of the original spaces, together with the limits of all such combinations.

**Corollary 4.1 (Decomposition of an element into orthogonal components).** Suppose $\Gamma$ and $\Gamma_1$ are Hilbert spaces with $\Gamma_1 \subseteq \Gamma$. For a given element $\gamma \in \Gamma$, there exists a unique pair $\gamma_1 \in \Gamma_1$ and $\gamma_2 \in \Gamma_2$ (with $\Gamma_2$ as in Theorem 4.1) such that

$$
\text{i. } \gamma = \gamma_1 + \gamma_2 \qquad \text{ii. } \gamma_1 \perp \gamma_2
$$

Here $\gamma_1$ is called the **projection** of $\gamma$ onto $\Gamma_1$. The projection has an optimization interpretation: it solves

$$
\min_{\xi \in \Gamma_1} \|\gamma - \xi\|^2
$$

i.e. the projection is the element of the subspace $\Gamma_1$ **closest** to $\gamma$. To see this, write

$$
\|\gamma - \xi\|^2 = \|\gamma_1 + \gamma_2 - \xi\|^2 = \|\gamma_1 - \xi\|^2 + \|\gamma_2\|^2 + \langle \gamma_1 - \xi, \gamma_2 \rangle + \langle \gamma_2, \gamma_1 - \xi \rangle = \|\gamma_1 - \xi\|^2 + \|\gamma_2\|^2
$$

where the last equality uses the orthogonality of $\Gamma_1$ and $\Gamma_2$ (note $\gamma_1 - \xi \in \Gamma_1$). Since $\xi$ affects the objective only through $\|\gamma_1 - \xi\|^2$, the solution is immediate: $\xi = \gamma_1$. This is the **projection theorem**, and it is why least squares is a projection.

---

## Orthonormal Bases

In $\mathbb{R}^n$, every element is a linear combination of the axis vectors $(1, 0, \dots)$, $(0, 1, \dots)$, etc. This property generalizes to all Hilbert spaces.

**Definition 4.11 (Orthonormal set).** A subset $\mathcal{S}$ of a Hilbert space is an orthonormal set if (i) each element of $\mathcal{S}$ is orthogonal to every other element, and (ii) the norm of every element equals $1$.

**Definition 4.12 (Orthonormal basis).** An orthonormal set $\mathcal{S}$ is an orthonormal basis for a Hilbert space $H$ if it is not a proper subset of any other orthonormal set in the space (maximality).

A set of elements **spans** a space if, starting from the set, linear combinations of the elements and their limits together constitute the space. Orthonormal bases span:

**Theorem 4.2.** An orthonormal basis of a given Hilbert space $\Gamma$ spans the space.

**Theorem 4.3 (Existence).** Every Hilbert space contains an orthonormal basis.

**Theorem 4.4 (Countable bases).** A Hilbert space contains a countable orthonormal basis if and only if the space is **separable** (i.e. it contains a countable dense subset). The spaces we study will always be separable.

Countable orthonormal bases are the "axes" of a Hilbert space, analogous to the unit vectors of $\mathbb{R}^k$. When $\mathcal{S}$ is a countable orthonormal basis for $\Gamma$, a typical element $\gamma$ can be represented (exactly, or arbitrarily well approximated — the difference being that a complete space also contains the limits of all such sums) as

$$
\gamma = \sum_{i=1}^{\infty} a_i s_i
$$

The coefficients are determined by projecting $\gamma$ onto each 1-dimensional subspace spanned by a basis element:

$$
a_i = \langle \gamma, s_i \rangle
$$

since $\langle s_i, s_i \rangle = 1$ by assumption. These are sometimes known as **Fourier coefficients**.

---

## Applications

**Euclidean space.** Consider $\mathbb{R}^k$. For vectors $x$ and $y$, the inner product is $\langle x, y \rangle = \sum_i x_i y_i$ with associated norm $\|x\| = \left( \sum_i x_i x_i \right)^{1/2}$. To see the decomposition theorem at work, take the Hilbert space generated around $(1, 0, 0, \dots, 0)$ (generate = add all linear combinations and their limits). Its orthogonal complement in $\mathbb{R}^k$ is the space spanned by $(0, 1, 0, \dots, 0), (0, 0, 1, \dots, 0), \dots, (0, 0, \dots, 0, 1)$.

**Example 4.1 — Pythagorean theorem.** Suppose an element $x$ of a Hilbert space is the sum of two orthogonal elements $y$ and $z$. Then the squared norm of $x$ obeys

$$
\langle x, x \rangle = \langle y + z, y + z \rangle = \langle y, y \rangle + \langle z, z \rangle
$$

since by orthogonality $\langle y, z \rangle = \langle z, y \rangle = 0$. This is the Pythagorean theorem for an arbitrary Hilbert space, not just $\mathbb{R}^2$ — and it is the engine behind variance decompositions in regression analysis.

**Example 4.2 — Cauchy–Schwarz inequality.** For any scalar $a$, the inner product of $x + ay$ with itself can be written

$$
\langle x + ay, x + ay \rangle = \langle x, x \rangle + \bar{a}\langle y, x \rangle + a\langle x, y \rangle + |a|^2 \langle y, y \rangle \geq 0
$$

Choose $a = -\langle x, y \rangle / \langle y, y \rangle$. Substituting,

$$
\langle x + ay, x + ay \rangle = \langle x, x \rangle - \frac{\langle x, y \rangle \langle y, x \rangle}{\langle y, y \rangle} - \frac{\langle x, y \rangle \langle y, x \rangle}{\langle y, y \rangle} + \frac{\langle x, y \rangle \langle y, x \rangle \langle y, y \rangle}{\langle y, y \rangle^2} = \langle x, x \rangle - \frac{\langle x, y \rangle \langle y, x \rangle}{\langle y, y \rangle} \geq 0
$$

Rearranging terms gives

$$
\langle x, x \rangle \langle y, y \rangle \geq |\langle x, y \rangle|^2
$$

with equality if and only if $x = ay$ (the "if and only if" claim is worth verifying yourself). This is the Cauchy–Schwarz inequality for an arbitrary Hilbert space; applied to random variables with $\langle x, y \rangle = \text{cov}(x, y)$, it bounds covariances by variances and hence bounds correlation in $[-1, 1]$.

---

## Difficulty checklist

- **The linear-space axioms are closure conditions.** The content of Definition 4.1 is that $\Gamma$ is closed under $+$ and scalar multiplication; axioms (i)–(viii) are bookkeeping. Spaces of random variables (with $\langle x, y \rangle = \text{cov}(x,y)$) qualify — this is the move that makes regression a projection.
- **Norm $\neq$ inner product.** Every inner product induces a norm via $\|\gamma\|^2 = \langle \gamma, \gamma \rangle$, but not every norm comes from an inner product. Hilbert space = Banach space *whose norm comes from an inner product*.
- **Completeness is a separate assumption.** A normed linear space need not contain the limits of its Cauchy sequences. Completeness is what guarantees that the minimizer in $\min_{\xi \in \Gamma_1} \|\gamma - \xi\|^2$ actually exists in $\Gamma_1$.
- **Inner-product conjugate symmetry.** $\langle \gamma_i, \gamma_j \rangle = \overline{\langle \gamma_j, \gamma_i \rangle}$, not plain symmetry, when scalars are complex. Forgetting the conjugate breaks the Cauchy–Schwarz proof.
- **Projection formula.** The component of $\gamma_i$ along $\gamma_j$ is $\frac{\langle \gamma_i, \gamma_j \rangle}{\langle \gamma_j, \gamma_j \rangle} \gamma_j$; the residual $\gamma_i - \frac{\langle \gamma_i, \gamma_j \rangle}{\langle \gamma_j, \gamma_j \rangle}\gamma_j$ is orthogonal to $\gamma_j$ by construction. Every later projection formula (including OLS) is this idea in matrix or moment form.
- **Direct sum has two clauses.** $\Gamma = \Gamma_1 \oplus \Gamma_2$ requires *both* that the spaces generate $\Gamma$ and that $\Gamma_1 \perp \Gamma_2$; the decomposition of an element $\gamma = \gamma_1 + \gamma_2$ is then **unique**.
- **Projection solves the minimization because cross terms vanish.** In $\|\gamma - \xi\|^2 = \|\gamma_1 - \xi\|^2 + \|\gamma_2\|^2$, orthogonality kills $\langle \gamma_1 - \xi, \gamma_2 \rangle$; without it the minimizer would not be $\gamma_1$.
- **Basis vs. spanning set.** An orthonormal basis is a *maximal* orthonormal set; only then does Theorem 4.2 (it spans the space) apply. Finite orthonormal sets need not span.
- **Countable basis $\Leftrightarrow$ separable.** Theorem 4.4 is an "if and only if"; nonseparable Hilbert spaces have only uncountable bases. Our spaces are separable, so $\gamma = \sum_i a_i s_i$ with $a_i = \langle \gamma, s_i \rangle$ is legitimate.
- **Cauchy–Schwarz equality case.** Equality in $\langle x,x\rangle\langle y,y\rangle \geq |\langle x,y\rangle|^2$ holds iff $x$ and $y$ are linearly dependent ($x = ay$). This is the abstract source of "correlation $= \pm 1$ iff perfectly linearly related."
