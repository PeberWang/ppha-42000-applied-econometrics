# Lecture 3 — Structures, Models and Identification

*PPHA 42000 Applied Econometrics I, Steven N. Durlauf, Fall 2026*

## The Identification Question

Identification is a fundamental question in any empirical exercise. Conceptually, to study identification is to ask whether, **in principle**, data can yield information about objects of interest. The definition is deliberately abstract about what the objects of interest are: a given data set, viewed through the lens of a statistical model, may allow some objects to be identified while others are not. Identification is therefore a property of a (data, model, parameter) triple — not of a data set alone.

---

## The Koopmans–Marschak Framework

Koopmans (1953) and Marschak (1950) provide the general definition of identification in statistical models. Consider a data set comprised of realizations of a set of random variables $Z$ with associated probability measure $\mu_z$. This measure may be interpreted as the **reduced form** of the data. Typically one implicitly assumes infinitely many realizations of an ergodic data-generating process, so that $\mu_z$ is itself revealed by the data — but this is not essential to defining identification.

Koopmans and Marschak consider the relationship between this probability measure and a **structure**

$$
(M, \alpha)
$$

where $M$ is a model of the data and $\alpha$ is a set of possible parameter values, i.e. $\alpha \in \mathcal{A}$ for some parameter space $\mathcal{A}$. From a model one can construct a probability measure $\mu_{z \mid M, \alpha}$ that describes probabilities conditional on the model and the parameter vector.

The identification problem, in essence, asks whether it is possible to learn about $\alpha$ from $\mu_z$ via comparison of this empirical measure with the theoretically informed measure $\mu_{z \mid M, \alpha}$. One can also ask whether $\mu_z$ allows **model-free** identification of a parameter.

---

## Two Definitions of Identification

**Definition 1 (Global identification).** A model is globally identified at the parameter value $\alpha$ if

$$
\forall \tilde{\alpha} \in \mathcal{A}, \quad \mu_{z \mid M, \alpha} = \mu_{z \mid M, \tilde{\alpha}} \Rightarrow \alpha = \tilde{\alpha}
$$

That is, no other parameter value in the space produces the same probability restrictions on the data.

**Definition 2 (Local identification).** A model is locally identified at $\alpha$ if there exists an $\varepsilon > 0$ such that

$$
\forall \tilde{\alpha} \text{ with } \|\alpha - \tilde{\alpha}\| \leq \varepsilon, \quad \mu_{z \mid M, \alpha} = \mu_{z \mid M, \tilde{\alpha}} \Rightarrow \alpha = \tilde{\alpha}
$$

If the global condition holds for **all** parameters $\alpha$, the structure is globally identified; an analogous generalization exists for local identification. When two distinct parameter values induce the same probability measure on the data, they are said to be **observationally equivalent** — the data cannot tell them apart.

**Three comments.**

1. Identification may hold for some parameter values but not for others — it is a pointwise property.
2. There may be identification for **subsets** of $\alpha$. This matters because we may not care about some aspects of the unknowns. Marschak (1953) originally emphasized that one should employ the *minimum* assumptions needed for identification; James Heckman calls this **Marschak's maxim**: estimators should be selected on the basis of their ability to answer well-posed economic problems with minimal assumptions.
3. Identification is **distinct from consistency** of estimators.

---

## Example 3.1 — The Linear Model

Consider the joint density of $K+1$ random variables, $y$ and $x_1, \dots, x_K$, denoted collectively as $x$. Assume all data have zero expected value. The set of structures of interest is

$$
y_i = x_i \beta + \varepsilon_i
$$

where $\varepsilon$ is $N(0,1)$, independent of $x$.

The data reveal three second moments: $\text{var}(y)$, $\text{var}(x)$, and $\text{cov}(x, y)$. The structure implies

$$
\text{cov}(x, y) = \text{var}(x)\, \beta \quad \Longrightarrow \quad \beta = \text{var}(x)^{-1}\, \text{cov}(x, y)
$$

when the inverse exists. The unknown parameters $\beta$ are therefore **globally identified if $\text{var}(x)^{-1}$ exists**. This is exactly the identification condition for ordinary least squares, applied to populations.

**Did knowledge of the density of the $\varepsilon$'s facilitate identification? No.** The normality assumption was unnecessary — an example of Marschak's maxim in action. The minimal assumption set changes across objects of interest.

---

## Example 3.2 — Identification Without Consistency

Consider ordinary least squares with a deterministically shrinking regressor:

$$
y_t = x_t \beta + \varepsilon_t, \quad x_t = \tfrac{1}{t}, \quad t = 1, \dots, L
$$

with $\varepsilon \sim N(0,1)$, independent of $x$. One knows the density of the OLS estimator $\hat{\beta}$. As $L \Rightarrow \infty$, the density converges to

$$
N\!\left(\beta,\ \frac{6}{\pi^2}\right)
$$

since $\lim_{L \Rightarrow \infty} \sum_{t=1}^{L} t^{-2} = \pi^2/6$ (the Basel problem, solved by Euler in 1734).

Note that the OLS estimator is **inconsistent**: it does not converge to the true $\beta$ with arbitrarily high probability. The additional data become less and less informative about $\beta$ because the variation of $1/t$ decreases as $t$ grows.

Yet $\beta$ is **still identified**: different $\beta$'s have different implications for the population covariance between regressor and regressand. Identification is a population-level logical property; consistency is about estimator behavior as samples grow.

---

## Beyond Binary Identification: Identification as the Limit of Imprecision

The linear examples hint at an important limit to treating identification as a binary property. Consider cases where $\text{var}(x)$ is invertible but the regressors are close to linearly dependent; the system then has very large standard errors.

Since $\text{var}(x)^{-1}$ is positive definite, it can be represented as

$$
\text{var}(x)^{-1} = \sum_{k=1}^{K} \lambda_k\, v_k v_k^\top
$$

where the $\lambda_k$'s are the eigenvalues and $v_k$ the associated eigenvectors. The eigenvalues of $\text{var}(x)^{-1}$ are the **inverses** of the eigenvalues of $\text{var}(x)$. If $\text{var}(x)$ has an eigenvalue of $0$, the matrix is not invertible; if the smallest eigenvalue is close to $0$, the inverse has arbitrarily large elements.

Hence **nonidentification is the limit of imprecision of information** — and ultimately it is the imprecision, not the binary label, that matters for empirical work.

---

## Partial Identification

With Charles Manski as the leading advocate (Manski 1995 is the readable overview), interest has developed in models that are **partially identified**. Intuitively: even when sets of parameters produce the same empirical implications (the same probability restrictions on the data), the data may still shrink the set of parameter values consistent with the data to a strict subset of the original parameter space.

Formally, the parameter $\alpha$ is partially identified if

$$
\forall \tilde{\alpha} \in \mathcal{A}, \quad \mu_{z \mid M, \alpha} = \mu_{z \mid M, \tilde{\alpha}} \Rightarrow \tilde{\alpha} \in \mathcal{A}' \subset \mathcal{A}
$$

---

## Example 3.3 — Treatment Effects

Consider measuring the effect of an afterschool program on the probability of graduating high school. Let $g_i = 1$ if student $i$ graduated, $g_i = 0$ otherwise. Participation was voluntary; denote participation by $h_i = 1$.

Suppose a policymaker considers making the program **mandatory** and wants to know $\mu(g_i = 1 \mid p_i)$, the probability a student graduates when treated under the mandatory regime ($p_i = 1$ for everybody; the notation emphasizes the new policy regime). This probability decomposes as

$$
\mu(g_i = 1 \mid p_i) = \mu(g_i = 1 \mid p_i = 1, h_i = 1)\,\mu(h_i = 1) + \mu(g_i = 1 \mid p_i = 1, h_i = 0)\,\mu(h_i = 0)
$$

Interpretation: a randomly selected student forced to enroll is either someone who would have volunteered or someone who would not have, and one must allow different success probabilities for each type.

**What the data know.** The volunteering probabilities $\mu(h_i = 1)$ and $\mu(h_i = 0) = 1 - \mu(h_i = 1)$ are identified from the data. And if making the program mandatory does not affect those who would have volunteered anyway, we can take

$$
\mu(g_i = 1 \mid p_i = 1, h_i = 1) = \mu(g_i = 1 \mid h_i = 1)
$$

as known. But $\mu(g_i = 1 \mid p_i = 1, h_i = 0)$ is a **counterfactual that has never been observed**: how would non-volunteers fare if forced? The data cannot speak to this probability. One could assume $\mu(g_i = 1 \mid p_i = 1, h_i = 0) = \mu(g_i = 1 \mid h_i = 1)$ — that volunteering status is irrelevant to treatment effects — but this is the classic **selection problem** and is usually untenable.

**Bounds via the law of probability.** The missing probability must lie between $0$ and $1$. Bounding $\mu(g_i = 1 \mid p_i)$ by the cases where $\mu(g_i = 1 \mid p_i = 1, h_i = 0)$ equals $0$ or $1$ gives

$$
\mu(g_i = 1 \mid h_i = 1)\,\mu(h_i = 1) \;\leq\; \mu(g_i = 1 \mid p_i) \;\leq\; \mu(g_i = 1 \mid h_i = 1)\,\mu(h_i = 1) + \mu(h_i = 0)
$$

The parameter is partially identified: the data plus minimal assumptions deliver an informative interval rather than a point.

**Monotone treatment effects.** Treating $0$ as the lower bound allows the possibility that $\mu(g_i = 1 \mid p_i = 1, h_i = 0) = 0$, i.e. that forcing non-volunteers into the program might harm them. Suppose instead that $\mu(g_i = 0 \mid h_i = 0) > 0$ and one assumes the program cannot hurt a non-volunteer:

$$
\mu(g_i = 1 \mid p_i = 1, h_i = 0) \geq \mu(g_i = 1 \mid h_i = 0)
$$

This is a **monotone treatment effect** assumption. The bounds then tighten to

$$
\mu(g_i = 1 \mid h_i = 1)\,\mu(h_i = 1) + \mu(g_i = 1 \mid h_i = 0)\,\mu(h_i = 0) \;\leq\; \mu(g_i = 1 \mid p_i) \;\leq\; \mu(g_i = 1 \mid h_i = 1)\,\mu(h_i = 1) + \mu(h_i = 0)
$$

The general lesson (Heckman's and Lewbel's perspective): credible empirical work proceeds by asking what can be learned under weak, defensible assumptions — bounds under monotonicity — rather than by imposing strong assumptions solely to obtain point identification.

---

## Bayesian Perspectives on Identification

The Bayesian perspective on parameters involves calculating the posterior

$$
\mu(\alpha \mid d) \propto \mu(d \mid \alpha)\, \mu(\alpha)
$$

Identification enters through the relationship between $\mu(d \mid \alpha)$ and $\mu(\alpha)$. The analogy to failure of identification occurs **when the data cannot update beliefs**: if the likelihood is flat in a subset of the $\alpha$'s, the posterior over that subset equals the prior, whatever the data. This is the Bayesian parallel of the identification failures discussed above.

---

## Difficulty checklist

- **Identification is defined relative to a model and a parameter of interest.** The same data set can identify some objects and not others; never say "the data identify" without specifying the structure $(M, \alpha)$ and the object.
- **Global vs. local identification.** Global: $\mu_{z \mid M, \alpha} = \mu_{z \mid M, \tilde{\alpha}} \Rightarrow \alpha = \tilde{\alpha}$ for all $\tilde{\alpha} \in \mathcal{A}$. Local: the same implication only in an $\varepsilon$-neighborhood. A model can be locally but not globally identified (e.g. multiple isolated roots).
- **Observational equivalence** is the operative failure mode: two parameter values producing the same probability measure on observables cannot be distinguished by any estimator, however clever.
- **Identification $\neq$ consistency.** Example 3.2 ($x_t = 1/t$) has $\beta$ identified but OLS inconsistent, because $\sum_t x_t^2 \to \pi^2/6 < \infty$ so new data add vanishing information. Identification is about the population mapping; consistency is about estimator limits.
- **Marschak's maxim is easy to misapply.** In Example 3.1, normality of $\varepsilon$ is unnecessary — but do not conclude distributional assumptions are *never* needed; the minimal set depends on the object of interest.
- **Eigenvalue logic.** $\text{var}(x)^{-1}$ exists iff no eigenvalue of $\text{var}(x)$ is zero; near-zero eigenvalues give arbitrarily large elements in the inverse. Nonidentification is the *limit* of imprecision — borderline cases are practically useless even when technically identified.
- **Partial identification is not failure.** Bounds like $\mu(g_i = 1 \mid h_i = 1)\mu(h_i = 1) \leq \mu(g_i = 1 \mid p_i) \leq \mu(g_i = 1 \mid h_i = 1)\mu(h_i = 1) + \mu(h_i = 0)$ are informative findings, not consolation prizes.
- **Treatment-effect bounds depend on the policy regime.** The decomposition $\mu(g_i = 1 \mid p_i) = \mu(g_i = 1 \mid p_i = 1, h_i = 1)\mu(h_i = 1) + \mu(g_i = 1 \mid p_i = 1, h_i = 0)\mu(h_i = 0)$ mixes an identified term with a pure counterfactual; keep track of which piece the data can see.
- **Monotonicity assumptions tighten bounds only in one direction.** Assuming $\mu(g_i = 1 \mid p_i = 1, h_i = 0) \geq \mu(g_i = 1 \mid h_i = 0)$ raises the lower bound; it does nothing to the upper bound.
- **Bayesian parallel.** Failure of identification = likelihood flat in $\alpha$ over a region $\Rightarrow$ posterior equals prior there. Priors can mask nonidentification by producing a proper posterior even when the data say nothing.
