# Lecture 2 — Data and Decisions

*PPHA 42000 Applied Econometrics I, Steven N. Durlauf, Fall 2026*

## The Decision Problem

These notes describe basic ideas in decision theory under uncertainty, which provide links between econometric analysis and policy evaluation. Basic decision theory is built from four elements:

- **Data:** realizations $d \in D$
- **Unknown(s):** $\theta$
- **Choices:** $c \in C$
- **Loss function:** $\ell(\theta, d, c)$

The goal is to construct a **decision rule** $c(d)$, a mapping from $D$ to the set of choices $C$. The rule may be stochastic, but we ignore this possibility for expositional purposes.

---

## Statistical Decision Theory (Bayesian)

The standard version of statistical decision theory is formulated along Bayesian lines: all unknowns are associated with posterior probability densities, and a decision rule is determined by minimization of expected loss. One implicitly chooses $c(d)$ by solving, for each value of $d$:

$$
\min_{c \in C} \int \ell(\theta, d, c) \, \mu(\theta \mid d) \, d\theta \tag{2.1}
$$

**Constructing the posterior.** The key question facing an expected-loss-minimizing policy analyst is how to construct $\mu(\theta \mid d)$. Observe that:

$$
\mu(\theta \mid d) = \frac{\mu(\theta, d)}{\mu(d)} = \frac{\mu(d \mid \theta) \, \mu(\theta)}{\mu(d)} \tag{2.2}
$$

Since $\mu(d)$ is not a function of $\theta$ (it is a normalization ensuring probabilities sum to 1), one can rewrite (2.2) as:

$$
\mu(\theta \mid d) \propto \mu(d \mid \theta) \, \mu(\theta) \tag{2.3}
$$

This is the classic statement that the **posterior** is proportional to the product of the **likelihood** $\mu(d \mid \theta)$ and the **prior** $\mu(\theta)$. The prior represents the information the analyst possesses about the unknown before the data. For some types of unknowns (e.g., parameters) under regularity conditions, as the number of observations grows the likelihood dominates the prior, in the sense that the prior no longer affects the posterior.

### Loss functions and "standard" statistical exercises

The decision-theoretic formulation does not at first glance appear related to standard exercises such as producing a point estimate, but such exercises are reproduced in their essentials by suitable loss functions. Consider quadratic loss:

$$
\ell(\theta, d, c) = (\theta - c(d))^2 \tag{2.4}
$$

The expected-loss-minimizing decision rule under (2.4) is the posterior mean:

$$
c(d) = \mathbb{E}(\theta \mid d)
$$

(Proof left as an exercise.) This posterior mean is the Bayesian analog to a point estimate in frequentist analysis. Hypothesis testing can be cast in the same language: let $C$ have two elements (accept or reject the hypothesis) and let $\theta$ encode the truth or falsity of the hypothesis.

---

## Priors

The assignment of priors is problematic: in many contexts the analyst has no good basis for constructing $\mu(\theta)$. First, for many problems one may simply have nothing to say about the ex ante relative probabilities of different realizations. Second, the prior knowledge one does possess may not be naturally quantified. We focus on constructing priors under **ignorance** — no basis for discriminating between values of the unknown — using the example of the mean $\kappa$ and variance $\sigma^2$ of a sequence of i.i.d. random variables.

### Principle of insufficient reason

The classical approach, attributed to Bernoulli and Laplace (though the name postdates them): in the absence of any reason to discriminate between two values of an unknown, assign them equal prior probabilities.

- For the mean $\kappa$ with support $(-\infty, \infty)$ under ignorance, the prior must be constant: a uniform density on the real line, $\mu(\kappa) \propto 1$. This prior is **improper** — it does not integrate to 1.
- For the standard deviation $\sigma$ (the parameter typically studied), the support under ignorance is $[0, \infty)$. The standard ignorant prior is $\mu(\sigma) \propto \frac{1}{\sigma}$, motivated by the fact that it assigns a uniform prior to $\log \sigma$, whose support is $(-\infty, \infty)$.

### The invariance problem and Jeffreys prior

Uniform priors raise a general problem. Consider a parameter vector $\theta$ and a transformation $\phi = f(\theta)$. Assuming all measures have densities, for a given prior $\mu(\theta)$ the induced prior for $\phi$ must be:

$$
\mu(\phi) = \mu\left(f^{-1}(\phi)\right) \left| \frac{d f^{-1}(\phi)}{d\phi} \right| \tag{2.5}
$$

So the uniform prior is **not invariant** across nonlinear transformations — which makes little sense if the prior truly captures ignorance, since ignorance about $\theta$ presumably implies ignorance about $\phi$. One way to impose invariance is **Jeffreys prior**. Recall Fisher's information matrix:

$$
I(\theta) = -\mathbb{E}\left( \frac{\partial^2 \log \mu(d \mid \theta)}{\partial \theta^2} \right) = -\int_D \frac{\partial^2 \log \mu(d \mid \theta)}{\partial \theta^2} \, \mu(d \mid \theta) \, dd \tag{2.6}
$$

Jeffreys prior is defined as:

$$
\mu(\theta) \propto I(\theta)^{1/2} \tag{2.7}
$$

This prior is invariant in form under reparameterization, but raises interpretability difficulties: the prior depends on the form of the likelihood, which seems odd since the prior is supposed to describe beliefs about parameters.

### Example 2.1. Priors and posteriors for the mean of normal random variables

Suppose $x_i \sim N(\kappa, \sigma^2)$, $i = 1, \ldots, K$, with $\sigma^2$ known, and the prior for $\kappa$ is $N(\kappa_0, \sigma_0^2)$. Algebraic manipulation reveals the posterior density for the mean:

$$
\mu(\kappa \mid x_1, \ldots, x_K, \sigma^2) = N\left( \left( \frac{1}{\sigma_0^2} + \frac{K}{\sigma^2} \right)^{-1} \left( \frac{\kappa_0}{\sigma_0^2} + \frac{K \bar{x}}{\sigma^2} \right), \; \left( \frac{1}{\sigma_0^2} + \frac{K}{\sigma^2} \right)^{-1} \right) \tag{2.8}
$$

**Intuition:** the posterior mean is a weighted average of the prior mean $\kappa_0$ and the sample average $\bar{x}$, with weights reflecting their respective variances — more weight attaches to the more precise piece of information about $\kappa$. As $K \to \infty$, the effect of the prior on the posterior disappears: the posterior mean converges to the sample mean.

---

## Decision Theory Without Probabilities

Now consider decision criteria when probabilities are not available for the unobservables of interest.

### Minimax

The minimax solution chooses $c$ so that:

$$
\min_{c \in C} \max_{\theta} \ell(\theta, d, c) \tag{2.9}
$$

Minimax is usually associated with Abraham Wald (1950); an important axiomatization is due to Gilboa and Schmeidler (1989). The criterion is operationally equivalent to an extremely risk-averse policymaker.

**Rawls' difference principle.** Individuals rank distributions of economic outcomes in a society without knowing which outcome will be theirs; Rawls (1971) argues they will choose the distribution maximizing the utility of the worst-off person. As Arrow (1973) pointed out, this is a minimax argument, and can be understood as the limit of a welfarist analysis. Suppose individual utility is $u_i = -\omega_i^{-\alpha}$, where $\omega_i$ is a scalar measure of individual $i$'s outcome, and the social state $\omega$ is assigned welfare $W = \left( \sum_i \omega_i^{-\alpha} \right)^{-\frac{1}{\alpha}}$. Then as $\alpha \to \infty$, $W \to \min_i \omega_i$. The limit implies an arbitrary degree of risk aversion — and since this holds for any monotonic transformation of the social welfare function, it applies to a utilitarian welfare calculation as well.

### Minimax regret

An alternative that avoids assigning probabilities is to minimize the maximum **regret**. Proposed by Savage (1951) as a less conservative alternative to minimax; Charles Manski is its leading current advocate (see Manski 2008 for a conceptual defense and Brock and Durlauf 2015 for criticism). Regret is formally defined as:

$$
r(\theta, d, c) = \ell(\theta, d, c) - \min_{c \in C} \ell(\theta, d, c) \tag{2.10}
$$

The minimax regret solution is:

$$
\min_{c} \max_{\theta} r(\theta, d, c) \tag{2.11}
$$

**Example 2.3. Minimax, minimax regret, and expected loss.** Two actions $c_1, c_2$; two states of the world $\theta_1, \theta_2$:

| | $\theta_1$ | $\theta_2$ | max loss | max regret |
|---|---|---|---|---|
| $c_1$ | 24 | 12 | 24 | 5 |
| $c_2$ | 19 | 15 | 19 | 3 |

By both the minimax and minimax regret criteria, the policymaker should choose $c_2$. A Bayesian would assign probability $p$ ($1-p$) to $\theta_1$ ($\theta_2$) and compute expected losses; if $1-p$ is close enough to 1, the policymaker should choose $c_1$.

**Example 2.4. Minimax and minimax regret producing different choices.**

| | $\theta_1$ | $\theta_2$ | max loss | max regret |
|---|---|---|---|---|
| $c_1$ | 24 | 23 | 24 | 8 |
| $c_2$ | 25 | 15 | 25 | 1 |

Here the minimax choice is $c_1$ whereas the minimax regret choice is $c_2$. This illustrates the intuitive appeal of minimax regret: the criterion focuses on differences in policy effects rather than absolute levels.

**Example 2.5. Minimax regret and the IIA violation.** Modify Example 2.3 by adding a third policy $c_3$:

| | $\theta_1$ | $\theta_2$ | max loss | max regret |
|---|---|---|---|---|
| $c_1$ | 24 | 12 | 24 | 12 |
| $c_2$ | 19 | 15 | 19 | 15 |
| $c_3$ | 50 | 0 | 50 | 31 |

The introduction of $c_3$ makes the minimax regret choice $c_1$: the mere availability of $c_3$ has changed the ranking of $c_1$ and $c_2$. This violates the **independence of irrelevant alternatives (IIA)**, often taken as a natural axiom of decisionmaking. The recognition that minimax regret violates IIA is due to Chernoff (1954). One interpretation is that context matters for decisionmaking — an idea found in behavioral economics — though it is unclear how closely context-dependent choice findings map to these IIA violations.

**Example 2.6. Another minimax regret IIA violation.** Compare the two-policy problem:

| | $\theta_1$ | $\theta_2$ | max loss | max regret |
|---|---|---|---|---|
| $c_1$ | 30 | 10 | 30 | 9 |
| $c_2$ | 21 | 21 | 21 | 11 |

with the three-policy problem:

| | $\theta_1$ | $\theta_2$ | max loss | max regret |
|---|---|---|---|---|
| $c_1$ | 30 | 10 | 30 | 20 |
| $c_2$ | 21 | 21 | 21 | 11 |
| $c_3$ | 10 | 30 | 30 | 20 |

$c_2$ is the minimax regret choice in the two-policy case; but with three policies, $c_2$ is no longer chosen, and if $c_3$ ($c_1$) were not available, $c_1$ ($c_3$) would be the solution. This seems especially paradoxical since $c_1$ and $c_3$ have symmetric structures.

**Axiomatic repair.** Stoye (2006a,b) shows that if one forgoes IIA, minimax regret can be derived under an axiom system including "independence of never-optimal alternatives": adding a policy that is not optimal in any state of nature cannot change one's choice among the original set.

### Hurwicz criterion

One can mix the approaches — for example, combining expected loss and minimax considerations:

$$
\min_{c \in C} \left[ \alpha \max_{\theta} \ell(\theta, d, c) + (1 - \alpha) \, \mathbb{E} \ell(\theta, d, c) \right] \tag{2.12}
$$

This was originally proposed by Hurwicz (1951).

---

## Decision Rules and Risk

The **risk function**, defined as:

$$
R(\theta, c(D)) = \int_D \ell(\theta, d, c(d)) \, \mu(d \mid \theta) \, dd \tag{2.13}
$$

is a metric for the performance of a rule across the data space, conditional on $\theta$; performance is weighted across the data space by probabilities.

A rule $c(\cdot)$ is **(strictly) dominant** if it produces (strictly) lower risk than any alternative rule for all values of $\theta$. A rule is **inadmissible** if there exists an alternative with lower risk for all $\theta$. There is no guarantee that a unique dominant rule exists.

The **Bayes risk** of a rule is:

$$
R = \int_{\Theta} R(\theta, c(d)) \, \mu(\theta) \, d\theta = \int_{\Theta} \int_D \ell(\theta, d, c(d)) \, \mu(d \mid \theta) \, \mu(\theta) \, dd \, d\theta \tag{2.14}
$$

Recall that $\mu(d \mid \theta) \mu(\theta) = \mu(\theta \mid d) \mu(d)$. Exchanging the order of integration (fine under uninteresting regularity assumptions) and substituting, the Bayes risk becomes:

$$
R = \int_D \int_{\Theta} \ell(\theta, d, c(d)) \, \mu(\theta \mid d) \, \mu(d) \, d\theta \, dd \tag{2.15}
$$

Minimizing the double integral evidently occurs by minimizing the inner integral at each data point — which is exactly the original Bayesian decision-theory solution. A general result: any decision rule that minimizes Bayes risk relative to a **proper** prior is admissible. A form of the converse also holds: any admissible rule minimizes Bayes risk under some prior.

---

## Bayes versus Frequentist Estimation

Bayesian methods construct the conditional probability $\mu(\theta \mid d)$ — a description of the uncertainty associated with unknowns $\theta$ given knowns $d$. Frequentist methods, in particular maximum likelihood, construct estimates $\hat{\theta}$ based on $\mu(d \mid \theta)$, in essence describing the ex ante uncertainty of knowns given unknowns.

For the decision problems described here, the Bayesian calculation is the policy-relevant one. Why then are frequentist methods so common? One reason is discomfort over priors (though, as discussed, there are ways to address the absence of principled prior choice). Another is computational: Bayesian methods can be difficult to implement, and the growth of computational power is a reason for their increasing popularity.

---

## Fiducial Inference

The difference between the approaches: Bayesian methods construct probabilities of unobservables based on observables, whereas frequentist methods analyze probabilities of observables based on unobservables. **Fiducial inference**, originally proposed by Ronald Fisher, attempts to bridge this distinction.

Consider a sequence of $N(\lambda, \sigma^2)$ random variables $x_i$, $i = 1, \ldots, K$, with variance known. The sample mean $\bar{x}$ has density $N\left(\lambda, \frac{\sigma^2}{K}\right)$. It is straightforward to show that $\bar{x} - \lambda$ has density $N\left(0, \frac{\sigma^2}{K}\right)$ — which is also the density of $\lambda - \bar{x}$. The quantity $\lambda - \bar{x}$ is a **pivotal quantity**: its distribution does not depend on $\lambda$. The fiducial argument is that since the deviation of $\lambda$ from $\bar{x}$ is defined by a pivotal quantity, this pivotal quantity characterizes the uncertainty associated with the parameter:

$$
\bar{x} - \lambda \sim N\left(0, \frac{\sigma^2}{K}\right) \;\; \text{implies} \;\; \lambda \sim N\left(\bar{x}, \frac{\sigma^2}{K}\right) \tag{2.16}
$$

The argument is generally rejected because, without justification, it transforms a nonrandom object ($\lambda$) into a random one. A variation of Fisher's ideas, **structural inference**, has been developed by Donald Fraser (see Fraser 1968); it is also controversial and has yet to impact empirical work, though unlike Fisher's original formulation it cannot be regarded as a failure.

---

## Some Final Thoughts

For policymakers, these approaches should be regarded as heuristics whose appropriateness cannot be disentangled from context: one might mix ambiguity aversion into a vaccine program decision differently than into a decision about nuclear weapons. In the 1950s John von Neumann argued for a preemptive US nuclear strike on the Soviet Union, on the grounds that the USSR might strike first or become too powerful later — a position justifiable as either minimax or minimax regret. The absurdity of the position holds at two levels: the moral level (one cannot engage in such analyses without taking normative considerations seriously), and the problem of assigning probability 1 to the least desirable events when they have unknown but negligible probabilities.

---

## Difficulty checklist

- **Posterior $\propto$ likelihood $\times$ prior:** $\mu(\theta \mid d) \propto \mu(d \mid \theta)\mu(\theta)$ works because $\mu(d)$ is just a normalizing constant; forgetting the normalization direction (which factor is the likelihood) is a common slip.
- **Quadratic loss $\Rightarrow$ posterior mean:** under $\ell(\theta,d,c) = (\theta - c(d))^2$, the optimal rule is $c(d) = \mathbb{E}(\theta \mid d)$; other loss functions (absolute error, asymmetric) give other functionals (median, quantiles).
- **Improper priors:** the uniform prior on $(-\infty, \infty)$ does not integrate to 1 — improper priors can be usable, but one must verify the resulting posterior is proper.
- **$\mu(\sigma) \propto 1/\sigma$ motivation:** it corresponds to a uniform prior on $\log \sigma$, not on $\sigma$ itself; keep track of which parameterization the uniformity claim applies to.
- **Uniform prior is not transformation-invariant:** by (2.5), nonlinear reparameterization of a uniform prior is generally non-uniform — "ignorance about $\theta$" does not automatically give "ignorance about $f(\theta)$".
- **Jeffreys prior trade-off:** $\mu(\theta) \propto I(\theta)^{1/2}$ restores invariance but makes the prior depend on the likelihood — conceptually odd for an object meant to encode beliefs about parameters alone.
- **Normal posterior weights:** in (2.8) the posterior mean weights $\kappa_0$ and $\bar{x}$ by *precisions* (inverse variances); as $K \to \infty$ the prior washes out. Do not weight by variances directly.
- **Minimax $\neq$ minimax regret:** minimax looks at absolute loss levels, regret at differences from the best feasible action per state; the two criteria can rank policies differently (Example 2.4).
- **Minimax regret violates IIA:** adding an irrelevant (even a terrible) alternative can flip the regret ranking (Examples 2.5–2.6); the axiomatization repair is "independence of never-optimal alternatives" (Stoye).
- **Rawls as minimax:** the difference principle is the $\alpha \to \infty$ limit of $W = \left( \sum_i \omega_i^{-\alpha} \right)^{-1/\alpha}$ — infinite risk aversion, not a separate ethical primitive.
- **Risk vs. Bayes risk:** the risk function $R(\theta, c(D))$ conditions on $\theta$ and averages over data; Bayes risk additionally averages over $\theta$ under the prior. Admissibility links them: Bayes-risk minimizers under proper priors are admissible, and every admissible rule is Bayes under some prior.
- **Bayes vs. frequentist direction:** Bayesian = probabilities of unknowns given observables, $\mu(\theta \mid d)$; frequentist = probabilities of observables given unknowns, $\mu(d \mid \theta)$. Decision problems call for the former.
- **Fiducial inference flaw:** treating a pivotal quantity's distribution as a distribution for the fixed parameter $\lambda$ converts a nonrandom object into a random one without justification — the standard reason the argument is rejected.
