# Lecture 1 — Probability Theory

*PPHA 42000 Applied Econometrics I, Steven N. Durlauf, Fall 2026*

## Why Probability Theory

Economic data, ranging from macroeconomic aggregates such as gross domestic product to microeconomic outcomes such as a binary indicator for high school graduation, are realizations of random variables. A random variable is a mapping from an underlying probabilistic mechanism into real or complex numbers. A stochastic process is a collection of random variables associated with index sets.

**Example 1.1. Coin tosses and payments.** A fair coin is tossed; a person receives 1 dollar on heads and loses 1 dollar on tails. Underlying this simple description is a combination of the toss itself and the translation of the toss into a numerical value.

The set of possible outcomes is the sample space:

$$
\Omega = \{H, T\}
$$

Individual elements $\omega \in \Omega$ are the possible realizations of the process. Probabilities are defined on subsets of $\Omega$; the mapping $\mu(\cdot)$ from these subsets to $[0,1]$ is called a probability measure. For this simple case, probabilities may be defined on all subsets — this is not always possible in more complicated environments. The collection of subsets on which probabilities are defined is denoted $\mathcal{F}$.

Since the coin is fair:

$$
\mu(\emptyset) = 0, \quad \mu(\Omega) = 1, \quad \mu(H) = 0.5, \quad \mu(T) = 0.5 \tag{1.1}
$$

Let $x(\omega)$ denote a map from $\Omega$ to $\{-1, 1\}$. Probabilities of the values of $x$ are induced by $(\Omega, \mathcal{F}, \mu)$:

$$
\mu(x = -1) = 0.5, \quad \mu(x = 1) = 0.5 \tag{1.2}
$$

The map $x$ is a random variable.

---

## Probability Spaces

A formal definition of probability requires specifying the events from which a realization is drawn, the sets of events over which probabilities may be defined, and the probabilities themselves. The object of study is $(\Omega, \mathcal{F}, \mu)$, called a probability space, where:

- $\Omega$ is the sample space; elements $\omega$ are realizations of the uncertainty with respect to which we assign probabilities.
- $\mathcal{F}$ is a $\sigma$-algebra of sets $f$, each a subset of $\Omega$; these are the objects to which probabilities are assigned.
- $\mu$ is a probability measure, assigning probabilities to elements of $\mathcal{F}$.

**Algebra of sets:** a collection of subsets of $\Omega$ such that (i) $\Omega$ is in the collection, and (ii) the collection is closed under finite complements and unions. A **$\sigma$-algebra** is an algebra of sets that is also closed under countable unions. Closure under countable unions ensures that limits of sequences of sets are included in $\mathcal{F}$.

**Probability measure:** a mapping $\mu(\cdot)$ from $\mathcal{F}$ to $[0,1]$ such that:

1. $\mu(f) \geq 0 \ \forall f \in \mathcal{F}$ (probabilities cannot be negative).
2. $\mu(\Omega) = 1$ (the probability that something happens is 1).
3. For any countable collection of disjoint sets, $\mu\left(\bigcup_i f_i\right) = \sum_i \mu(f_i)$ (countable additivity: for mutually exclusive events, the probability that either happens is the sum of the individual probabilities).

Probability densities are special cases of probability measures, as are discrete probabilities; densities derive from measures and need not always exist, just as not every function is differentiable.

**Why the $\sigma$-algebra?** For probabilities to be defined coherently in the sense made precise by the definition of a probability measure, assignment of probabilities may not be possible for all subsets of $\Omega$. Restricting probabilities to a $\sigma$-algebra of subsets avoids these technical problems: in general it is impossible to define probabilities over all possible subsets while preserving the required properties.

**Borel $\sigma$-algebra:** the standard example, denoted for $\mathbb{R}$, is the smallest $\sigma$-algebra of subsets of $\mathbb{R}$ that contains all half-open intervals of the form $(a, b]$. One can conceptualize it as generated around the half-open intervals by following the inclusion rules that define a $\sigma$-algebra. It is a very rich collection.

---

## Random Variables

For a probability space $(\Omega, \mathcal{F}, \mu)$, suppose there exists a mapping $x$ from $\Omega$ to $\mathbb{R}$ such that for every member $B$ of the Borel sets of $\mathbb{R}$, the set defined by $\{\omega : x(\omega) \in B\}$ is an element of $\mathcal{F}$. Then $x(\omega)$ is a **random variable** (the measurability requirement). The probability space $(\Omega, \mathcal{F}, \mu)$ determines the probabilities associated with $x$.

This is a natural way to think about socioeconomic data. Behaviors such as high school completion are coded into numbers; in other cases, such as wages, the outcome is itself a real number whose randomness is induced by the randomness of underlying factors.

---

## Stochastic Processes

A stochastic process is a collection of indexed random variables, with an index $a \in A$ distinguishing the individual random variables.

**Time series:** a collection of random variables $\{x_t\}$, $t \in \mathbb{Z}$. Natural representations for macroeconomic data (e.g., monthly unemployment). A key feature: the index is interpretable in terms of distance — it is meaningful to say that for annual observations, $x_{1970}$ is closer to $x_{1971}$ than to $x_{2015}$.

**Strict stationarity:** the time series $\{x_t\}$ is strictly stationary if for all $\tau$ and any collection of indices $t_1, \ldots, t_k$:

$$
\mu(x_{t_1}, \ldots, x_{t_k}) = \mu(x_{t_1+\tau}, \ldots, x_{t_k+\tau}) \tag{1.3}
$$

**Second-order stationarity:** $\mathbb{E}(x_t) = \mathbb{E}(x_{t+\tau})$ and $\mathbb{E}(x_t x_{t+\tau})$ is determined by $\tau$, the difference between the locations of the random variables in time. Second-order stationarity is typically sufficient for linear time series models, since first and second moments are the objects of interest.

**Violations of stationarity** take many forms. Heteroskedasticity of errors in a linear regression can induce a stationarity violation, but the interpretation of the regression coefficients is not affected per se. In other cases parameters of interest may not be constant — consider time series on the Russian economy from January 1992 to 2000, or the Chinese economy before versus after 1980. These data exhibit nonstationarity driven by the emergence of a new type of economy; they cannot be well studied with constant-parameter models, nor does standard econometric methodology apply.

**Cross-sections:** a second common example. For a group of high school seniors enrolled in the same school in the same year, the outcomes graduate/not-graduate are a collection of random variables $x_i$ indexed by student $i$. Unlike the time series case, the index has no interpretation: designating one student "7" and another "8" does not make them closer than students "48" and "100". Panel data combine the two, giving $x_{it}$.

---

## Construction of Time Series via a Measure-Preserving Transformation

One can conceptualize time series as generated by a common probability space, using a measure-preserving transformation for strictly stationary series.

**Definition 1.2. Measure-preserving transformation.** $T$ is a mapping from $\Omega$ to $\Omega$ such that for all $f \in \mathcal{F}$:

$$
\mu(T^{-1} f) = \mu(f)
$$

One can use $T$ to define a sequence of random variables $x_t = x(T^t \omega)$. Such a sequence is strictly stationary. To prove this, recall that strict stationarity requires, for all $a, b$:

$$
\mu(\omega : a \leq x_t \leq b) = \mu(\omega : a \leq x(T^t \omega) \leq b) \tag{1.4}
$$

Let $C$ denote this set of $\omega$'s. Since $T$ is measure preserving, for $D = TC$:

$$
\mu(C) = \mu(D) \tag{1.5}
$$

However,

$$
\mu(D) = \mu(\omega : a \leq x(T^{t+1}\omega) \leq b) = \mu(\omega : a \leq x_{t+1} \leq b) \tag{1.6}
$$

which is the desired result.

---

## Ergodicity

The properties of the measure-preserving transformation determine much about what may be learned about the probability structure of a time series from a given sample path realization. The key property is ergodicity, which requires the notion of an invariant set.

**Definition 1.3. Invariant set.** A set $S$ is invariant under a measure-preserving transformation $T$ if $TS = S$. A set $S$ is almost invariant if $\mu(S \triangle TS) = 0$, where $A \triangle B = A \cup B - A \cap B$ is the symmetric difference.

**Definition 1.4. Ergodicity.** A measure-preserving transformation is ergodic if for every invariant and/or almost-invariant set $S$, either $\mu(S) = 0$ or $\mu(\Omega - S) = 0$. In words: the only subsets of $\Omega$ in which the system can be trapped (never leave) are those with probability 1 and probability 0.

**Theorem 1.1.** If $x(\omega)$ is a random variable with finite expectation, then for each $\omega$, outside a set of measure 0, there exists a map $\hat{x}(\omega)$ such that:

$$
\lim_{n \to \infty} \frac{1}{n} \sum_{j=0}^{n-1} x(T^j \omega) = \hat{x}(\omega) \tag{1.7}
$$

and $\hat{x}(\omega)$ is finite with probability 1.

**Theorem 1.2.** If $T$ is ergodic, then:

$$
\hat{x}(\omega) = \mathbb{E}(x) \quad (\mu = 1) \tag{1.8}
$$

Together, these theorems say that except for a set of $\omega$'s with probability 0, the average of every sample path realization converges to the expected value of $x_t$. The ergodic theorem is one way to understand the conditions for a law of large numbers to hold. The intuition: the random process generating $x_t$ cannot get "stuck" in a subset of $\Omega$ with probability less than 1; if it could, the data would not reveal values of $x_t$ that could have occurred (with positive probability) but do not.

**Contrast of two examples.** Consider the time series

$$
x_t = \varepsilon_t \tag{1.9}
$$

where $\varepsilon_t$ is a sequence of independent $N(0,1)$ random variables. The sample average $\bar{x}_T = \frac{1}{T} \sum_{t=1}^{T} x_t$ is $N(0, T^{-1})$ and so converges with probability 1 to the common expected value 0 as $T$ grows. In contrast, suppose

$$
x_t = \varepsilon \tag{1.10}
$$

where $\varepsilon$ is a single $N(0, \sigma^2)$ random variable. Then $\bar{x}_T = \varepsilon \sim N(0, \sigma^2)$: there is no convergence to the expected value 0, since the different possible values are never seen — only the single realization of the random variable.

**Empirical implications.** It is common to introduce dummy variables to account for violations of second-order stationarity (e.g., a 9/11 dummy); the ergodic theorem indicates why the effects of such one-time events cannot be consistently estimated. More broadly, one must consider whether one studies an ergodic or nonergodic environment. Questions about the long-run effects of unique historical events — the intensity of serfdom on regional development in Russia, the Great Leap Forward on current Chinese outcomes, the causes of World War I — involve one-time events and counterfactuals for which there are no data. This is a textbook case of nonergodicity: the causal inference literature, whose deep idea is to identify counterfactuals with observed data, is essentially mute here. Such questions remain meaningful; they require theory and qualitative information, and this does not render the answers noncredible. (With apologies to Hamlet: there are more things in heaven and earth than are dreamt of in statistical visions of the world.)

---

## Exchangeability

### Definitions

To motivate, consider the standard linear model:

$$
y_i = x_i \beta + \varepsilon_i \tag{1.11}
$$

It is common to assume the errors are independent and identically distributed (i.i.d.). Such assumptions do not naturally correspond to the substantive social science knowledge one brings to an empirical exercise; exchangeability does. The definition employs the permutation operator $\rho(\cdot)$, which rearranges a set of integers one-to-one (every index appears exactly once in the permuted set).

**Definition 1.5. Exchangeability.** A collection of random variables $\varepsilon_i$ is exchangeable if for every finite subset $\varepsilon_{i_1}, \ldots, \varepsilon_{i_N}$ and every permutation operator $\rho(\cdot)$:

$$
\mu(\varepsilon_{i_1} \leq a_1, \ldots, \varepsilon_{i_N} \leq a_N) = \mu(\varepsilon_{\rho(i_1)} \leq a_1, \ldots, \varepsilon_{\rho(i_N)} \leq a_N) \tag{1.12}
$$

Exchangeability means there are symmetries in the joint probabilities: the subscripts identifying particular observations do not matter. Clearly, if the collection $\varepsilon_i$ is i.i.d., it is exchangeable. Exchangeability, however, does not imply independence, as is easily seen for the case

$$
\varepsilon_i = \eta_i + \xi \tag{1.13}
$$

where $\eta_i$ is i.i.d. across $i$ and $\xi$ is a single random variable independent of $\eta_i$ for all $i$.

**Definition 1.6. Conditional exchangeability.** Suppose each $\varepsilon_i$ is associated with a vector $X_i$ (whose elements may be realizations of a stochastic process). The collection $\varepsilon_i$ is conditionally exchangeable if for any finite collection $\varepsilon_{i_1}, \ldots, \varepsilon_{i_N}$ and every permutation operator $\rho(\cdot)$:

$$
\mu(\varepsilon_{i_1} \leq a_1, \ldots, \varepsilon_{i_N} \leq a_N \mid X_{i_1} = X, \ldots, X_{i_N} = X) = \mu(\varepsilon_{\rho(i_1)} \leq a_1, \ldots, \varepsilon_{\rho(i_N)} \leq a_N \mid X_{i_1} = X, \ldots, X_{i_N} = X) \tag{1.14}
$$

One can modify this to require that the $X_i$ are not equal but lie in a common set $\mathcal{X}$; conditional exchangeability with respect to this set holds if for every permutation operator $\rho(\cdot)$:

$$
\mu(\varepsilon_{i_1} \leq a_1, \ldots, \varepsilon_{i_N} \leq a_N \mid X_{i_1} \in \mathcal{X}, \ldots, X_{i_N} \in \mathcal{X}) = \mu(\varepsilon_{\rho(i_1)} \leq a_1, \ldots, \varepsilon_{\rho(i_N)} \leq a_N \mid X_{i_1} \in \mathcal{X}, \ldots, X_{i_N} \in \mathcal{X}) \tag{1.15}
$$

This latter definition suggests an important use of conditional exchangeability: partitioning observations into exchangeable classes. For example, if $X_i$ encodes an individual's ethnicity and gender, conditional exchangeability of errors within ethnicity/gender combinations leads naturally to partitioning the data along these lines.

### Why exchangeability?

Exchangeability pushes a researcher to think about model specification in terms of the implied residuals. Omitted variables are easily interpreted as an exchangeability violation: conditional on the omitted variable, the errors in the misspecified regression are no longer exchangeable. Parameter heterogeneity also leads to nonexchangeability.

Brock and Durlauf (2001) argue that exchangeability links substantive social science knowledge to error structure: this knowledge can be used to evaluate the plausibility of exchangeability. Good empirical practice is to question whether the errors in a model are exchangeable and, if not, whether the violation invalidates the purposes for which the regression is used. This cannot be done algorithmically; as with empirical work generally, it requires judgment by the analyst. (See also Draper et al. 1993 on the primacy of exchangeability in statistical model building, and McCullagh 2005.)

Exchangeability is not necessary for a statistical model to be interpretable as a behavioral structure. Heteroskedasticity violates exchangeability but does not affect the interpretability of regression coefficients per se. Rather, exchangeability is a criterion by which a researcher can evaluate modeling choices.

### DeFinetti's theorem

Remarkably, exchangeability also provides a justification for the assumption that data are i.i.d. The first formulation of the result (for binary random variables) is due to Bruno de Finetti.

**Theorem 1.3. DeFinetti's Theorem.** Suppose $\varepsilon_i$ is an infinite real-valued exchangeable sequence. Then there exists a probability measure $Q$ on the space of scalar random variable probability measures $p$ such that:

$$
\mu\left(\bigcup_i \varepsilon_i\right) = \int \prod_i p(\varepsilon_i) \, dQ(p) \tag{1.16}
$$

In words: the probability measure describing any infinite exchangeable sequence can be written as a mixture of i.i.d. probability measures. Each sample path realization obeys one of the probability measures, so each sample path behaves as an i.i.d. sequence. Note the theorem does not apply to finite exchangeable sequences — the standard counterexample is draws without replacement from an urn with finitely many red and black balls: the color sequence is exchangeable, but the probability of the next draw depends on what has already been removed. And even though DeFinetti's theorem justifies treating data as i.i.d., this is conditional on the realization of $p$, so the sequence may not contain information allowing $Q$ to be recovered.

---

## Difficulty checklist

- **$\sigma$-algebra vs. algebra:** an algebra requires closure under finite complements and unions; a $\sigma$-algebra additionally requires closure under countable unions. Confusing the two misses why limits of sets are measurable.
- **Why $\mathcal{F}$ may not be all subsets:** probabilities cannot in general be coherently assigned to every subset of $\Omega$; the $\sigma$-algebra restriction is a necessity, not a convenience.
- **Borel $\sigma$-algebra generation:** it is the *smallest* $\sigma$-algebra containing the half-open intervals $(a, b]$ — not the set of intervals itself.
- **Random variable = measurable map:** $x : \Omega \to \mathbb{R}$ requires $\{\omega : x(\omega) \in B\} \in \mathcal{F}$ for every Borel set $B$; the randomness lives in $\omega$, not in $x$.
- **Strict vs. second-order stationarity:** strict stationarity invaries the entire joint distribution under time shifts; second-order stationarity only fixes the mean and makes autocovariances depend on the gap $\tau$ alone.
- **Heteroskedasticity is a stationarity violation** but does not by itself affect the interpretation of regression coefficients — do not conflate statistical violation with interpretive failure.
- **Measure-preserving transformation direction:** the defining property is $\mu(T^{-1}f) = \mu(f)$, on the *preimage*, not $\mu(Tf) = \mu(f)$.
- **Ergodicity trap intuition:** ergodicity means the only invariant sets have measure 0 or 1; if a process could be trapped in a subset with intermediate probability, sample paths would not reveal the full distribution and the LLN fails.
- **$x_t = \varepsilon_t$ vs. $x_t = \varepsilon$:** with a fresh draw each period the sample mean converges; with one permanent draw it does not. This is the canonical ergodic/nonergodic contrast.
- **One-time events:** dummy variables for unique events (9/11) cannot be consistently estimated — the ergodic theorem is the formal reason.
- **i.i.d. $\Rightarrow$ exchangeable, but not conversely:** $\varepsilon_i = \eta_i + \xi$ with shared $\xi$ is exchangeable but dependent. Reversing the implication is a classic error.
- **DeFinetti conditions:** the theorem requires an *infinite* exchangeable sequence; finite sequences (urn draws without replacement) need not be mixtures of i.i.d. measures, and the mixing measure $Q$ may not be recoverable from a single realized path.
