# Lecture 5 — Linear Regression

*PPHA 42000 Applied Econometrics I, Steven N. Durlauf, Fall 2026*

## Population-Level Interpretation

These notes give a **population-level** interpretation of linear regression: they describe linear relationships between random variables themselves, not between sample realizations of random variables. The question is what a linear projection of one random variable onto others means, and when the projection coefficients coincide with structural (behavioral) parameters.

---

## The Linear Regression Model and the Orthogonality Condition

The canonical linear regression model is

$$
y = x\beta + \varepsilon
$$

where $y$ is a scalar random variable, $x$ is a $1 \times K$ vector of random variables, $\beta$ is a $K \times 1$ vector of coefficients, and $\varepsilon$ is a scalar random variable. All random variables are assumed to be zero mean (relaxable without loss of generality). The important substantive assumption is the **orthogonality condition**

$$
\forall k, \quad \mathbb{E}(x_k \varepsilon) = 0
$$

This is the standard orthogonality (exogeneity) condition for linear models; it is an assumption about the world, not a property of a procedure.

---

## Recovering $\beta$ from Population Moments

From the perspective of empirical work, the variance–covariance matrix of realizations of $(y, x)$ constitutes the statistics from which estimates of $\beta$ are constructed. Premultiplying the model by $x^\top$ gives

$$
x^\top y = x^\top x\, \beta + x^\top \varepsilon
$$

Taking expectations of both sides (and applying the orthogonality condition) yields

$$
\mathbb{E}(x^\top y) = \mathbb{E}(x^\top x)\, \beta \quad \Longleftrightarrow \quad \text{cov}(x^\top, y) = \text{var}(x)\, \beta
$$

If the $K \times K$ matrix $\text{var}(x)$ is invertible, i.e. $\text{var}(x)^{-1}$ exists, then

$$
\beta = \text{var}(x)^{-1}\, \text{cov}(x^\top, y)
$$

which is the **population version of the OLS formula**.

**Interpretation of invertibility.** Nonexistence of $\text{var}(x)^{-1}$ means there exists a vector $c$ such that $c^\top \text{var}(x)\, c = 0$, i.e. $\text{var}(xc) = 0$ — some linear combination of the regressors is degenerate. So invertibility of $\text{var}(x)$ is exactly the requirement that the components of $x$ be **linearly independent**: none can be expressed as a linear combination of the others.

---

## Best Linear Prediction as Hilbert-Space Projection

Contrast the identification discussion above with the following minimization problem:

$$
\min_{\pi \in \mathbb{R}^K}\ \text{var}(y - x\pi)
$$

This is the **best linear predictor** (linear projection) problem, and the Hilbert space formalism of Lecture 4 is one route to solving it. Form the minimal Hilbert space around the random variables $(y, x)$ with metric $\|\cdot\| = \text{var}(\cdot)^{1/2}$ and inner product

$$
\langle u, v \rangle = \text{cov}(u, v)
$$

The construction starts with the $K+1$ random variables, all their weighted linear combinations, and the limits thereof — a Hilbert space of random variables denoted $\mathcal{H}_{y,x}$. Analogously, construct the Hilbert space around the elements of $x$ alone, denoted $\mathcal{H}_x$.

Since $\mathcal{H}_x \subseteq \mathcal{H}_{y,x}$, the Hilbert-space decomposition theorem gives

$$
\mathcal{H}_{y,x} = \mathcal{H}_x \oplus \mathcal{H}_\perp
$$

where $\mathcal{H}_x$ and $\mathcal{H}_\perp$ are orthogonal. Because $\mathcal{H}_{y,x}$ is built around $K+1$ random variables, there exist an element $\hat{y} \in \mathcal{H}_x$ (the projection of $y$ onto $\mathcal{H}_x$) and an element $\eta \in \mathcal{H}_\perp$ such that

$$
y = \hat{y} + \eta
$$

The projection theorem says this decomposition solves the minimization problem, and except for pathological cases $\hat{y} = x\pi$, which characterizes the solution: the **best linear predictor** of $y$ given $x$.

---

## Key Subtleties of the Projection

**Uniqueness under linear dependence.** The solution $\hat{y} = x\pi$ is **unique even if the elements of $x$ are linearly dependent** — even though the coefficient vector $\pi$ then is not. The data can place restrictions on $\pi$ even under linear dependence: the fact that one cannot find a unique combination of the $x$'s closest to $y$ does not mean all combinations are equally close. This is a form of **partial identification** (of $\pi$, not of $\hat{y}$).

**Projection coefficients vs. behavioral parameters.** Does the solution $\pi$ equal $\beta$ in the structural model? **Yes, if the orthogonality condition holds** — this is a sufficient condition (why is it not necessary?). The point is important: suppose the regression equation is a behavioral relationship derived from economic theory. If the theory does not imply orthogonality, then projecting $y$ onto $x$ — which always solves the best-linear-prediction problem — **does not produce the behavioral parameters**. It produces the parameters that render the regressors orthogonal to the component of $y$ lying in $\mathcal{H}_\perp$. A regression is always a projection; whether the projection recovers structure depends on assumptions about the world.

---

## Omitted Variables

Assume the structural model and orthogonality hold, but partition $x$ into $x_1$ and $x_2$ and rewrite the model as

$$
y = x_1 \beta_1 + x_2 \beta_2 + \varepsilon
$$

Suppose one omits $x_2$ and projects $y$ onto $x_1$ alone, i.e. solves

$$
\min_{\pi \in \mathbb{R}^{K_1}}\ \text{var}(y - x_1 \pi)
$$

What will $\pi$ equal? Assume for simplicity that $x_1$ and $x_2$ are scalars. By the Hilbert-space projection theorem, $x_2$ can itself be written as a projection on $x_1$ plus an orthogonal residual:

$$
x_2 = x_1 \gamma + \nu, \qquad \text{cov}(x_1, \nu) = 0
$$

Substituting into the structural model,

$$
y = x_1 \beta_1 + (x_1 \gamma + \nu)\beta_2 + \varepsilon = x_1(\beta_1 + \gamma \beta_2) + \nu \beta_2 + \varepsilon
$$

Since $\text{cov}(x_1, \nu\beta_2 + \varepsilon) = 0$, the short-projection coefficient is

$$
\pi = \beta_1 + \gamma \beta_2
$$

Hence $\pi \neq \beta_1$ unless $\gamma\beta_2 = 0$. The product $\gamma\beta_2$ vanishes if $\gamma = 0$ — the included and excluded variables are uncorrelated, $\text{cov}(x_1, x_2) = 0$ — and/or if $\beta_2 = 0$ — the omitted variable is irrelevant. So $\pi \neq \beta_1$ exactly when the omitted variable **both matters in the structural equation and covaries with the included regressor**. For the multivariate case, the analogy is that $\gamma$ is a matrix and $\nu$ a vector.

**Signing the bias.** It is easy to attach a sign to $\gamma\beta_2$ when $x_1$ and $x_2$ are scalars — the sign is determined by the signs of $\gamma$ and $\beta_2$. This is **not** possible in multivariate contexts.

---

## Measurement Error

**Errors in the regressor: attenuation.** Suppose the relationship of interest is

$$
y = x\beta + \varepsilon, \qquad \text{cov}(x, \varepsilon) = 0, \quad \mathbb{E}(x) = \mathbb{E}(\varepsilon) = 0
$$

but $x$ is measured with error and therefore unobservable. The available data are $y$ and

$$
x^* = x + \eta, \qquad \text{cov}(x, \eta) = \text{cov}(\eta, \varepsilon) = 0
$$

(classical measurement error). What does the regression on the mismeasured regressor deliver? Write

$$
y = x^* \pi + \varsigma
$$

and analyze it via the auxiliary regression $y = x\beta_1 + x^*\beta_2 + \varepsilon$, where the true model implies $\beta_1 = \beta$ and $\beta_2 = 0$. From the projection

$$
x = x^*\gamma + \xi
$$

it is immediate that

$$
\gamma = \frac{\text{var}(x)}{\text{var}(x) + \text{var}(\eta)}
$$

Therefore the regression coefficient on the mismeasured regressor obeys

$$
\pi = \frac{\text{var}(x)}{\text{var}(x) + \text{var}(\eta)}\, \beta < \beta
$$

This is the classic **attenuation** (errors-in-variables) result: in a bivariate regression, measurement error in the regressor that is orthogonal to everything else shrinks the regression coefficient toward $0$ relative to the true coefficient. The multivariate case is more complicated and admits no such clean proposition.

**Errors in the dependent variable: no bias.** In contrast, suppose only the dependent variable is mismeasured: one observes

$$
y^* = y + \nu
$$

Then the relationship between observables is

$$
y^* = x\beta + \varepsilon + \nu
$$

If $\text{cov}(x, \nu) = 0$, the linear projection of $y^*$ onto $x$,

$$
y^* = x\pi + \omega
$$

produces $\pi = \beta$, since $\text{cov}(x, \varepsilon + \nu) = 0$. **Measurement error in the dependent variable causes no bias** (it only inflates the error variance) — an asymmetry with the regressor case worth remembering.

---

## Difficulty checklist

- **Population vs. sample.** Everything here is about random variables, not data: $\beta = \text{var}(x)^{-1}\text{cov}(x^\top, y)$ is a statement about population moments. Sample OLS is the analogue, not the definition.
- **The orthogonality condition is an assumption, not a theorem.** $\mathbb{E}(x_k\varepsilon) = 0$ for all $k$ must come from economic reasoning about the behavioral model; a projection *constructs* orthogonality between regressors and the residual, which is not the same thing.
- **Every regression is a projection; not every projection is structure.** Projecting $y$ onto $x$ always solves $\min_\pi \text{var}(y - x\pi)$ and always yields $\hat{y} = x\pi$ with $y - \hat{y} \perp \mathcal{H}_x$. The coefficients equal behavioral $\beta$ only if the orthogonality condition holds.
- **Sufficiency, not necessity, of orthogonality.** $\mathbb{E}(x\varepsilon) = 0$ guarantees $\pi = \beta$, but $\pi = \beta$ can occur in special configurations without full orthogonality — do not confuse the directions.
- **Uniqueness of $\hat{y}$ vs. uniqueness of $\pi$.** Under perfect multicollinearity, the fitted value $\hat{y} = x\pi$ is still unique even though $\pi$ is not. Restrictions on $\pi$ survive linear dependence — a partial-identification phenomenon.
- **Invertibility of $\text{var}(x)$ = linear independence of regressors.** $\exists c: c^\top\text{var}(x)c = 0 \Leftrightarrow \text{var}(xc) = 0 \Leftrightarrow$ a degenerate linear combination exists. This is the population version of the full-rank condition.
- **OVB formula and its two channels.** $\pi = \beta_1 + \gamma\beta_2$: bias requires *both* $\gamma \neq 0$ (omitted covaries with included) and $\beta_2 \neq 0$ (omitted matters structurally). Sign is determinable only in the scalar case.
- **Projection of the omitted variable.** The derivation runs through $x_2 = x_1\gamma + \nu$ with $\text{cov}(x_1, \nu) = 0$ — the omitted variable is itself decomposed by projection, which is why the residual $\nu\beta_2$ is innocuous in the short regression.
- **Measurement error is asymmetric.** Error in the regressor attenuates: $\pi = \frac{\text{var}(x)}{\text{var}(x) + \text{var}(\eta)}\beta$ (bivariate case). Error in the dependent variable with $\text{cov}(x,\nu)=0$ leaves $\pi = \beta$ unbiased. Do not apply attenuation logic to $y$-side error.
- **Attenuation factor is a variance ratio.** $\gamma = \text{var}(x)/(\text{var}(x) + \text{var}(\eta)) \in (0,1]$; it equals $1$ only when $\text{var}(\eta) = 0$. Larger measurement-error variance $\Rightarrow$ more shrinkage toward zero.
