# Numerical Integration of ODEs

Notes from the online course at MIT OpenCourseWare

## First order ODEs

Let us consider the general first order ODE of the unknown variable $u(t)$:

$$
    u_t = f(u(t), t) \qquad u(0) = u_0 \qquad 0 < t < T
$$

The time interval $[0, T]$ is discretized with a timestep $\Delta t$ and the numerical approximation of $u(n\Delta t)$ is denoted $v^n$. The numerical scheme is denoted by:

$$
    v^{n + 1} = N(v^i, t^i, \Delta t)
$$

where $i$ spans the set of all indices required by the numerical scheme.

In the case of a $s$-multistep method

$$
    v^{n + 1} + \sum_{i=1}^s \alpha_i v^{n + 1 - i} = \sum_{i = 0}^s \beta_i f^{n + 1 - i}
$$

### Local order of accuracy and consistency

The local truncation error (LTE) is defined as

$$
    \tau = N(v^i, t^i, \Delta t) - u^{n + 1}
$$

The local order of accuracy $p$ is

$$
    | \tau | = \mc{O}(\Delta t^{p + 1})
$$

The global error $e(T) = u(T) - v^{T/\Delta t}$ is expected to behave as

$$
\begin{aligned}
    &e(T) = \sum_{n = 1}^{T/\Delta t} \Delta e^n \\
    \implies &e(T) = \sum_n \mc{O}(\Delta t^{p + 1}) \\
    \implies &e(T) = \mc{O}(\Delta t^p)
\end{aligned}
$$

Consistency requires that the scheme be at least first order ($p = 1$) for the numerical method to capture the model.

### Convergence and global order of accuracy

A finite difference method for solving

$$
    u_t = f(u(t), t) \qquad u(0) = u_0 \qquad 0 < t < T
$$

is convergent if

$$
    \max_{n \in \iint{0}{T/\Delta t}} |v^n - u(n\Delta t)| \to 0 \qquad \text{as} \qquad \Delta t \to 0
$$

The method is said to have global order of accuracy $p$ is

$$
    \max_{n \in \iint{0}{T/\Delta t}} |v^n - u(n\Delta t)| \leq \mc{O}(\Delta t^p)
$$

for any $f(u, t)$ that has $p$ continuous derivatives.

### Stability

Let us consider the unforced problem associated to an $s$-multistep method

$$
    v^{n + 1} + \sum_{i=1}^s \alpha_i v^{n + 1 - i} = 0
$$

A multistep method is stable if all solutions to the previous are bounded as $n\to +\infty$

In practice we assume that the solution has the following form

$$
    v^n = v^0 z^n
$$

We solve the unforced problem in $z$, determine the roots of the polynomial and check whether or not $|z| \leq 1$. For a multistep method $|z| = 1$ is stable only if it is a simple root.

### Dahlquist Equivalence Theorem

A stable and consistent $s$-multistep method is convergent (similar to Lax Theorem).

### Systems of ODEs and Eigenvalue stability



<!-- Katex macros definitions -->
$$
    \gdef\mrm#1{\mathrm{#1}}
    \gdef\pdv#1#2{\frac{\partial #1}{\partial #2}}
    \gdef\vb#1{\mathbf{#1}}
    \gdef\veps{\varepsilon}
    \gdef\P{\mathbb{P}}
    \gdef\R{\mathbb{R}}
    \gdef\N{\mathbb{N}}
    \gdef\mc#1{\mathcal{#1}}
    \gdef\iint#1#2{\llbracket #1, #2 \rrbracket}
$$
