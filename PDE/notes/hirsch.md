# Numerical Computation of Internal and External Flows notes

## Chapter 7: Consistency, Stability and Error Analysis of Numerical Schemes

### Basic concepts and definitions

#### Consistency

A scheme is said to be consistent if it tends to the differential equation when time and space steps tend to zero

#### Stability

Stability is the requirement that all errors must remain bounded when the iteration process advances. Denoting by $u_i^n$ the computed solution and $\bar{u}_i^n$ the exact solution of the numerical scheme, we define the error as $\bar{\veps}_i^n = u_i^n - \bar{u}_i^n$ and require:

$$
    \lim_{n\to +\infty} |\bar{\veps}_i^n| \leq K 
$$
at fixed $\Delta t$ with $K$ independent of $n$.

#### Convergence

Convergence is a condition on the numerical solution. When time and space steps tend to zero, the numerical solution $u_i^n$ must converge to the exact solution $\tilde{u}_i^n$ of the differential equation. Defining the error as:

$$
    \tilde{\veps}_i^n = u_i^n - \tilde{u}_i^n
$$

we require that:

$$
    \lim_{\Delta x \to 0, \Delta t \to 0}\max_{n \in \intint{0}{T/\Delta t}, i \in \intint{0}{L_x/\Delta x}} |\tilde{\veps}_i^n| = 0
$$

### Von Neumann analysis

#### Fourier decomposition of the solution

Let us consider a domain of length $L$ in the $x$-axis with $N + 1$ points and spacing $\Delta x$. Reflecting this domain with respect to the origin the minimum and maximum wavelength and wave numbers are:

$$
\begin{aligned}
    &\lambda_\mrm{max} = 2L \implies k_\mrm{min} = \frac{\pi}{L} \\
    &\lambda_\mrm{max} = 2\Delta x \implies k_\mrm{max} = \frac{\pi}{\Delta x} \\
\end{aligned}
$$

Hence
$$
    k_j = \frac{j\pi}{L}, \: j \in \intint{-N}{N}
$$

And we can decompose the solution into spatial modes:

$$
    u_i^n = \sum_{j=-N}^N V_j^n e^{\I k_j x_i} = \sum_{j=-N}^N V_j^n e^{\I i \cdot k_j \Delta x} = \sum_{j=-N}^N V_j^n e^{\I i \cdot \phi_j}
$$

where $\phi_j = k_j \Delta x \in [-\pi, \pi]$ is the phase angle which remaps spatial frequencies into the range $[-\pi, \pi]$.

#### Amplification factor

The Von Neumann stability condition requires that the amplitudes do not grow in time, we define the amplification factor:

$$
    G = \frac{V^{n+1}}{V^n}
$$

and require that

$$
    |G(\phi_j)| \leq 1 \qquad \forall j \in \intint{-N}{N} 
$$

In practice for a numerical scheme we set:

$$
    u_{i+m}^{n+k} = V^{n+k} e^{\I(i+m)\phi}
$$

and find the amplification factor depending on $\phi$.

### Spectral analysis of numerical schemes

For the exact solution of the differential equation $\tu$ let us consider its mode of wave number $k$ where $\tilde{\omega} = \tilde{\omega}(k)$ is the exact dispersion relation:

$$
    \tu_k(x, t) = \hat{V}(k) e^{\I(kx - \tilde{\omega}t)}
$$

where

$$
    \hat{V}(k) = \frac{1}{2L} \int_{-L}^{L} u(x, 0)e^{-\I kx} dx
$$

Considering now the computed solution with an approximate dispersion relation $\omega = \omega(k)$:

$$
    (u_i^n)_k = \hat{V}(k) e^{\I kx_i} e^{-\I \omega t^n}
$$

We can then identify the exact and numerical amplification factors:

$$
\begin{aligned}
    \tilde{G} &= e^{-\I \tilde{\omega} \Delta  t} = |\tilde{{G}}| e^{-\I\tilde{\Phi}} \\
    G &= e^{-\I \omega \Delta  t} = |\tilde{{G}}| e^{-\I\Phi} 
\end{aligned}
$$

From these two amplification factors the diffusion and dispersion errors are:

$$
\begin{aligned}
    \veps_D &= \frac{|G|}{|\tilde{G}|} \\
    \veps_\phi &= \frac{\tilde{\Phi}}{\Phi}
\end{aligned}
$$

## Chapter 8: General properties and High-Resolution Numerical Schemes

### Two level explicit schemes

We consider a general two level explicit scheme of the general form:

$$
    u_i^{n + 1} = \sum_j b_j u_{i+j}^n
$$

The range of $j$ is called the support of the schemes and is separated between $ju$ upwind points and $jd$ downwind points so that the total number of support points is $M=ju + jd + 1$.

We recall the first order upwind scheme and second order Lax-Wendroff schemes for the linear convection equation:

$$
\begin{aligned}
    &\mrm{FOU} \qquad u_i^{n+1} = u_i^n - \sig(u_i^n - u_{i-1}^n) \\
    &\mrm{LW} \qquad u_i^{n+1} = u_i^n - \frac{\sig}{2}(u_{i+1}^n - u_{i-1}^n) + \frac{\sig^2}{2}(u_{i+1}^n + u_{i-1}^n - 2u_i^n)
\end{aligned}
$$

where $(b_{-1}, b_0) = (\sig, 1 - \sig)$ for the FOU and $(b_{-1}, b_0, b_1) = (\sig(\sig + 1)/2, 1 - \sig^2, \sig(\sig - 1)/2)$ for LW.

We will now search for a general formulation of two level explicit schemes. Starting from Taylor expansions in time and space:

$$
\begin{aligned}
    u_i^{n+1} &= \sum_{m = 0}^{+\infty} \frac{\Delta t^m}{m!} \left(\pmdv{u}{t}{m}\right) \\
    u_{i+j}^n &= \sum_{m = 0}^{+\infty} \frac{(j \Delta x)^m}{m!} \left(\pmdv{u}{x}{m}\right) \\
\end{aligned}
$$

The first requirement is that a constant function should be solution of the problem:

$$
    \sum_j b_j = 1
$$

Reinjecting the Taylor expansions into the initial formula:

$$
    \Delta t \pdv{u}{t} + \sum_{m=2}^{+\infty} \frac{\Delta t^m}{m!} \left(\pmdv{u}{t}{m}\right) = \sum_j b_j \Delta x \pdv{u}{x} + \sum_{m = 2}^{+\infty} \frac{(j \Delta x)^m}{m!} \left(\pmdv{u}{x}{m}\right)
$$

For the linear convection equation:

$$
    \pdv{u}{t} + a \pdv{u}{x} = 0 \implies \sum_j j \cdot b_j = - \sig 
$$

This yields:

$$
\begin{aligned}
    &\pdv{u}{t} = -a \pdv{u}{x} + \frac{1}{\Delta t}\sum_{m = 2}^{+\infty} \frac{(j \Delta x)^m}{m!} \left(\pmdv{u}{x}{m}\right) - \sum_{m=2}^{+\infty} \frac{\Delta t^{m-1}}{m!} \left(\pmdv{u}{t}{m}\right) \\
\implies &\pdv{}{t} = -a \pdv{}{x} + \frac{1}{\Delta t}\sum_{m = 2}^{+\infty} \frac{(j \Delta x)^m}{m!} \left(\pmdv{}{x}{m}\right) - \sum_{m=2}^{+\infty} \frac{\Delta t^{m-1}}{m!} \left(\pmdv{}{t}{m}\right)
\end{aligned}
$$

Let us recall a simple formula where $y$ and $z$ are first order terms compared to $x$:

$$
    (x + y + z)^m = x\left(1 + \frac{y}{x} + \frac{z}{x} \right)^m = x^m + mx^{m-1} y + mx^{m-1} z + \mrm{HOT}
$$

Developing the derivative operator using the formula above:

$$
    \pdv{u}{t} + a \pdv{u}{x} = \sum_{m = 2}^{+\infty} \alpha_m \Delta x^{m-1} \left(\pmdv{u}{x}{m}\right) - \sum_{m=2}^{+\infty} \frac{(-\sig)^{m-1}}{(m-1)!} \sum_{n=2}^{+\infty} \alpha_n \Delta x^{m+n-2} \left(\pmdv{u}{x}{n+m-1}\right)
$$

where 

$$
    \alpha_p = \frac{1}{p!} \left[\sum_j b_j j^p - (-\sig)^p\right] \frac{\Delta x}{\Delta t}
$$

Then requiring the scheme to be of order $p$ yields the condition:

$$
    \pdv{u}{t} + a \pdv{u}{x} = \mc{O}(\Delta x^p) \implies \sum_j b_j j^p = (-\sig)^m \:\: \mrm{for}\: m \in \intint{0}{p}
$$

When the scheme is of order $p$ then the numerical scheme writes:

$$
    \pdv{u}{t} + a \pdv{u}{x} = \sum_{m=p}^{+\infty} a_{m+1} \Delta x^m \left(\pmdv{u}{x}{m+1}\right)
$$

and the $a_j$ are linked with the $\alpha_j$:

$$
\begin{aligned}
    a_{p+1} &= \alpha_{p+1} \\
    a_{p+2} &= \alpha_{p+2} + \sig \alpha_{p+1} \\
    a_{p+3} &= \alpha_{p+3} + \sig \alpha_{p+2} - \frac{\sig^2}{2}\alpha_{p+1} \\
    a_{p+4} &= \alpha_{p+4} + \sig \alpha_{p+3} - \frac{\sig^2}{2}\alpha_{p+2} + \frac{\sig^3}{6}\alpha_{p+1}\\
\end{aligned}
$$

### One parameter family schemes on the support $(i - 2, i - 1, i, i + 1)$

We consider second order schemes on the support $(i - 2, i - 1, i, i + 1)$. Thus the $b_j$ verify the relations:

$$
\begin{aligned}
b_{-2} + b_{-1} + b_0 + b_1 &= 1 \\
-2b_{-2} - b_{-1} + b_1 &= -\sig \\
4b_{-2} + b_{-1} + b_1 &= \sig^2
\end{aligned}
$$

Solving the system:

$$
\begin{aligned}
b_{-2} &= - \gamma\\
b_{-1} &= \frac{\sig(\sig + 1)}{2} + 3 \gamma\\
b_0 &= 1 - \sig^2 - 3 \gamma\\
b_1 &= \frac{\sig(\sig - 1)}{2} + \gamma
\end{aligned}
$$

Hence the scheme can be written as:

$$
    S = S_\mrm{LW} + \gamma H(-1, 3, -3, 1)
$$

Multiple schemes can then be derived:
$$
\begin{aligned}
    &S_\mrm{LW}(0, \frac{\sig(\sig + 1)}{2}, 1 - \sig^2, \frac{\sig(\sig - 1)}{2}) \\
    &S_\mrm{SOU}(\frac{\sig(\sig - 1)}{2}, \sig(2 - \sig), \frac{(1 - \sig)(2 - \sig)}{2}, 0) \\
    &S_\mrm{Fromm}(\frac{\sig(\sig - 1)}{4}, \frac{\sig(5 - \sig)}{4}, \frac{(1- \sig)(4 + \sig)}{4}, \frac{\sig(\sig - 1)}{2}) \\
    &S_\mrm{3}(\frac{\sig(\sig^2 - 1)}{6}, \frac{\sig(2 - \sig)(\sig + 1)}{2}, \frac{\sig(2 - \sig)(1 - \sig^2)}{2}, \frac{\sig(\sig - 1)(2 - \sig^2)}{6})
\end{aligned}
$$