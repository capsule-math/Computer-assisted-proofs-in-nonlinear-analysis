### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Periodic orbits"

using Markdown
using InteractiveUtils

# ╔═╡ a40ea48d-24a9-452f-a2b4-3e53d3de3763
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ b4c62765-aeb3-45cc-83b4-02bcc51d6c98
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
    using PlutoTeachingTools, PlutoUI, RadiiPolynomial, Plots
end

# ╔═╡ 044aabf3-b817-4266-93a9-a4cabce1ec4e
html"""
<style>
main {
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ 2ae8d07b-a05e-4136-a121-6317641916fe
TableOfContents(title = "Table of Contents"; indent = true, depth = 4, aside = true)

# ╔═╡ c33dc650-3f94-11ef-398a-8bbc4a2b69b8
md"""
We continue to do CAPs to find solutions of ODEs, but rather than looking at an initial value problem, we now look for periodic solutions.
The setting is again infinite dimensional, with some small differences:
- We use Fourier series rather than Taylor series.
- The Fourier coefficients are two-sided sequences of complex numbers.
- The convolution product is given by an infinite series.
- The problem of finding a closed orbit is global in the sense that there is no equivalent of "short time" integration.
"""

# ╔═╡ fa5da4c1-cbcc-413c-8807-b9d2b40d5cc0
md"# Motivating example: a forced, damped, nonlinear pendulum"

# ╔═╡ 3674e20c-4c06-4a5d-8d7c-24e408b6d979
md"""
Consider the second-order ODE

```math
u''(t) + \beta_1 u'(t) + \beta_2 u(t) - u(t)^2 = \beta_3 \cos(t),
```

where $\beta_1, \beta_2, \beta_3$ are parameters.
We aim to find a periodic solution $u(t)$.
Since the (minimal) period of the forcing term is $2\pi$, we look for a periodic solution of period $2 \pi$ as a Fourier series $u(t) = \sum_{k \in \mathbb{Z}} x_k e^{i k t}$.

We adopt the same strategy as before, organized into the following steps:
1. Identify an appropriate Banach space.
2. Reformulate the problem as a zero-finding problem $F(x) = 0$.
3. Compute an approximate zero $\bar{x}$.
4. Construct an approximate inverse $A$ of $DF(\bar{x})$.
5. Derive computable formulas for the bounds $Y, Z_1, Z_2$, and verify the radii polynomial inequalities.
"""

# ╔═╡ 1bfcc76a-266a-4784-b76b-3871445e603e
md"""
###### Why is this a "simple" example?

A more natural problem to look at is to find a periodic orbit in a system of ODEs.
In such an *autonomous* system of ODEs the period of a periodic solution is unknown a priori, and determining the period $\tau$ is thus part of the problem, in addition to finding the unknown $\tau$-periodic function $u(t)$.
Furthermore, when $u(t)$ is a periodic solution, so is $u(t-t_0)$ for any time shift $t_0$. In particular, periodic solutions are not isolated unless an additional phase condition is imposed.
A final complication is that for a *system* of equations the bookkeeping is more involved than for a single equation.
Although these difficulties can be overcome, the notation can get in the way of understanding.
Since we want to focus here on how to deal with periodic functions, we choose to look at a scalar, non-autonomous problem where the period is given.
"""

# ╔═╡ 64fffe67-78ae-4261-a149-61046fcf3bb4
md"## Sequence space"

# ╔═╡ 8f3e96dc-0cdc-411e-8f94-e9e47488abf3
md"""
The space we are interested in is the space of periodic analytic functions.
THe first step is to discretize this space as the sequence space ($\nu \ge 1$)

```math
\ell^1_{\nu, \mathbb{Z}}
\overset{\text{def}}{=}
\left\{
x \in \mathbb{C}^\mathbb{Z} \, : \, \| x \|_\nu \overset{\text{def}}{=} \sum_{k \in \mathbb{Z}} |x_k| \nu^{|k|} < \infty
\right\}.
```

We will choose the value of $\nu$ later.
"""

# ╔═╡ 766164f3-eef3-4798-bf7c-afdca52d3ba0
md"""
!!! lemma "Lemma"
	Let $\omega > 0$ and $\nu \ge 1$.
	If $x \in \ell^1_{\nu, \mathbb{Z}}$, then $\sum_{k \in \mathbb{Z}} x_k e^{i k \omega t}$ converges uniformly.
	If $\nu > 1$, then term-by-term differentiation up to any order of the Fourier series is justified.
"""

# ╔═╡ fe626b98-e13b-4af0-8c7b-b1442934d56b
md"""
!!! lemma "Lemma"
	 $\ell^1_{\nu, \mathbb{Z}}$ is a Banach space.
"""

# ╔═╡ f1cad266-547f-45fb-9391-cc0a45b65137
md"### Conjugate symmetry"

# ╔═╡ 9ab6ce17-0637-4008-be29-e0bf4fda4287
md"""
We denote the complex conjugate of $z \in \mathbb{C}$ by $z^*$.
For a real-valued periodic function $u$, its Fourier coefficients $x \in \ell^1_{\nu, \mathbb{Z}}$ satisfy $x_{-k} = x_k^*$ for all $k \in \mathbb{Z}$.
Conversely, if $x_{-k} = x_k^*$ for all $k \in \mathbb{Z}$, then the Fourier series $\sum_{k \in \mathbb{Z}} x_k e^{i k \omega t}$ is real-valued.

For $x \in \ell^1_{\nu, \mathbb{Z}}$, we denote by $x^\dagger \in \ell^1_{\nu, \mathbb{Z}}$ the conjugate given by $(x^\dagger)_k \overset{\text{def}}{=} x_{-k}^*$.
With this notation, $x = x^\dagger$ is equivalent to the Fourier series $\sum_{k \in \mathbb{Z}} x_k e^{i k \omega t}$ being real-valued.
"""

# ╔═╡ 069bbc98-d32c-4871-a317-7b44b1e45251
md"### Convolution product"

# ╔═╡ 47e39d6a-231e-40d1-906e-eeafe1b1f807
md"""
!!! lemma "Lemma (Convolution product)"
	Let $\omega > 0$.
	If $u(t) = \sum_{k \in \mathbb{Z}} x_k e^{i k \omega t}$ and $v(t) = \sum_{k \in \mathbb{Z}} y_k e^{i k \omega t}$, then $u(t) v(t) = \sum_{k \in \mathbb{Z}} (x * y)_k e^{i k \omega t}$ where

	```math
	(x * y)_k \overset{\text{def}}{=} \sum_{l \in \mathbb{Z}} x_{k-l} y_l, \qquad k \in \mathbb{Z}.
	```
"""

# ╔═╡ f91abfa7-8a8f-4204-95d9-637c622dd720
md"""
The product and the norm play together nicely: they give the space $\ell^1_{\nu, \mathbb{Z}}$ the structure of a Banach algebra, as expressed by the following lemma.
"""

# ╔═╡ 8f6431ed-346c-4c51-83c4-ea7ec12a8d60
md"""
!!! lemma "Lemma (Banach algebra property)"
	For any $\nu \ge 1$ and $x, y \in \ell^1_{\nu, \mathbb{Z}}$ we have $\| x * y \|_\nu \le \| x \|_\nu \| y \|_\nu$.
"""

# ╔═╡ 7796a035-7ba7-4ad2-b109-e1f1f966eda5
md"### Finite dimensional projection"

# ╔═╡ c32d9839-7cb7-4018-bab1-b2211f05b71e
md"""
We introduce a projection operator on a finite number of modes ($2K+1$ coefficients) $\Pi_K : \ell^1_{\nu, \mathbb{Z}} \to \ell^1_{\nu, \mathbb{Z}}$ given by

```math
(\Pi_K x)_k
\overset{\text{def}}{=}
\begin{cases}
x_k, & |k| \le K, \\
0, & |k| > K,
\end{cases}
```

together with its complement $\Pi_{> K} \overset{\text{def}}{=} I - \Pi_K$, where $I$ is the identity operator.
"""

# ╔═╡ 872359fa-dcb9-47da-9612-f46d6291fad5
md"## Zero-finding problem"

# ╔═╡ 8d7fb15d-32d3-4fbb-9d06-c8ec1ed7d00e
md"""
Substituting $u(t) = \sum_{k \in \mathbb{Z}} x_k e^{i k t}$ into the ODE leads to the infinite set of equations

```math
\lambda_k x_k + (x * x)_k + \beta_3 c_k = 0, \qquad k \in \mathbb{Z},
```

where

```math
\lambda_k \overset{\text{def}}{=} k^2 - i \beta_1 k - \beta_2
```

and $c$ represents the Fourier transform of $\cos(t)$, that is

```math
c_k
\overset{\text{def}}{=}
\begin{cases}
\frac{1}{2}, & k = \pm 1, \\
0, & k \ne \pm 1.
\end{cases}
```

It is convenient to introduce the diagonal linear operator

```math
(\Lambda x)_k
\overset{\text{def}}{=}
\lambda_k x_k.
```

When $\beta_1, \beta_2 \ne 0$, then the inverse $\Lambda^{-1}$, given by $(\Lambda^{-1} h)_k = \lambda_k^{-1} h_k$ for all $h \in \ell^1_{\nu, \mathbb{Z}}$, is well-defined as a bounded operator on $\ell^1_{\nu, \mathbb{Z}}$

Therefore, the problem of finding a periodic solution to the ODE corresponds to finding a zero of the mapping $F : \ell^1_{\nu, \mathbb{Z}} \to \ell^1_{\nu, \mathbb{Z}}$ given by

```math
F(x) \overset{\text{def}}{=} x + \Lambda^{-1} (x * x + \beta_3 c).
```
"""

# ╔═╡ ff2ca8bd-00b9-4830-a128-afdb24634896
md"""
!!! note "Remark"
	Dividing by $\lambda_k$ is a choice that simplifies the linear term, which is beneficial for some of the algebra.
	In terms of estimates it merely shifts the difficulty to the nonlinear term.
	Working with the alternative zero-finding problem $G(x) \overset{\text{def}}{=} \Lambda x + x*x + \beta_3 c$ is also a perfectly valid choice.
"""

# ╔═╡ 76961562-23b1-4aae-b761-6f95bdbbbab3
md"""
An implementation, using the RadiiPolynomial library, of the truncated zero-finding problem $\Pi_K \circ F \circ \Pi_K$ is given below.
"""

# ╔═╡ df46f5c1-dbc6-4cf1-8e97-36b22c6b4c51
function Λ⁻¹(β, space)
	∂ = project(Derivative(1), space, space)
	∂² = ∂ * ∂
	return -inv(∂² + β[1] * ∂ + β[2]*I)
end

# ╔═╡ f0e42405-8068-455b-b17c-ede3e259f01d
F(x, β, c, space) = x + Λ⁻¹(β, space) * (x*x + β[3] * c)

# ╔═╡ f5cc03ca-09eb-4419-9d51-78c956ee73d9
md"""
The Fréchet derivative of $F$ reads

```math
DF(x) = I + \Lambda^{-1} \mathcal{M}_{2x},
```

where

```math
\mathcal{M}_{2x} h = 2x * h, \qquad \forall h \in \ell^1_{\nu, \mathbb{Z}}.
```

In other words, $\mathcal{M}_{2x}$ is the multiplication operator associated with $2x$.

As a bounded linear operator acting on infinite sequences in $\ell^1_{\nu, \mathbb{Z}}$, $\mathcal{M}_{2x}$ can be visualized as a matrix with an infinite number of rows and columns:

```math
\mathcal{M}_{2x} =
\begin{pmatrix}
\ddots & \ddots & \ddots \\
\ddots & 2x_0    & 2x_{-1} & 2x_{-2} \\
\ddots & 2x_1    & 2x_0    & 2x_{-1} & \ddots \\
       & 2x_2    & 2x_1    & 2x_0    & \ddots \\
       &        & \ddots & \ddots & \ddots
\end{pmatrix}.
```

An implementation of $\Pi_K \circ DF(\Pi_K x) \circ \Pi_K$ of the truncated problem is given below.
"""

# ╔═╡ 58ad1864-9446-489c-b175-bb6bed64b385
DF(x, β, domain, codomain) = I + Λ⁻¹(β, codomain) * project(Multiplication(2x), domain, codomain)

# ╔═╡ 4412ca19-ade5-41f1-a178-f65b030e144f
md"## Numerical zero"

# ╔═╡ 4a9daa68-86ed-4327-98fc-03d436885d2e
md"""
We apply Newton's method to $\Pi_K \circ F \circ \Pi_K$ to find $\bar{x} \in \Pi_K \ell^1_{\nu, \mathbb{Z}}$ such that $\Pi_K F (\bar{x}) \approx 0$.
If the truncation dimension $K$ is large enough we may hope that also $\Pi_{> K} F (\bar{x}) \approx 0$, since we expect the coefficients of $\bar{x}$ to decrease for moderately large $k$.
"""

# ╔═╡ 095cd71f-1df2-4375-b400-9df58b585ea2
β = [0.1, 4.0, 1.0]

# ╔═╡ 98f5889c-2233-4b7a-b195-2380c1450d13
c = Sequence(Fourier(1, 1.0), [0.5, 0.0, 0.5])

# ╔═╡ e1f6ef9b-5818-415a-8945-97b65288a433
K = 10

# ╔═╡ 3fdf17d9-f6ac-46ba-a42d-8739a688794d
initial_guess = zeros(ComplexF64, Fourier(K, 1.0));

# ╔═╡ 60a73d6b-b56e-405d-8cda-ad4ac3c57257
F_DF(x) = (F(x, β, c, space(x)), DF(x, β, space(x), space(x)))

# ╔═╡ f05c0ebc-a3c6-49c6-a131-0eafd00ad098
x_bar, success = newton(F_DF, initial_guess)

# ╔═╡ d5800e45-8809-46d6-9f41-7565e5b1cbfd
begin
	plot(LinRange(0, 2π, 201), t -> real(x_bar(t)); linewidth = 3, label = "x_bar")
	xlims!(0, 2π)
	xlabel!("t")
end

# ╔═╡ 5f3d91ee-14c2-4b5e-8e2a-a607d086eb51
md"""
###### Initializing Newton's method

Here we were a bit lucky that a simple guess sufficed for Newton to converge.
For other parameter values, however, we may need to study the solutions of the ODE in more detail using numerical integration.
"""

# ╔═╡ a74b307d-20a3-424c-9e65-23b3bae465aa
md"## Constructing $A$"

# ╔═╡ 24e1a5af-569e-4ae8-a3d7-19b7962db48b
md"""
We use the projection to split $A$, the approximate inverse of $DF(\bar{x})$, into a finite part and a tail part:

```math
A
\overset{\text{def}}{=}
A_K \Pi_K + \Pi_{>K},
```

where $A_K$ is a linear map (constructed below) on $\Pi_K \ell^1_{\nu, \mathbb{Z}}$.
By construction, we can write $A$ in the block diagonal form

```math
A =
\begin{pmatrix}
A_K & 0 \\
0 & \Pi_{>K}
\end{pmatrix}.
```

Next, since the operator $\Pi_K DF(\bar{x}) \Pi_K$ is a finite-dimensional square matrix, we define $A_K$ to be its numerical inverse.
"""

# ╔═╡ e37811fa-cd67-424e-8707-25ad8534d293
A_K = inv(DF(x_bar, β, space(x_bar), space(x_bar)))

# ╔═╡ 0b4e2571-cbec-4a4d-aac6-dd8b32e3f474
md"""
!!! note "Remark"
	It is not necessary for the truncation dimensions of $\bar{x} \in \Pi_K \ell^1_{\nu, \mathbb{Z}}$ and $A_K$ to be the same.
	For clarity, we do not introduce two different truncation dimension parameters.
    Essentially, the truncation dimension for $\bar{x}$ controls the size of the residue, i.e. $Y$, while the truncation dimension for the approximate inverse $A_K$ controls the contractivity of the fixed-point operator, i.e. $Z_1$.
	In principle, these dimensions can be managed rather independently.
"""

# ╔═╡ 0d75d077-53c3-4f03-99d3-af796d973d66
md"""
## Verifying the contraction
"""

# ╔═╡ 35df9eed-5f22-4e26-b119-ea7354d2c762
md"""
!!! lemma "Lemma"
	Let $r > 0$, $\bar{x}$, $\tilde{x}$ such that $\tilde{x}$ is the unique zero of $F$ such that $\|\tilde{x} - \bar{x}\| \le r$.
	If $\bar{x}^\dagger = \bar{x}$, then $\tilde{x}^\dagger = \tilde{x}$.
"""

# ╔═╡ 34e0b4cb-147e-4959-a19b-6060d9f17bdc
conjugacy_symmetry!(x_bar)

# ╔═╡ 6df91aa3-ad36-4b9e-ba4a-4dfb46eb1d27
conj(x_bar[K:-1:-K]) == x_bar[-K:K]

# ╔═╡ 83f772ea-aeec-49a3-baa2-d2fd8771c0a7
md"""
We now verify that $T(x) \overset{\text{def}}{=} x - A F(x)$ is a contraction around $\bar{x} \in \Pi_K \ell^1_{\nu, \mathbb{Z}}$.
To rigorously compute the bounds $Y, Z_1, Z_2$, we enclose all relevant data into intervals.
"""

# ╔═╡ d1dc009f-4063-4c50-808e-2c2e5794dbe4
iβ = [interval(1)/interval(10), interval(4), interval(1)]

# ╔═╡ ee580d32-f357-4cad-87b6-0798b8423013
ic = interval(c)

# ╔═╡ aecf22cc-e97a-44ee-b2ea-b9b3f3a0ff42
ix_bar = interval(x_bar)

# ╔═╡ 6da19359-b4a6-4e6e-823a-0fc3af25f8d3
iA_K = interval(A_K)

# ╔═╡ b6cfde9d-357b-43e9-82ec-f88e8cb7ba8e
md"""
We choose a value for the weight $\nu$ for the Banach space $\ell^1_{\nu, \mathbb{Z}}$.
"""

# ╔═╡ 16678eaf-9b67-4fe0-ae0b-82f7bbc5e229
ν = interval(1.1) # does not have to be exactly 11/10

# ╔═╡ 58eceacb-d59b-4f4d-afdd-f616a241f2ca
X = Ell1(GeometricWeight(ν))

# ╔═╡ fcd0b92b-9b16-451e-a54a-6ee0158fa67b
md"""
We are now in the Newton-Kantorovich context and thus need to derive the appropriate bounds.
"""

# ╔═╡ 7c04c657-a86b-492f-b317-ace25f5f2f17
md"##### Computing $Y$"

# ╔═╡ eb2ae1fd-5062-436b-ac50-eb3b169bca42
md"""
This concerns the bound on the residue $A F(\bar{x})$.
Determining the residue requires only a finite computation.
Indeed, since $(\bar{x} * \bar{x})_k = 0$ for $|k| > 2K$ we infer that $\Pi_{>2K} F(\bar{x}) = 0$, and $F(\bar{x})$ can be computed explicitly.
In turn, this implies that $A F(\bar{x}) = \Pi_{2K} A \Pi_{2K} F(\bar{x})$.
We conclude that $\|A F(\bar{x})\|_X = \|\Pi_{2K} A \Pi_{2K} F(\bar{x})\|_X $ requires only a finite computation, which we can perform with interval arithmetic to get a rigorous upper bound $Y \ge \|A F(\bar{x})\|_X$.
"""

# ╔═╡ 6307aecd-b4f7-49ea-b060-4b9236a5ae10
function bound_Y(x_bar, β, c, A, X)
	K = order(x_bar)
	space = Fourier(2K, frequency(x_bar))

	Fx = F(x_bar, β, c, space)

	Id = interval(I)
	ext_A = complex(project(Id, space, space))
	ext_A[-K:K,-K:K] = coefficients(A)

	return norm(ext_A * Fx, X)
end

# ╔═╡ 9bf0e55e-de5f-4ec9-92e9-5169d3d14302
Y = bound_Y(ix_bar, iβ, ic, iA_K, X)

# ╔═╡ 0c454883-ff97-4874-99b3-cc4bc86c211f
md"##### Computing $Z_1$"

# ╔═╡ 907f8fb5-859b-4751-b483-2006ba358d5c
md"""
This concerns the bound on the operator $I - A DF(\bar{x})$.
In the analysis we will split the bound in a part which is computable by "brute force" and an estimate of the tail.
The triangle inequality is often helpful for such arguments, but in our case the characterization of the operator norm gives an additional tool.
In particular we will use that for any $K' \in \mathbb{N}$ we have that

```math
\| M \|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})} = \max \left( \| M \Pi_{K'} \|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})}, \| M \Pi_{>K'} \|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})} \right).
```

It follows from the expression for $DF(x)$ that $DF(\bar{x})$ is an finite bandwidth operator.

It follows from the block-diagonal structure of $A$ that this implies that

```math
A DF(\bar{x}) \Pi_{> 2K} =  (A_K \Pi_K + \Pi_{>K}) DF(\bar{x}) \Pi_{>2K} = \Pi_{>K} DF(\bar{x}) \Pi_{>2K},
```

and we see that this expression does not involve $A_K$.

Furthermore, it also follows from the block-diagonal structure of $A$ and the finite bandwidth of $DF(\bar{x})$ that

```math
A DF(\bar{x}) \Pi_{2K} = A \Pi_{3K} DF(\bar{x}) \Pi_{2K} = \Pi_{3K} A DF(\bar{x}) \Pi_{2K},
```

which can be represented as a matrix and requires only a finite computation to evaluate.

By choosing $K' = 2K$ we can thus split the computation of the norm in two parts:

```math
\begin{aligned}
\|I-A &DF(\bar{x})\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})} \\
&= \max \left( \|\Pi_{3K} [I-A DF(\bar{x})] \Pi_{2K}\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})}, \| [I-ADF(\bar{x})] \Pi_{>2K} \|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})} \right).
\end{aligned}
```

Provided that $K \ge \sqrt{\beta_2}$, we conclude that

```math
\|I-ADF(\bar{x})\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})} = \|\Pi_{3K+1}[I-ADF(\bar{x})]\Pi_{2K+1}\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})}.
```

The latter operator can be represented by a finite matrix and we can compute the operator norm $\|I-ADF(\bar{x})\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})}$ using interval arithmetic to obtain the bound $Z_1$.
"""

# ╔═╡ ac7b5711-b2ef-44dc-a272-c5d3f4e9cb74
function bound_Z₁(x_bar, β, A, X)
	K = order(x_bar)
	domain = Fourier(2K+1, frequency(x_bar))
	codomain = Fourier(3K+1, frequency(x_bar))

	DFx = DF(x_bar, β, domain, codomain)

	Id = interval(I)
	ext_A = complex(project(Id, codomain, codomain))
	ext_A[-K:K,-K:K] = coefficients(A)

	return opnorm(Id - ext_A * DFx, X)
end

# ╔═╡ 47ef9015-eae8-4b29-807f-d223754701d3
Z₁ = bound_Z₁(ix_bar, iβ, iA_K, X)

# ╔═╡ 1d43be44-58d5-43a5-a850-791077493097
md"##### Computing $Z_2$"

# ╔═╡ 6f346c67-879a-4427-afd4-6acdeb17af21
md"""
This concerns the bound on the operator $A[DF(x) - DF(\bar{x})]$.
Since this is a Lipschitz bound on the first derivative, one approach is to bound the second derivative.
However, here we just do a direct computation to see that $A[DF(x) - DF(\bar{x})]h = 2 A \Lambda^{-1}[(x - \bar{x}) * h]$, hence by the Banach algebra property and the definition of the operator norm we arrive at

```math
\|A [DF(x) - DF(\bar{x})] h\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})} \le 2 \|A \Lambda^{-1}\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})}\|x-\bar{x}\|_{\ell^1_{\nu, \mathbb{Z}}} \|h\|_{\ell^1_{\nu, \mathbb{Z}}}.
```

We thus need to bound $2 \|A \Lambda^{-1}\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})}$.

Since the operator $\Lambda^{-1}$ is diagonal, it has a natural restriction to $\Pi_K \ell^1_{\nu, \mathbb{Z}}$.
Assume $N \ge \sqrt{\beta_2}$.
Then $\|A \Lambda^{-1}\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})} = \max\left( \|A_K \Lambda^{-1} \|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})}, |\lambda_{K+1}|^{-1} \right)$ (Hint: consider $\|A \Lambda^{-1} e_k\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{Z}})}$ for $|k|\le K$ and for $|k| > K$ separately).
"""

# ╔═╡ 8c31a57c-0f2a-49a1-8567-c49409a9ab85
R = Inf

# ╔═╡ 271decf5-68cb-4c8f-9452-6354477f6131
function bound_Z₂(β, A, X)
	K = order(domain(A))
	return 2max(opnorm(Λ⁻¹(β, domain(A)) * A, X), abs((K+1)^2 - im * β[1] * (K+1) - β[2]))
end

# ╔═╡ 0d7146b3-3015-416b-aae8-697ac0850603
Z₂ = bound_Z₂(iβ, iA_K, X)

# ╔═╡ 27019122-b944-4adc-b453-02af8eb378ac
md"##### Finishing the CAP"

# ╔═╡ 7329dca5-8aca-4112-b950-d9ec6182a0ed
md"""
To determine the values of $r$ that satisfy the radii polynomial inequalities, we use the `interval_of_existence` function from RadiiPolynomial.
This function returns an interval of existence, within which all contained values satisfy the contraction conditions.
"""

# ╔═╡ 60b43603-9937-4869-b5fc-f74c1c18e710
ie = interval_of_existence(Y, Z₁, Z₂, R)

# ╔═╡ 238c9111-7dda-40e8-8297-9776fe77f7d9
r = inf(ie)

# ╔═╡ 77b91c27-ddb0-4115-b0d0-a78f2729dadc
md"""
To guarantee that the fixed-point of $T$ is a zero of $F$, we must check that $A$ is injective.
This property comes as a by-product of the contraction argument.
Indeed, $\|I - A DF(\bar{x})\| \le Z_1 < 1$ ensures that $A DF(\bar{x})$ is invertible, which implies that $A$ is surjective.
By construction, $A_K$ is surjective and $A$ is injective if and only if $A_K$ is injective.
Since $A_K$ can be represented by a finite square matrix, it follows that $A$ must also be injective.

Therefore, we have proved the existence of a zero $\tilde{x}$ of $F$ such that

```math
\tilde{x} = \bar{x} + \gamma,
```

where $\gamma \in \ell^1_{\nu, \mathbb{Z}}$ is not known explicitly, but satisfies $\| \gamma \|_\nu \le r$.
Our choice of norm also gives us a $C^0$-error bound for the function

```math
\begin{aligned}
\max_{t \in [0, 2\pi]} | \sum_{k \in \mathbb{Z}} \tilde{x}_k e^{i k t} - \sum_{|k| \le K} \bar{x}_k e^{i k t} |
&= \max_{t \in [0, 2\pi]} | \sum_{k \in \mathbb{Z}} (\tilde{x}_k - \bar{x}_k) e^{i k t} | \\
&\le \| \tilde{x} - \bar{x} \|_\nu \\
&\le r.
\end{aligned}
```
"""

# ╔═╡ 26d06cfe-24cc-41ec-87d2-8e489865b5b8
begin
	rigorous_eval(x, r, t) = interval(real(x(t)), r; format = :midpoint)

	ts = LinRange(0, 2π, 101)

	plot([[interval(ts[i], ts[i+1]), rigorous_eval(ix_bar, r, interval(ts[i], ts[i+1]))] for i = 1:length(ts)-1];
		label = "", color = :palegreen3)

	plot!(ts, t -> real(x_bar(t));
		label = "x_bar", color = :royalblue1, linewidth = 3)

	xlims!(0, 2π)
	xlabel!("t")
end

# ╔═╡ Cell order:
# ╟─a40ea48d-24a9-452f-a2b4-3e53d3de3763
# ╠═b4c62765-aeb3-45cc-83b4-02bcc51d6c98
# ╟─044aabf3-b817-4266-93a9-a4cabce1ec4e
# ╟─2ae8d07b-a05e-4136-a121-6317641916fe
# ╟─c33dc650-3f94-11ef-398a-8bbc4a2b69b8
# ╟─fa5da4c1-cbcc-413c-8807-b9d2b40d5cc0
# ╟─3674e20c-4c06-4a5d-8d7c-24e408b6d979
# ╟─1bfcc76a-266a-4784-b76b-3871445e603e
# ╟─64fffe67-78ae-4261-a149-61046fcf3bb4
# ╟─8f3e96dc-0cdc-411e-8f94-e9e47488abf3
# ╟─766164f3-eef3-4798-bf7c-afdca52d3ba0
# ╟─fe626b98-e13b-4af0-8c7b-b1442934d56b
# ╟─f1cad266-547f-45fb-9391-cc0a45b65137
# ╟─9ab6ce17-0637-4008-be29-e0bf4fda4287
# ╟─069bbc98-d32c-4871-a317-7b44b1e45251
# ╟─47e39d6a-231e-40d1-906e-eeafe1b1f807
# ╟─f91abfa7-8a8f-4204-95d9-637c622dd720
# ╟─8f6431ed-346c-4c51-83c4-ea7ec12a8d60
# ╟─7796a035-7ba7-4ad2-b109-e1f1f966eda5
# ╟─c32d9839-7cb7-4018-bab1-b2211f05b71e
# ╟─872359fa-dcb9-47da-9612-f46d6291fad5
# ╟─8d7fb15d-32d3-4fbb-9d06-c8ec1ed7d00e
# ╟─ff2ca8bd-00b9-4830-a128-afdb24634896
# ╟─76961562-23b1-4aae-b761-6f95bdbbbab3
# ╠═df46f5c1-dbc6-4cf1-8e97-36b22c6b4c51
# ╠═f0e42405-8068-455b-b17c-ede3e259f01d
# ╟─f5cc03ca-09eb-4419-9d51-78c956ee73d9
# ╠═58ad1864-9446-489c-b175-bb6bed64b385
# ╟─4412ca19-ade5-41f1-a178-f65b030e144f
# ╟─4a9daa68-86ed-4327-98fc-03d436885d2e
# ╠═095cd71f-1df2-4375-b400-9df58b585ea2
# ╠═98f5889c-2233-4b7a-b195-2380c1450d13
# ╠═e1f6ef9b-5818-415a-8945-97b65288a433
# ╠═3fdf17d9-f6ac-46ba-a42d-8739a688794d
# ╠═60a73d6b-b56e-405d-8cda-ad4ac3c57257
# ╠═f05c0ebc-a3c6-49c6-a131-0eafd00ad098
# ╟─d5800e45-8809-46d6-9f41-7565e5b1cbfd
# ╟─5f3d91ee-14c2-4b5e-8e2a-a607d086eb51
# ╟─a74b307d-20a3-424c-9e65-23b3bae465aa
# ╟─24e1a5af-569e-4ae8-a3d7-19b7962db48b
# ╠═e37811fa-cd67-424e-8707-25ad8534d293
# ╟─0b4e2571-cbec-4a4d-aac6-dd8b32e3f474
# ╟─0d75d077-53c3-4f03-99d3-af796d973d66
# ╟─35df9eed-5f22-4e26-b119-ea7354d2c762
# ╠═34e0b4cb-147e-4959-a19b-6060d9f17bdc
# ╠═6df91aa3-ad36-4b9e-ba4a-4dfb46eb1d27
# ╟─83f772ea-aeec-49a3-baa2-d2fd8771c0a7
# ╠═d1dc009f-4063-4c50-808e-2c2e5794dbe4
# ╠═ee580d32-f357-4cad-87b6-0798b8423013
# ╠═aecf22cc-e97a-44ee-b2ea-b9b3f3a0ff42
# ╠═6da19359-b4a6-4e6e-823a-0fc3af25f8d3
# ╟─b6cfde9d-357b-43e9-82ec-f88e8cb7ba8e
# ╠═16678eaf-9b67-4fe0-ae0b-82f7bbc5e229
# ╠═58eceacb-d59b-4f4d-afdd-f616a241f2ca
# ╟─fcd0b92b-9b16-451e-a54a-6ee0158fa67b
# ╟─7c04c657-a86b-492f-b317-ace25f5f2f17
# ╟─eb2ae1fd-5062-436b-ac50-eb3b169bca42
# ╠═6307aecd-b4f7-49ea-b060-4b9236a5ae10
# ╠═9bf0e55e-de5f-4ec9-92e9-5169d3d14302
# ╟─0c454883-ff97-4874-99b3-cc4bc86c211f
# ╟─907f8fb5-859b-4751-b483-2006ba358d5c
# ╠═ac7b5711-b2ef-44dc-a272-c5d3f4e9cb74
# ╠═47ef9015-eae8-4b29-807f-d223754701d3
# ╟─1d43be44-58d5-43a5-a850-791077493097
# ╟─6f346c67-879a-4427-afd4-6acdeb17af21
# ╠═8c31a57c-0f2a-49a1-8567-c49409a9ab85
# ╠═271decf5-68cb-4c8f-9452-6354477f6131
# ╠═0d7146b3-3015-416b-aae8-697ac0850603
# ╟─27019122-b944-4adc-b453-02af8eb378ac
# ╟─7329dca5-8aca-4112-b950-d9ec6182a0ed
# ╠═60b43603-9937-4869-b5fc-f74c1c18e710
# ╠═238c9111-7dda-40e8-8297-9776fe77f7d9
# ╟─77b91c27-ddb0-4115-b0d0-a78f2729dadc
# ╟─26d06cfe-24cc-41ec-87d2-8e489865b5b8
