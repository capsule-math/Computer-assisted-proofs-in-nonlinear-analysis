### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Taylor integration"

using Markdown
using InteractiveUtils

# ╔═╡ db5e06bb-9eea-40d1-9328-da7e9cf77cb3
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ b21e396d-9088-404e-a992-3994e9b0eb7a
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
    using PlutoTeachingTools, PlutoUI, RadiiPolynomial, Plots
end

# ╔═╡ 551ab699-e4c5-4748-b567-9a29b8783b33
html"""
<style>
main {
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ 5ed309bf-dc26-4008-b199-f1de5078130e
TableOfContents(title = "Table of Contents"; indent = true, depth = 4, aside = true)

# ╔═╡ 3e086655-fc49-42bb-803e-27b82c7a7577
md"""
In this lecture we will learn how to solve initial value problems for ODEs using Taylor series.
"""

# ╔═╡ 86c84f47-a64b-43f9-a7e0-2eeb78b0d41c
md"# Motivating example: solving the logistic equation"

# ╔═╡ ab841d58-1385-4d5b-89aa-7df1e9667978
md"""
Consider the [initial value problem](https://en.wikipedia.org/wiki/Initial_value_problem) (IVP) for the [logistic equation](https://en.wikipedia.org/wiki/Logistic_function#Logistic_differential_equation):

```math
\begin{cases}
u'(t) = u(t) (1 - u(t)), & t \ge 0, \\
u(0) = ϕ.
\end{cases}
```

Given the initial condition $\phi = 1/2$ and $\nu > 0$, we will search for a function $u$ as a power series $u(t) = \sum_{n \ge 0} x_n t^n$, defined for all $t \in [-\nu, \nu]$.

Of course, this ODE is simple enough to be solvable by hand.
One can verify that $u(t) = 1/(1 + e^{-t})$ is the exact solution.
A numerical plot is displayed below.
"""

# ╔═╡ 5fe08415-9185-4054-8dec-a657cfb4484b
begin
	plot(LinRange(-2.5, 2.5, 51), t -> 1/(1+exp(-t)); label = "1/(1+exp(-t))", linewidth = 3)
	xlims!(-2.5, 2.5)
	xlabel!("t")
end

# ╔═╡ 7cd532e8-f370-47e2-a2c1-3a5004cb8e8c
md"""
We adopt the contraction strategy:
1. Identify an appropriate Banach space.
2. Reformulate the problem as a zero-finding problem $F(x) = 0$.
3. Compute an approximate zero $\bar{x}$.
4. Construct an approximate inverse $A$ of $DF(\bar{x})$.
5. Derive computable formulas for the bounds $Y, Z_1, Z_2$, and verify the radii polynomial inequalities.
"""

# ╔═╡ 2be11063-c9dd-4c36-af2b-2d80fad17fcb
md"## Sequence space"

# ╔═╡ 529d0a8a-c650-4b70-bd1c-bdd0621e06df
md"""
We are interested in the space of analytic functions with a radius of convergence strictly greater than some fixed $\nu > 0$.
If $\nu$ is chosen too large, the computer-assisted proof (CAP) will fail.

The first step is to discretize this space as a sequence space.
However, if the radius of convergence of the function lies in the interval $(0, 1)$, the coefficients of the series grow, which leads to numerical instability.
To mitigate this, it is good practice to re-scale the problem so that the radius of convergence is at least $1$, and we may as well scale it to $1$.
We rescale the time variable by introducing

```math
w(s) = u(\nu s), \qquad \text{for all } s \in [-1, 1].
```

Then, the sequence space we use is

```math
\ell^1_\mathbb{N} \overset{\text{def}}{=} \left\{ x \in \mathbb{R}^\mathbb{N} \, : \, \| x \|_1 \overset{\text{def}}{=} \sum_{n \ge 0} |x_n| < \infty \right\}.
```

The properties of this space were discussed in Module 2.
Here, we briefly recall them.
Suppose $u(t) = \sum_{n \ge 0} x_n t^n$ and $v(t) = \sum_{n \ge 0} y_n t^n$.
Their product is given by

```math
u(t) v(t) = \sum_{n \ge 0} (x * y)_n t^n,
```

where the convolution operation $x * y$ is given by

```math
(x * y)_n \overset{\text{def}}{=} \sum_{l=0}^n x_{n-l} y_l, \qquad n \ge 0.
```

Moreover, the convolution satisfies the inequality

```math
\| x * y \|_1 \leq \| x \|_1 \| y \|_1.
```
"""

# ╔═╡ 42e2418b-55bc-4e17-8cb7-445f8eee2cc6
md"""
As discussed in Module 2, infinite-dimensional problems are addressed by decomposing the sequence space into a finite-dimensional subspace and the complementary infinite-dimensional space.
To formalize this, we define the projection operator $\Pi_N : \ell^1_\mathbb{N} \to \ell^1_\mathbb{N}$ which extracts the first $N+1$ coefficients of a sequence.
Specifically, for all $x \in \ell^1_\mathbb{N}$, we have

```math
(\Pi_N x)_n
\overset{\text{def}}{=}
\begin{cases}
x_n, & n \le N, \\
0, & n > N.
\end{cases}
```

The complementary projection, which removes the first $N+1$ modes, is defined as $\Pi_{> N} \overset{\text{def}}{=} I - \Pi_N$, where $I$ is the identity operator.
"""

# ╔═╡ c454f8a4-5553-494a-8d08-7c7d756805a0
md"## Zero-finding problem"

# ╔═╡ 0691e2f7-110c-42f7-95f6-056699578782
md"""
Substituting $w(t) = \sum_{n \ge 0} x_n t^n$ into the ODE, with initial condition $w(0) = \phi$, leads to the recurrence relation

```math
x_0 = \phi, \qquad x_n = \frac{\nu}{n} (x_{n-1} - (x * x)_{n-1}), \quad n \ge 1.
```

We can re-write this set of equations more concisely as

```math
x = \phi + \nu \Lambda (x - x * x),
```

where

```math
\Lambda x \overset{\text{def}}{=}
\begin{cases}
0, & n = 0, \\
\displaystyle \frac{1}{n} x_{n-1}, & n \ge 1.
\end{cases}
```

We slightly abuse notation and identify $\phi$ with its sequence of Taylor coefficients $(\phi, 0, \dots) \in \ell^1_\mathbb{N}$.

Therefore, the problem of finding the solution $u$ corresponds to finding a zero of the mapping $F : \ell^1_\mathbb{N} \to \ell^1_\mathbb{N}$ given by

```math
F(x) \overset{\text{def}}{=}
x - \phi - \nu \Lambda (x - x * x).
```
"""

# ╔═╡ 4ee023ac-e513-41f8-bfde-9d2daec71c57
md"""
This zero-finding problem is of course identical to integrating first the ODE $u(t) = \phi + \int_0^t u(s)(1-u(s)) \, ds$, and plugging the power series ansatz into it.
In this sense, $\Lambda$ is nothing more than the representation of the integral operator in $\ell^1_\mathbb{N}$.

An implementation, using the RadiiPolynomial library, of the zero-finding problem $F$ is given below.
"""

# ╔═╡ 4e6473ec-02c6-4ed5-ab05-8f0086853313
f(x) = x * (1 - x)

# ╔═╡ d53d988f-da50-4839-aae6-d7bdac46d4ac
F(x, ϕ, ν, space) = x - ϕ - ν * project(integrate(f(x)), space)

# ╔═╡ 9516114e-de9b-4705-9197-52737d050602
md"""
The Fréchet derivative is given by

```math
DF(x) = I - \nu \Lambda \circ (I - \mathcal{M}_{2x}).
```
"""

# ╔═╡ abc4b3dd-518c-4239-87bb-e9688b34a8db
Df(x) = 1 - 2x

# ╔═╡ 7e5be3bf-ea25-4269-b05f-8263d5ede515
DF(x, ν, domain, codomain) = I -
    ν * project(Integral(1), domain, codomain) * project(Multiplication(Df(x)), domain, codomain)

# ╔═╡ e6b00d47-adf6-4935-b3c8-5b30ece1e2f7
md"""
###### An other contraction argument

We can actually take advantage of the explicit recurrence relation to solve the IVP.
Using interval arithmetic, it’s possible to compute the Taylor coefficients $\tilde{x}_0, \dots, \tilde{x}_N$ exactly, up to any desired order $N > 0$.
Once we have these, we can define the truncated sequence $\hat{x} = (\tilde{x}_0, \dots, \tilde{x}_N, 0, \dots)$, where the terms beyond $N$ are set to zero.
To estimate the truncation error, observe that the recurrence relation naturally leads to a fixed-point operator $T$.
By focusing on the complementary subspace $\Pi_{> N} \ell^1_\mathbb{N}$, we can define $G : \Pi_{> N} \ell^1_\mathbb{N} \to \Pi_{> N} \ell^1_\mathbb{N}$ by

```math
G(h) \overset{\text{def}}{=} T(\hat{x} + h).
```

Intuitively, $h$ represents the correction to the tail of $\hat{x}$.
The key idea is that, for a sufficiently large $N$, this operator $G$ becomes contracting around $0 \in \Pi_{> N} \ell^1_\mathbb{N}$.

There are two main differences with the CAPs done so far:
- we do not construct a Newton-like fixed-point operator, and
- we must solve **exactly** a finite number of equations.

That said, we do not pursue this approach in this lecture, as it is specific to power series.
Our broader goal is to develop techniques that can be applied to other kinds of series.
In particular, [Chebyshev polynomials](https://en.wikipedia.org/wiki/Chebyshev_polynomials) are especially relevant in the context of IVPs.
"""

# ╔═╡ 3a6bae8b-933e-4384-81a8-e3d3bb5459c7
md"## Numerical zero"

# ╔═╡ eef70ab3-fb71-4cfe-bca6-4ab0d578aca3
N = 20 # order of the Taylor series

# ╔═╡ 4a173aa3-1672-4a21-a9d0-e64a9336144a
ϕ = 0.5 # initial condition

# ╔═╡ 8ead7dda-965c-4540-a77d-d37ed4d8aea9
ν = 2.5 # time rescaling

# ╔═╡ 5266d1d6-45c5-497b-9d2a-22a6f4423efa
initial_guess = zeros(Taylor(N));

# ╔═╡ fc74e5bf-ccfa-460b-a69c-29c50d9400e7
F_DF(x) = (F(x, ϕ, ν, Taylor(N)), DF(x, ν, Taylor(N), Taylor(N)))

# ╔═╡ 6edbe03d-77b5-4c8f-884a-fe008a061228
x_bar, success = newton(F_DF, initial_guess)

# ╔═╡ 25e484d2-e228-45a8-9a8f-53ec54326fd9
begin
	plot(LinRange(-ν, ν, 101), t -> x_bar(t/ν); linewidth = 3, label = "x_bar")
	xlims!(-ν, ν)
	xlabel!("t")
end

# ╔═╡ 3395bf15-3e17-4039-bba1-839933a9f481
md"""
###### Initializing Newton's method

In this example, Newton converged with a simple guess.
In general, however, we need to study the solutions of the ODE in more detail using numerical integration.
"""

# ╔═╡ 641caade-ff17-4484-af95-988b92977bb5
md"## Constructing $A$"

# ╔═╡ d55fb36b-a0a8-41cd-a50f-9495ede7f2ad
md"""
We use the projection to split $A$, the approximate inverse of $DF(\bar{x})$, into a finite part and a tail part:

```math
A
\overset{\text{def}}{=}
A_N \Pi_N + \Pi_{>N},
```

where $A_N$ is a linear map (constructed below) on $\Pi_N \ell^1_\mathbb{N}$.
By construction, we can write $A$ in the block diagonal form

```math
A =
\begin{pmatrix}
A_N & 0 \\
0 & \Pi_{>N}
\end{pmatrix}.
```

Next, since the operator $\Pi_N DF(\bar{x}) \Pi_N$ is a finite-dimensional square matrix, we define $A_N$ to be its numerical inverse.
"""

# ╔═╡ 29ce245e-e838-42f1-a524-6f8059b0291b
A_N = inv(DF(x_bar, ν, Taylor(N), Taylor(N)))

# ╔═╡ 97e52e3e-fd1c-4dff-8102-fd8b509d45f6
md"## Verifying the contraction"

# ╔═╡ 43ec73f0-778b-4571-a04f-56832939194e
md"""
We now verify that $T(x) \overset{\text{def}}{=} x - A F(x)$ is a contraction around $\bar{x} \in \Pi_N \ell^1_\mathbb{N}$.
To rigorously compute the bounds $Y, Z_1, Z_2$, we enclose all relevant data into intervals.
"""

# ╔═╡ 9ec0cbd7-e74d-42b4-b911-9b164a0c7abc
iϕ = interval(0.5) # 0.5 is exactly representable as a Float64

# ╔═╡ ddcbee38-f3ec-47ff-81fd-7cae07482929
iν = interval(2.5) # 2.5 is exactly representable as a Float64

# ╔═╡ bc8c3312-e7d8-413e-9c4c-95bf46ba802e
ix_bar = interval(x_bar)

# ╔═╡ 5b83f453-55a7-4a38-8e74-aea6af53a837
iA_N = interval(A_N)

# ╔═╡ 2e0c0dc6-031a-46a8-9c37-27ffab80db2c
md"##### Computing $Y$"

# ╔═╡ 10898ff1-2a1f-4369-8f3f-0f096eca6b75
function bound_Y(x_bar, ϕ, ν, A)
	N = order(x_bar)
	space = Taylor(2N+1)

	Fx = F(x_bar, ϕ, ν, space)

	Id = interval(I)
	ext_A = project(Id, space, space)
	ext_A[0:N,0:N] = coefficients(A)

	return norm(ext_A * Fx, 1)
end

# ╔═╡ 34152734-45ff-4dfe-bc25-3d60ba4c9d19
Y = bound_Y(ix_bar, iϕ, iν, iA_N)

# ╔═╡ e3ae92f8-e06e-4931-b424-38c893001471
md"##### Computing $Z_1$"

# ╔═╡ 82b4d717-835e-4577-b4bd-410f5fecab41
function bound_Z₁(x_bar, ν, A)
	N = order(x_bar)
	space = Taylor(N+1)

	DFx = DF(x_bar, ν, space, space)

	Id = interval(I)
	ext_A = project(Id, space, space)
	ext_A[0:N,0:N] = coefficients(A)

	return opnorm(Id - ext_A * DFx, 1)
end

# ╔═╡ 4bce0a13-ca2f-46c4-9ef3-12f07d525fbf
Z₁ = bound_Z₁(ix_bar, iν, iA_N)

# ╔═╡ 765ee153-3a73-4f0e-81ff-a9852461e4ef
md"##### Computing $Z_2$"

# ╔═╡ 09c77d32-4541-458e-86a0-e4fb3c324e93
R = Inf

# ╔═╡ b2d7aa0e-f9aa-4dc4-aff7-4474e9a6c9b1
function bound_Z₂(ν, A)
	N = order(domain(A))
	Λ = project(Integral(1), domain(A), domain(A), Interval{Float64})
	return 2ν * max(opnorm(A * Λ, 1), inv(interval(N)))
end

# ╔═╡ 031ab1e8-d479-4197-919c-20416f19aa9f
Z₂ = bound_Z₂(iν, iA_N)

# ╔═╡ 7b38e714-75bf-4a56-87f1-f8db5db517e3
md"##### Finishing the CAP"

# ╔═╡ b53f9f0a-c4d8-4102-ac55-0bc4cea3cbdf
md"""
To determine the values of $r$ that satisfy the radii polynomial inequalities, we use the `interval_of_existence` function from RadiiPolynomial.
This function returns an interval of existence, within which all contained values satisfy the contraction conditions.
"""

# ╔═╡ 7fc34159-4b00-4e5e-b969-b4550b38ab53
ie = interval_of_existence(Y, Z₁, Z₂, R)

# ╔═╡ 67008bc9-20ae-435b-9228-438c022fa9f2
r = inf(ie)

# ╔═╡ 559d9d44-3501-4a7c-ad4c-f81d633f1a15
md"""
To guarantee that the fixed-point of $T$ is a zero of $F$, we must check that $A$ is injective.
This property comes as a by-product of the contraction argument.
Indeed, $\|I - A DF(\bar{x})\| \le Z_1 < 1$ ensures that $A DF(\bar{x})$ is invertible, which implies that $A$ is surjective.
By construction, $A_N$ is surjective and $A$ is injective if and only if $A_N$ is injective.
Since $A_N$ can be represented by a finite square matrix, it follows that $A$ must also be injective.

Therefore, we have proved the existence of a zero $\tilde{x}$ of $F$ such that

```math
\tilde{x} = \bar{x} + \gamma,
```

where $\gamma \in \ell^1_\mathbb{N}$ is not known explicitly, but satisfies $\| \gamma \|_1 \le r$.
Our choice of norm also gives us a $C^0$-error bound for the function

```math
\begin{aligned}
\sup_{t \in [-\nu,\nu]} | u(t) - \sum_{n=0}^N \bar{x}_n (t/\nu)^n |
&= \sup_{t \in [-1,1]} | w(t) - \sum_{n=0}^N \bar{x}_n t^n | \\
&= \sup_{t \in [-1,1]} | \sum_{n \ge 0} \tilde{x}_n t^n - \sum_{n=0}^N \bar{x}_n t^n | \\
&= \sup_{t \in [-1,1]} | \sum_{n \ge 0} (\tilde{x}_n - \bar{x}_n) t^n | \\
&\le \| \tilde{x} - \bar{x} \|_1 \\
&\le r.
\end{aligned}
```
"""

# ╔═╡ fb6bd3c5-6622-457f-b57f-bdd0d395d7ab
begin
	rigorous_eval(x, r, t) = interval(x(t), r; format = :midpoint)

	ts = LinRange(-ν, ν, 101)

	plot([[interval(ts[i], ts[i+1]), rigorous_eval(ix_bar, r, interval(ts[i]/ν, ts[i+1]/ν))] for i = 1:length(ts)-1];
		label = "", color = :palegreen3)

	plot!(ts, t -> 1/(1+exp(-t));
		label = "1/(1+exp(-t))", color = :salmon, linewidth = 3)

	plot!(ts, t -> x_bar(t/ν);
		label = "x_bar", color = :royalblue1, linewidth = 3, linestyle = :dash)

	xlims!(-ν, ν)
	xlabel!("t")
end

# ╔═╡ Cell order:
# ╟─db5e06bb-9eea-40d1-9328-da7e9cf77cb3
# ╠═b21e396d-9088-404e-a992-3994e9b0eb7a
# ╟─551ab699-e4c5-4748-b567-9a29b8783b33
# ╟─5ed309bf-dc26-4008-b199-f1de5078130e
# ╟─3e086655-fc49-42bb-803e-27b82c7a7577
# ╟─86c84f47-a64b-43f9-a7e0-2eeb78b0d41c
# ╟─ab841d58-1385-4d5b-89aa-7df1e9667978
# ╟─5fe08415-9185-4054-8dec-a657cfb4484b
# ╟─7cd532e8-f370-47e2-a2c1-3a5004cb8e8c
# ╟─2be11063-c9dd-4c36-af2b-2d80fad17fcb
# ╟─529d0a8a-c650-4b70-bd1c-bdd0621e06df
# ╟─42e2418b-55bc-4e17-8cb7-445f8eee2cc6
# ╟─c454f8a4-5553-494a-8d08-7c7d756805a0
# ╟─0691e2f7-110c-42f7-95f6-056699578782
# ╟─4ee023ac-e513-41f8-bfde-9d2daec71c57
# ╠═4e6473ec-02c6-4ed5-ab05-8f0086853313
# ╠═d53d988f-da50-4839-aae6-d7bdac46d4ac
# ╟─9516114e-de9b-4705-9197-52737d050602
# ╠═abc4b3dd-518c-4239-87bb-e9688b34a8db
# ╠═7e5be3bf-ea25-4269-b05f-8263d5ede515
# ╟─e6b00d47-adf6-4935-b3c8-5b30ece1e2f7
# ╟─3a6bae8b-933e-4384-81a8-e3d3bb5459c7
# ╠═eef70ab3-fb71-4cfe-bca6-4ab0d578aca3
# ╠═4a173aa3-1672-4a21-a9d0-e64a9336144a
# ╠═8ead7dda-965c-4540-a77d-d37ed4d8aea9
# ╠═5266d1d6-45c5-497b-9d2a-22a6f4423efa
# ╠═fc74e5bf-ccfa-460b-a69c-29c50d9400e7
# ╠═6edbe03d-77b5-4c8f-884a-fe008a061228
# ╟─25e484d2-e228-45a8-9a8f-53ec54326fd9
# ╟─3395bf15-3e17-4039-bba1-839933a9f481
# ╟─641caade-ff17-4484-af95-988b92977bb5
# ╟─d55fb36b-a0a8-41cd-a50f-9495ede7f2ad
# ╠═29ce245e-e838-42f1-a524-6f8059b0291b
# ╟─97e52e3e-fd1c-4dff-8102-fd8b509d45f6
# ╟─43ec73f0-778b-4571-a04f-56832939194e
# ╠═9ec0cbd7-e74d-42b4-b911-9b164a0c7abc
# ╠═ddcbee38-f3ec-47ff-81fd-7cae07482929
# ╠═bc8c3312-e7d8-413e-9c4c-95bf46ba802e
# ╠═5b83f453-55a7-4a38-8e74-aea6af53a837
# ╟─2e0c0dc6-031a-46a8-9c37-27ffab80db2c
# ╠═10898ff1-2a1f-4369-8f3f-0f096eca6b75
# ╠═34152734-45ff-4dfe-bc25-3d60ba4c9d19
# ╟─e3ae92f8-e06e-4931-b424-38c893001471
# ╠═82b4d717-835e-4577-b4bd-410f5fecab41
# ╠═4bce0a13-ca2f-46c4-9ef3-12f07d525fbf
# ╟─765ee153-3a73-4f0e-81ff-a9852461e4ef
# ╠═09c77d32-4541-458e-86a0-e4fb3c324e93
# ╠═b2d7aa0e-f9aa-4dc4-aff7-4474e9a6c9b1
# ╠═031ab1e8-d479-4197-919c-20416f19aa9f
# ╟─7b38e714-75bf-4a56-87f1-f8db5db517e3
# ╟─b53f9f0a-c4d8-4102-ac55-0bc4cea3cbdf
# ╠═7fc34159-4b00-4e5e-b969-b4550b38ab53
# ╠═67008bc9-20ae-435b-9228-438c022fa9f2
# ╟─559d9d44-3501-4a7c-ad4c-f81d633f1a15
# ╟─fb6bd3c5-6622-457f-b57f-bdd0d395d7ab
