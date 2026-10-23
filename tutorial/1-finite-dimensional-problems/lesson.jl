### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Finite-dimensional problems"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 92b0e241-aec3-455a-98e6-8c56fadf5e7c
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ 77f07d6e-f87e-4edc-bb77-3773b48e8ca6
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
    using PlutoTeachingTools, PlutoUI, Plots, IntervalArithmetic, RadiiPolynomial
end

# ╔═╡ a5169b4b-f0d4-4f66-8198-0aec5f8e135a
html"""
<style>
main {
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ b3845641-1537-4a27-8550-1eff30900a6b
TableOfContents(title = "Table of Contents"; indent = true, depth = 4, aside = true)

# ╔═╡ ed6f44f4-3696-494a-b8e3-b30184f5bb06
md"# Motivating example: a chaotic dynamical system"

# ╔═╡ c33dc650-3f94-11ef-398a-8bbc4a2b69b8
md"""
Consider the sequence defined by

$\begin{align}
x_{n+1} = \mu x_n (1-x_n),
\end{align}$

where $\mu$ is a parameter in $[0,4]$ and the initial condition $x_1$ is in $[0,1]$.
This is known as the [logistic map](https://en.wikipedia.org/wiki/Logistic_map).

The slider below can be used to vary the value of $\mu$.
"""

# ╔═╡ 9e0bfed7-7ccc-42b2-9b1a-15420d85c356
@bind μ Slider(0:0.1:4; default = 3.9)

# ╔═╡ e3ca081f-0e6a-4dcf-9436-ca4f34467379
μ

# ╔═╡ 79138e2e-6123-4923-8329-e5c17499d785
begin
	x1 = 0.4
	maxiter = 50
	xs = zeros(maxiter+1)
	xs[1] = x1
	for n = 1:maxiter
		xs[n+1] = μ*xs[n]*(1-xs[n])
	end
	plot(1:maxiter+1, xs; marker = (:circle, 5), legend = false)
	xlabel!("n")
	ylims!(-0.1, 1.1)
end

# ╔═╡ 1d9286fa-647a-4743-806f-9cea5aab6be1
md"""
As $\mu$ increases, the dynamics become progressively more complex.
This behaviour arises from a cascade of [period doubling bifurcations](https://en.wikipedia.org/wiki/Period-doubling_bifurcation), ultimately leading to [chaos](https://en.wikipedia.org/wiki/Chaos_theory).
A well-known result for discrete dynamical systems states that the existence of a period-3 orbit implies chaos.
"""

# ╔═╡ 1e56532e-9886-4a41-9ba0-1a619a85e014
md"""
!!! theorem "Theorem"
	Let $f : [0,1] \to [0,1]$ be continuous, and consider the dynamical system defined by $x_{n+1} = f(x_n)$.
	If there exists an orbit of period-3, i.e. $x_1 \in [0,1]$ such that $x_1 \ne x_2 \ne x_3$ and $x_4 = x_1$, then the system is chaotic.
	In particular, there exist orbits of any period.[^Sha64][^LY75]
"""

# ╔═╡ 8d6c631f-832a-48a9-891f-650f53a0ac6c
md"""
Let us fix a value of $\mu$ for which the dynamics seem complicated, say $\mu = 3.9$.
Our goal is to demonstrate the existence of a period-3 orbit.
We will begin by numerically approximating such an orbit and then rigorously establish the existence of a true period-3 orbit nearby.

While there are various methods to prove the existence of a period-3 orbit (both with and without computer assistance), the approach we present will extend naturally to more challenging infinite-dimensional problems.
"""

# ╔═╡ 429e5261-d7d0-47a3-94b5-8e7905a5bd53
md"# Contraction mapping and the Newton-Kantorovich theorem"

# ╔═╡ a27b45ca-77d1-408e-98e9-6b5248aca85e
md"""
Given a problem and an approximate solution $\bar{x}$, our goal is to define a fixed-point operator $T$ such that
- the solutions to our problem are in one-to-one correspondence with fixed-points of $T$,
-  $T$ is a contraction on a small neighborhood of $\bar{x}$.

The specific form of $T$ may vary depending on the nature of the problem, but we can outline a general strategy for constructing such an operator.
First, we provide sufficient conditions that can be verified in practice to demonstrate that an operator $T$ is contracting on a small neighborhood of $\bar{x}$.
Although this first module focuses on finite-dimentional problems, we state these conditions in the more general setting of a [Banach space](https://en.wikipedia.org/wiki/Banach_space), as they will play an important role in the next modules.
"""

# ╔═╡ a801e2c7-c907-40d2-9976-bba19c590b2b
md"""
!!! theorem "Theorem (Contraction mapping)"
	Let $\mathcal{X}$ be a Banach space, $\bar{x} \in \mathcal{X}$, and $T : \mathcal{X} \to \mathcal{X}$ a continuously differentiable map.
	Assume that there exist a constant $Y$ and a non-decreasing map $Z : [0, \infty) \to [0, \infty)$ such that

	```math
	\begin{align}
	\| T(\bar{x}) - \bar{x} \| &\le Y, \\
	\| DT(x) \| &\le Z(\| x-\bar{x} \|), \qquad \forall x \in \mathcal{X}.
	\end{align}
	```

	If there exists $r > 0$ such that

	```math
	\begin{align}
	Y + \int_0^r Z(s) \, ds &\le r, \\
	Z(r) &< 1,
	\end{align}
	```

	then $T$ has a unique fixed point $\tilde{x}$ such that $\| \tilde{x} - \bar{x} \| \le r$.
"""

# ╔═╡ 6caaf267-3015-4a66-975f-ec1c3a6f86f8
md"""
!!! note "Remark"
	- The bound $Z$ is actually only needed locally, i.e. we only need $Z(s)$ for $s \le r$.
	  Therefore, one can fix a priori some $R > 0$, and weaken the assumption by only asking for $Z$ to be defined on $[0, R]$ and to satisfy $\| DT(x) \| \le Z(\| x-\bar{x} \|)$ for all $x \in \mathcal{X}$ such that $\| x-\bar{x} \| \le R$.
	  Of course, we are then only allowed to consider $r \in (0, R]$.
	- In practice, studying carefully how $\| DT(x) \|$ depends on $\| x-\bar{x} \|$ helps us to get better bounds.
	  However, we can sometimes get away with coarser estimates.
	  Specifically, we derive a constant $Z_*$ such that $\| DT(x) \| \le Z_*$ for all $x \in \mathcal{X}$ such that $\| x-\bar{x} \| \le R$, and then use $Z(s) = Z_*$ for all $s \in [0, R]$.
	  The conditions on $r$ then simplify, and, if $Z_* < 1$, we can take any $r \ge \frac{Y}{1-Z_*}$.
"""

# ╔═╡ 77380230-67f2-4794-bb51-39ceeac7bef1
md"""
For a dynamical system $x_{n+1} = f(x_n)$, a natural fixed-point operator for period-3 orbits is given by $T(x) \overset{\text{def}}{=} f^3(x)$, which represents $f$ composed with itself three times.
However, in general such a $T$ has no reason to be contracting.

To develop a more general procedure for constructing a fixed-point problem, one can start by finding a map $F$ whose isolated zeros are solutions to our problem.
For this example, we consider the mapping $F : \mathbb{R}^3 \to \mathbb{R}^3$ given by

```math
F(x_1, x_2, x_3) \overset{\text{def}}{=}
\begin{pmatrix}
f(x_1) - x_2 \\
f(x_2) - x_3 \\
f(x_3) - x_1
\end{pmatrix}.
```

Given an approximate solution $\bar{x}$, one can consider the Newton-like operator $T(x) \overset{\text{def}}{=} x - DF(\bar{x})^{-1} F(x)$.
This operator should be contracting in a neighborhood of $\bar{x}$, as $DT(\bar{x}) = 0$.

For finite-dimensional problems of moderate size, computing the inverse of $DF(\bar{x})$ is possible.
However, in infinite-dimensional problems, working with $DF(\bar{x})^{-1}$ can be challenging.
To address this, we often consider an approximate inverse $A \approx DF(\bar{x})^{-1}$ and define the fixed-point operator $T(x) \overset{\text{def}}{=} x - AF(x)$.
In the following results, we will operate under this framework, where $A$ is not specified and could be $DF(\bar{x})^{-1}$ or another suitable approximation.
"""

# ╔═╡ 08dab1b4-c78d-4e85-89b6-cf9b485030e0
md"""
!!! corollary "Corollary (Newton-Kantorovich)"
	Let $\mathcal{X}$ and $\mathcal{Y}$ be two Banach spaces, $\bar{x} \in \mathcal{X}$, $F : \mathcal{X} \to \mathcal{Y}$ a continuously differentiable map, $A : \mathcal{Y} \to \mathcal{X}$ am injective linear map, and $R > 0$.
	Assume that there exist constants $Y$, $Z_1$ and $Z_2$ such that

	```math
	\begin{align}
	\| AF(\bar{x}) \| &\le Y, \\
	\| I-ADF(\bar{x}) \| &\le Z_1, \\
	\| A(DF(x)-DF(\bar{x})) \| &\le Z_2 \| x-\bar{x} \|, \qquad \forall x \in \mathcal{X} \text{ such that } \| x-\bar{x} \| \le R.
	\end{align}
	```

	If there exists $r \in (0,R]$ such that

	```math
	\begin{align}
	Y + Z_1 r + \frac{1}{2} Z_2 r^2 &\le r, \\
	Z_1 + Z_2 r &< 1,
	\end{align}
	```

	then $F$ has a unique zero $\tilde{x}$ such that $\| \tilde{x} - \bar{x} \| \le r$.
"""

# ╔═╡ 3dc209e6-d7c0-464b-a477-0c46ede9483e
md"""
This corollary is a simplified version of the [Newton-Kantorovich Theorem](https://en.wikipedia.org/wiki/Kantorovich_theorem) [^Ort68].
Variations of this result can be found throughout the literature and are commonly used in computer-assisted proofs (CAPs) (see e.g. [^BL15] [^NPW19] and the references therein).
"""

# ╔═╡ f0d6ab5b-13cf-4408-aaae-93f959ce29c4
md"# Interval arithmetic"

# ╔═╡ 24ddc2c5-7976-4c67-8525-f6bb5d2e04a5
md"""
When we want to apply the contraction mapping theorem with $\bar{x}$ as an approximate solution obtained using the computer, we also need to compute values for the bounds $Y$ and $Z$.
However, these calculations typically rely on floating-point arithmetic, which introduces rounding errors.
For instance, if we execute `Y = norm(A * F(x̄), 1)` on the computer, we cannot be certain that $Y$ is an upper bound for $\|A F(\bar{x})\|_1$.

[Interval arithmetic](https://en.wikipedia.org/wiki/Interval_arithmetic) provides a robust method to control rounding errors and ensuring guaranteed results from computer calculations.
In this section, we offer a brief description of interval arithmetic; for more detailed discussions, please refer to [^Moo79] and [^Tuc11].
"""

# ╔═╡ 3b123092-1241-4a2a-bb61-36c7f68a0b35
md"## Overview"

# ╔═╡ 66ca51ad-0f2e-4f3c-8dd0-e66649b224d9
md"""
The central concept of interval arithmetic is to represent any real number by an interval, whose endpoints are floating-point numbers.
More precisely, for $a \in \mathbb{R}$, we define the interval $[\underline{a}, \overline{a}]$ where $\underline{a}, \overline{a}$ are floating-point numbers satisfying $\underline{a} \le a \le \overline{a}$.
While we give up the hope of representing numbers exactly, we recover guaranteed information: the real number $a$ is contained within the interval stored in the computer.

The rules of interval arithmetic ensure that this containment property is preserved for arithmetic operations.
For instance, consider two intervals $[a]$ and $[b]$ of the form $[\underline{a}, \overline{a}]$ and $[\underline{b}, \overline{b}]$, which enclose the reals numbers $a$ and $b$ respectively.
When we perform the operation $[c] = [a] + [b]$, the resulting interval $[c] = [\underline{c}, \overline{c}]$ is computed as follows:
- for $\underline{c}$, we take $\underline{a} + \underline{b}$ **rounded downward**, and
- for $\overline{c}$, we take $\overline{a}+\overline{b}$ **rounded upward**.
In particular, the real number $c = a + b$ is contained within the interval $[c]$.

Similarly, if $[c]$ is an interval $[\underline{c}, \overline{c}]$ and we compute $[d] = e^{[c]}$, the resulting interval $[d]$ contains $e^c$.

In this course, we use the Julia library [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl), though many other interval arithmetic libraries are available in different programming languages.
Note that IntervalArithmetic is automatically included when using [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl).
This library provides standard arithmetic operations as well as implementations of elementary functions (such as `exp`, `log`, `cos`, etc.) that comply with the rules of interval arithmetic.
"""

# ╔═╡ f8eab26f-c893-4d2e-b4e8-6b59f33cbc9c
a = interval(2, 4)

# ╔═╡ 556163e0-80d2-4f4e-b198-b8e0922438db
cos(exp(sqrt(a)))

# ╔═╡ 97fadca1-fcf9-46b8-aa07-c615ca0deb7f
md"""
!!! note "Remark"
	Even when doing computer-assisted proofs, it’s not necessary to use interval arithmetic for every computation.
	For finding an approximate solution $\bar{x}$, floating-point calculations are perfectly fine.
	It’s only when computing the bounds for the contraction theorem, and verifying the radii polynomial inequalities, that interval arithmetic must be used.
"""

# ╔═╡ 0bf1cce8-6caf-40d5-93a6-b495b764d943
md"## Beware of typed-in floating-point numbers"

# ╔═╡ c09b05d4-0325-486a-a433-6b32fbb7dfc7
md"""
When you execute a command like `x = 0.1` in Julia, the variable `x` is assigned the closest floating-point number approximating $1/10$.
Since $1/10$ cannot be represented exactly in base 2, the value stored in `x` is not precisely equal to the typed-in number `0.1`.
In fact, it turns out to be slightly larger than $1/10$, but it could just as easily have been smaller (as seen with `x = 0.3`).
For more information on this topic, you can refer to this [webpage](https://0.30000000000000004.com).
"""

# ╔═╡ 2d68d26d-6e10-407f-8227-515ddadd9599
x = 0.1

# ╔═╡ 78f17262-7ffc-4d74-835c-f17863817b9a
x > 1//10

# ╔═╡ c44fd902-7655-43b5-82cc-53c2ecc4f77b
ix = interval(1)/interval(10)

# ╔═╡ 4c15bb79-2591-4f2a-9243-ff811de70df7
in_interval(1//10, ix)

# ╔═╡ 68243410-d2df-4606-82a5-53d515403e40
md"# Back to period-3 implies chaos"

# ╔═╡ 9a0e21bd-da27-48e6-9a08-b20b3cb805bd
md"""
Let’s recap what we have covered so far:
- Proving the existence of a period-3 orbit is sufficient to establish chaotic dynamics.
- We reformulated the task of finding a period-3 orbit as solving an equation $F = 0$.
- We derived sufficient conditions, in terms of bound estimates $Y, Z_1, Z_2$, to determine the radius of a ball, centered at a numerical approximation $\bar{x}$, that contains a (locally unique) true zero $\tilde{x}$ of $F$.
- We introduced interval arithmetic to rigorously estimate the bounds $Y$, $Z_1$, and $Z_2$ with computer assistance.

Now, we walk through the implementation of a computer-assisted proof to demonstrate a period-3 orbit for $\mu = 3.9$.

We start by defining the logistic map in Julia:
"""

# ╔═╡ b3121e23-03c4-46f2-89b7-c6b39f8efdad
f(x, μ) = μ * x * (1 - x) # logistic map

# ╔═╡ aeda376d-b1d1-4613-9964-05c24daf4b9a
md"""
and the zero-finding problem:
"""

# ╔═╡ c31f9913-0f16-4e74-b395-97038706be3c
function F(x, μ)
    v = zeros(eltype(x), 3)
    v[1] = f(x[1], μ) - x[2]
	v[2] = f(x[2], μ) - x[3]
	v[3] = f(x[3], μ) - x[1]
	return v
end

# ╔═╡ 129c8b50-1ba3-427e-a306-5be57c9755a3
md"""
To apply the Newton-Kantorovich theorem we also need the first derivate:
"""

# ╔═╡ d2f64f7e-2acb-4bf6-a961-3e2d8a87c409
Df(x, μ) = μ * (1 - 2x) # derivative of the logistic map w.r.t. x

# ╔═╡ 8591da5a-0b9a-4bff-8adc-fe1d3d9fca09
function DF(x, μ) # derivative of the zero-finding problem w.r.t. x
	M = zeros(eltype(x), 3, 3)
	M[1,1] = Df(x[1], μ)
	M[1,2] = -1
	M[2,2] = Df(x[2], μ)
	M[2,3] = -1
	M[3,1] = -1
	M[3,3] = Df(x[3], μ)
	return M
end

# ╔═╡ fa0be782-c09e-45b7-9cbb-d0b3161c3317
md"## Numerical zero"

# ╔═╡ cf53d067-ec00-4c06-b4ae-65ed7aa19d34
md"""
We apply [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method) to $F$ using the `newton` function provided in the RadiiPolynomial library.
"""

# ╔═╡ ef6afc18-f779-4370-8ccc-0a9494341166
initial_guess = [0.5, 1.0, 0.1]

# ╔═╡ a71e7ddf-e4d7-4200-abfb-b248883200e1
F_DF(x) = (F(x, μ), DF(x, μ)) # returns `F` and `DF` at `μ = 3.9`

# ╔═╡ b1ff2f19-31b7-4064-9ac0-e600ccd937d9
x_bar, success = newton(F_DF, initial_guess)

# ╔═╡ 8cc5d9d9-6630-43eb-85b7-2e658da6ec03
md"""
The `newton` function returns two outputs: the result of the iterations and a value indicating whether the method converged -- `true` if it did (as it does here) and `false` otherwise.
"""

# ╔═╡ 81d97321-ff79-42d6-adc6-3d66f91c438b
md"""
###### Initializing Newton's method

Finding a good starting point for Newton's method is not always easy.
Here we were a bit lucky that a simple guess worked well.
For other parameter values, however, we may need to take a closer look at the iterations of the map $f$.
"""

# ╔═╡ a2c9c052-1628-4ccd-8ac6-e2b16e0fc687
md"## Constructing $A$"

# ╔═╡ 32b5d274-3d0f-41aa-b1cf-85dfe40d0f76
md"""
As mentioned earlier, $A$ doesn’t need to be exactly $DF(\bar{x})^{-1}$; it only needs to be a reasonably close approximation.
For finite-dimensional problems, this inverse can be efficiently computed using built-in numerical algorithms.
"""

# ╔═╡ 16ec745d-0353-471f-b177-c25c5440e860
A = inv(DF(x_bar, μ))

# ╔═╡ 8ba3efdb-bff7-4211-bcd0-1537ae0420ed
md"## Verifying the contraction"

# ╔═╡ d9117644-f54f-42ce-94df-5933e57d765c
md"""
We now verify that $T(x) \overset{\text{def}}{=} x - A F(x)$ is a contraction around $\bar{x} \in \mathbb{R}^3$.
To rigorously compute the bounds $Y, Z_1, Z_2$, we enclose all relevant data into intervals.
"""

# ╔═╡ 787a3624-7a81-4f9f-9bca-f28e914a8f07
# since 39/10 is not representable as a floating-point number
iμ = interval(39)/interval(10)

# ╔═╡ a182a2c2-2e20-44f5-92ba-e85805050170
ix_bar = interval(x_bar)

# ╔═╡ f6759116-c48c-4a54-b073-8c2a0adc97cf
iA = interval(A)

# ╔═╡ 761a7427-6561-4754-8815-b1bf491a17da
md"##### Computing $Y$"

# ╔═╡ 7a037025-57ef-4b62-a880-8a2538096a9e
Y = norm(iA * F(ix_bar, iμ), 1)

# ╔═╡ 7a451708-9427-425a-93e0-0e66cad2c80b
md"##### Computing $Z_1$"

# ╔═╡ b7d36207-5f2c-41de-b142-7e8ee8700f79
Z₁ = opnorm(I - iA * DF(ix_bar, iμ), 1)

# ╔═╡ 7c1876f9-ccb7-4084-841f-316b1c1f6308
md"##### Computing $Z_2$"

# ╔═╡ 1a196987-95f2-48e0-9648-66373dce0ed5
md"""
```math
\|A(DF(x) - DF(\bar{x}))\|_1 \le \|A\|_1 \|DF(x) - DF(\bar{x})\|_1 = 2 \mu \|A\|_1 \| x - \bar{x} \|_1, \qquad \forall x \in \mathbb{R}^3.
```

So $Z_2 \ge 2 \mu \|A\|_1$ is a constant.
In other words, we have a Lipschitz control for $DF$ on the whole space $\mathbb{R}^3$, and we can freely choose $R = \infty$.
"""

# ╔═╡ 0b1f112c-d8b0-4903-966c-8ac6f4be69dc
R = Inf

# ╔═╡ 6f38f2e6-170c-4fea-b260-e55ca3fc3849
Z₂ = interval(2) * iμ * opnorm(iA, 1)

# ╔═╡ f7db0f90-473b-4265-8918-1b7acb08ebd2
md"##### Finishing the CAP"

# ╔═╡ afadd46e-458b-4975-925c-2add34ec75af
md"""
To determine the values of $r$ that satisfy the radii polynomial inequalities, we use the `interval_of_existence` function from RadiiPolynomial.
This function returns an interval of existence, within which all contained values satisfy the contraction conditions.
"""

# ╔═╡ 980ecce4-f02d-4976-b953-587b31009adb
ie = interval_of_existence(Y, Z₁, Z₂, R)

# ╔═╡ 958265d4-f604-47a4-8195-64ba21275f6d
md"""
The largest value in the interval represents the maximum radius of the ball centered at $\bar{x}$ within which the true zero $\tilde{x}$ of $F$ is proven to be unique.
On the other hand, the smallest value provides the sharpest error bound.
"""

# ╔═╡ 9eb669e1-f191-46cf-a3cb-b441b7d7f9e5
r = inf(ie)

# ╔═╡ c1a3ff98-5846-4601-87f7-b94b4f1fdf17
md"""
In fact, we cannot yet conclude that $\tilde{x}$ is a zero of $F$ since we did not check that $A$ was injective.
This property actually comes as a by-product of the contraction argument.
Indeed, $\|I - A DF(\bar{x})\| \le Z_1 < 1$ ensures that $A DF(\bar{x})$ is invertible, which implies that $A$ is surjective.
Since $A$ is a square matrix, it follows that $A$ must also be injective.

Finally, we verify that the true zero $\tilde{x}$ of $F$ is a period-3 orbit, that is $\tilde{x}_1 \ne \tilde{x}_2 \ne \tilde{x}_3$.
Since we have found $r$ such that $\| \tilde{x} - \bar{x} \|_1 \le r$, it follows that $\tilde{x}_j \in [\bar{x}_j - r, \bar{x}_j + r]$ for $j = 1, 2, 3$.
We then confirm that the three intervals $[\bar{x}_j - r, \bar{x}_j + r]$, for $j = 1, 2, 3$, are mutually disjoint:
"""

# ╔═╡ 611dc902-dd8d-4ec8-a2c0-1e9ac949a128
enclosure_x_tilde_1 = interval(ix_bar[1], r; format = :midpoint)

# ╔═╡ 61e71daf-8a6d-4c27-bccd-91f5e4db992c
enclosure_x_tilde_2 = interval(ix_bar[2], r; format = :midpoint)

# ╔═╡ 5f7e039e-b1a7-40a2-a076-ee56215da40c
enclosure_x_tilde_3 = interval(ix_bar[3], r; format = :midpoint)

# ╔═╡ 0f3d6fcf-2794-4d09-ae35-9dc5db3d0fa8
isdisjoint_interval(enclosure_x_tilde_1, enclosure_x_tilde_2, enclosure_x_tilde_3)

# ╔═╡ 1a9d9d5b-22bd-4fea-b480-817d80c63f8e
md"""
###### A final note on interval arithmetic

In the code provided in this notebook we used the Julia interval arithmetic library [IntervalArithmetic.jl](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/). We took a bit of liberty by taking advantage of automatic conversion of floats and integers to intervals, when variables of these types are combined in elementary operations. The interval arithmetic library is not entirely happy about this:
"""

# ╔═╡ 26a665d6-77c1-4bb4-b2df-04e05a80ee80
isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)

# ╔═╡ 0a2208b6-f90c-4587-978f-33643fd4bdba
md"# References"

# ╔═╡ b8d56ba6-01da-4604-8dad-3e63ec203fd4
md"""
[^BL15]: J. B. van den Berg and J.-P. Lessard. Rigorous numerics in dynamics. *Notices Of The American Mathematical Society*, 62(9), 2015.
[^LY75]: T.-Y. Li and J. A. Yorke. Period three implies chaos. *The American Mathematical Monthly*, 82(10):985--992, 1975.
[^Moo79]: R. E. Moore. *Methods and applications of interval analysis*. SIAM, 1979
[^Ort68]: J. M. Ortega. The Newton-Kantorovich theorem. *The American Mathematical Monthly*, 75(6):658--660, 1968.
[^NPW19]: M. T. Nakao, M. Plum, and Y. Watanabe. *Numerical Verification Methods and Computer-Assisted Proofs for Partial Differential Equations*. Springer Singapore, 2019.
[^Sha64]: A. N. Sharkowskii. Co-existence of the cycles of a continuous mapping of the line into itself. *Ukrainian Mathematical Journal*, 16(1), 1964.
[^Tuc11]: W. Tucker. *Validated numerics: a short introduction to rigorous computations*. Princeton University Press, 2011.
"""

# ╔═╡ Cell order:
# ╟─92b0e241-aec3-455a-98e6-8c56fadf5e7c
# ╠═77f07d6e-f87e-4edc-bb77-3773b48e8ca6
# ╟─a5169b4b-f0d4-4f66-8198-0aec5f8e135a
# ╟─b3845641-1537-4a27-8550-1eff30900a6b
# ╟─ed6f44f4-3696-494a-b8e3-b30184f5bb06
# ╟─c33dc650-3f94-11ef-398a-8bbc4a2b69b8
# ╟─9e0bfed7-7ccc-42b2-9b1a-15420d85c356
# ╠═e3ca081f-0e6a-4dcf-9436-ca4f34467379
# ╟─79138e2e-6123-4923-8329-e5c17499d785
# ╟─1d9286fa-647a-4743-806f-9cea5aab6be1
# ╟─1e56532e-9886-4a41-9ba0-1a619a85e014
# ╟─8d6c631f-832a-48a9-891f-650f53a0ac6c
# ╟─429e5261-d7d0-47a3-94b5-8e7905a5bd53
# ╟─a27b45ca-77d1-408e-98e9-6b5248aca85e
# ╟─a801e2c7-c907-40d2-9976-bba19c590b2b
# ╟─6caaf267-3015-4a66-975f-ec1c3a6f86f8
# ╟─77380230-67f2-4794-bb51-39ceeac7bef1
# ╟─08dab1b4-c78d-4e85-89b6-cf9b485030e0
# ╟─3dc209e6-d7c0-464b-a477-0c46ede9483e
# ╟─f0d6ab5b-13cf-4408-aaae-93f959ce29c4
# ╟─24ddc2c5-7976-4c67-8525-f6bb5d2e04a5
# ╟─3b123092-1241-4a2a-bb61-36c7f68a0b35
# ╟─66ca51ad-0f2e-4f3c-8dd0-e66649b224d9
# ╠═f8eab26f-c893-4d2e-b4e8-6b59f33cbc9c
# ╠═556163e0-80d2-4f4e-b198-b8e0922438db
# ╟─97fadca1-fcf9-46b8-aa07-c615ca0deb7f
# ╟─0bf1cce8-6caf-40d5-93a6-b495b764d943
# ╟─c09b05d4-0325-486a-a433-6b32fbb7dfc7
# ╠═2d68d26d-6e10-407f-8227-515ddadd9599
# ╠═78f17262-7ffc-4d74-835c-f17863817b9a
# ╠═c44fd902-7655-43b5-82cc-53c2ecc4f77b
# ╠═4c15bb79-2591-4f2a-9243-ff811de70df7
# ╟─68243410-d2df-4606-82a5-53d515403e40
# ╟─9a0e21bd-da27-48e6-9a08-b20b3cb805bd
# ╠═b3121e23-03c4-46f2-89b7-c6b39f8efdad
# ╟─aeda376d-b1d1-4613-9964-05c24daf4b9a
# ╠═c31f9913-0f16-4e74-b395-97038706be3c
# ╟─129c8b50-1ba3-427e-a306-5be57c9755a3
# ╠═d2f64f7e-2acb-4bf6-a961-3e2d8a87c409
# ╠═8591da5a-0b9a-4bff-8adc-fe1d3d9fca09
# ╟─fa0be782-c09e-45b7-9cbb-d0b3161c3317
# ╟─cf53d067-ec00-4c06-b4ae-65ed7aa19d34
# ╠═ef6afc18-f779-4370-8ccc-0a9494341166
# ╠═a71e7ddf-e4d7-4200-abfb-b248883200e1
# ╠═b1ff2f19-31b7-4064-9ac0-e600ccd937d9
# ╟─8cc5d9d9-6630-43eb-85b7-2e658da6ec03
# ╟─81d97321-ff79-42d6-adc6-3d66f91c438b
# ╟─a2c9c052-1628-4ccd-8ac6-e2b16e0fc687
# ╟─32b5d274-3d0f-41aa-b1cf-85dfe40d0f76
# ╠═16ec745d-0353-471f-b177-c25c5440e860
# ╟─8ba3efdb-bff7-4211-bcd0-1537ae0420ed
# ╟─d9117644-f54f-42ce-94df-5933e57d765c
# ╠═787a3624-7a81-4f9f-9bca-f28e914a8f07
# ╠═a182a2c2-2e20-44f5-92ba-e85805050170
# ╠═f6759116-c48c-4a54-b073-8c2a0adc97cf
# ╟─761a7427-6561-4754-8815-b1bf491a17da
# ╠═7a037025-57ef-4b62-a880-8a2538096a9e
# ╟─7a451708-9427-425a-93e0-0e66cad2c80b
# ╠═b7d36207-5f2c-41de-b142-7e8ee8700f79
# ╟─7c1876f9-ccb7-4084-841f-316b1c1f6308
# ╟─1a196987-95f2-48e0-9648-66373dce0ed5
# ╠═0b1f112c-d8b0-4903-966c-8ac6f4be69dc
# ╠═6f38f2e6-170c-4fea-b260-e55ca3fc3849
# ╟─f7db0f90-473b-4265-8918-1b7acb08ebd2
# ╟─afadd46e-458b-4975-925c-2add34ec75af
# ╠═980ecce4-f02d-4976-b953-587b31009adb
# ╟─958265d4-f604-47a4-8195-64ba21275f6d
# ╠═9eb669e1-f191-46cf-a3cb-b441b7d7f9e5
# ╟─c1a3ff98-5846-4601-87f7-b94b4f1fdf17
# ╠═611dc902-dd8d-4ec8-a2c0-1e9ac949a128
# ╠═61e71daf-8a6d-4c27-bccd-91f5e4db992c
# ╠═5f7e039e-b1a7-40a2-a076-ee56215da40c
# ╠═0f3d6fcf-2794-4d09-ae35-9dc5db3d0fa8
# ╟─1a9d9d5b-22bd-4fea-b480-817d80c63f8e
# ╠═26a665d6-77c1-4bb4-b2df-04e05a80ee80
# ╟─0a2208b6-f90c-4587-978f-33643fd4bdba
# ╟─b8d56ba6-01da-4604-8dad-3e63ec203fd4
