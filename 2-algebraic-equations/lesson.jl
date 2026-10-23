### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Algebraic equations"

using Markdown
using InteractiveUtils

# ╔═╡ ec5584b5-8b4d-4149-97c7-54fe0ee7fb8b
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ 51bbeeac-7ab4-4afc-9c2d-3f32b65c9b6d
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
    using PlutoTeachingTools, PlutoUI, RadiiPolynomial, Plots
end

# ╔═╡ 2af9b21a-2dcf-497f-bd7c-9ebbf80b9428
html"""
<style>
main {
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ 11810156-e8c9-4145-af2a-93e0757a8cc9
TableOfContents(title = "Table of Contents"; indent = true, depth = 4, aside = true)

# ╔═╡ 44fc05a5-a0c1-4428-adee-edb5b8b869e3
md"# Motivating example: the inverse function"

# ╔═╡ 9162d799-c8a8-44cd-8921-b6fa4c6f291f
md"""
Given $v(t) = 2 + t$, consider the problem of finding the function $u : [-1, 1] \to \mathbb{R}$ such that

```math
u(t) v(t) = 1, \qquad t \in [-1, 1].
```

This constitutes the first example of infinite-dimensional problems.

We will search for a function $u$ as a power series $u(t) = \sum_{n \ge 0} x_n t^n$, defined for all $t \in [-1, 1]$.
A numerical plot of $(2+t)^{-1}$ is displayed below.
"""

# ╔═╡ 37827087-01ad-4fe7-bcb6-df2239ea0bca
begin
	plot(LinRange(-1, 1, 51), t -> 2 + t; label = "2+t", linewidth = 3)
	plot!(LinRange(-1, 1, 51), t -> inv(2 + t); label = "1/(2+t)", linewidth = 3)
	xlims!(-1, 1)
	xlabel!("t")
end

# ╔═╡ d794367f-d346-416f-b8f1-1f442f876e20
md"""
We follow the same strategy outlined in Module 1, with some key adaptations.
1. Identify an appropriate Banach space (an additional step compared to the previous module on finite-dimensional problems).
2. Reformulate the problem as a zero-finding problem $F(x) = 0$.
3. Compute an approximate zero $\bar{x}$.
4. Construct an approximate inverse $A$ of $DF(\bar{x})$.
5. Derive computable formulas for the bounds $Y, Z_1, Z_2$, and verify the radii polynomial inequalities.
"""

# ╔═╡ f6da61f7-64a3-4c33-b9ba-7917ffb67b2a
md"## Sequence space"

# ╔═╡ 65cebd36-0474-42ee-b983-70c5f9c1614f
md"""
The space we are interested in is the one of analytic functions with radius of convergence strictly greater than $1$.
The first step is to discretize this space as the sequence space

```math
\ell^1_\mathbb{N} \overset{\text{def}}{=} \left\{ x \in \mathbb{R}^\mathbb{N} \, : \, \| x \|_1 \overset{\text{def}}{=} \sum_{n \ge 0} |x_n| < \infty \right\}.
```
"""

# ╔═╡ 2dc02de6-f250-44db-b347-7d4d65c61d89
md"""
!!! lemma "Lemma"
	 $x \in \ell^1_\mathbb{N}$ if and only if $\sum_{n \ge 0} x_n t^n$ converges uniformly on $[-1, 1]$.
"""

# ╔═╡ f9d764b9-024a-4845-98ad-c3e79a4a5cd9
md"""
!!! lemma "Lemma"
	 $\ell^1_\mathbb{N}$ is a Banach space.
"""

# ╔═╡ 56bf8063-cc76-47ee-9074-f6a59466ec60
md"### Cauchy product"

# ╔═╡ 3cb18cba-18e5-4e6f-8b16-eccd1dcef9fd
md"""
!!! lemma "Lemma (Cauchy product)"
	If $u(t) = \sum_{n \ge 0} x_n t^n$ and $v(t) = \sum_{n \ge 0} y_n t^n$, then $u(t) v(t) = \sum_{n \ge 0} (x * y)_n t^n$ where

	```math
	(x * y)_n \overset{\text{def}}{=} \sum_{l=0}^n x_{n-l} y_l, \qquad n \ge 0.
	```
"""

# ╔═╡ 7390523a-ecf1-4b7e-a3e6-cce50d43083b
md"""
The product and the norm play together nicely: they give the space $\ell^1_\mathbb{N}$ the structure of a Banach algebra, as expressed by the following lemma.
"""

# ╔═╡ 2a4448c5-e14a-4791-926f-172f3cb6e1f0
md"""
!!! lemma "Lemma (Banach algebra property)"
	For any $x, y \in \ell^1_\mathbb{N}$ we have $\| x * y \|_1 \le \| x \|_1 \| y \|_1$.
"""

# ╔═╡ faec7172-c732-470c-88a0-8de5d92075b0
md"### Finite dimensional projection"

# ╔═╡ 4290872b-923b-45dc-81b0-3e280a5a77c1
md"""
We introduce a projection operator on a finite number of modes ($N+1$ coefficients) $\Pi_N : \ell^1_\mathbb{N} \to \ell^1_\mathbb{N}$ given by

```math
(\Pi_N x)_n
\overset{\text{def}}{=}
\begin{cases}
x_n, & n \le N, \\
0, & n > N,
\end{cases}
```

together with its complement $\Pi_{> N} \overset{\text{def}}{=} I - \Pi_N$, where $I$ is the identity operator.

The projections let us  break down the problem into two parts: a finite-dimensional part, handled by the computer, and an infinite-dimensional part (often called the *tail*), handled entirely by pen-and-paper analysis.
"""

# ╔═╡ 8356c027-ead5-4073-9d0e-b25f41bf9543
md"""
###### Why is this so complicated?

The strange looking subscript $> N$ is in fact meant to indicate clearly which index values are involved.
It avoids writing down infinite matrices or formulas with a lot of indices.
"""

# ╔═╡ c11ce6e3-4e23-46bc-8bbd-217908d966af
md"## Zero-finding problem"

# ╔═╡ b6655a7c-32d0-473a-ac1d-ddb96e31de81
md"""
Writing $v(t) = \sum_{n \ge 0} y_n t^n$ with $y \overset{\text{def}}{=} (2, 1, 0, \dots) \in \ell^1_\mathbb{N}$, and substituting $u(t) = \sum_{n \ge 0} x_n t^n$ into the algebraic equation leads to the infinite set of equations

```math
(x * y)_0 = 1, \qquad (x * y)_n = 0, \quad n \ge 1.
```

Therefore, the problem of finding the inverse of $v$ corresponds to finding a zero of the mapping $F : \ell^1_\mathbb{N} \to \ell^1_\mathbb{N}$ given by

```math
F(x) \overset{\text{def}}{=} x * y - 1.
```
"""

# ╔═╡ 266d8eb2-fc1f-42ca-802e-90a4a358b313
md"""
!!! note "Remark"
    Here we slightly abuse notation by identifying the constant $1$ with its sequence of Taylor coefficients $(1, 0, \dots) \in \ell^1_\mathbb{N}$.
    This is consistent with the convention that $1$ denotes the identity element of the algebra $(\ell^1_\mathbb{N}, *)$.
"""

# ╔═╡ 0dab4ba7-4ad0-4ea8-a4c9-5892427969e7
md"""
An implementation, using the RadiiPolynomial library, of the zero-finding problem $F$ is given below.
"""

# ╔═╡ 4f37b759-371e-4042-82d0-720523718021
F(x, y) = x * y - 1

# ╔═╡ 17013b31-0c79-4dd0-a6e3-5946528d69f6
md"## Numerical zero"

# ╔═╡ 151f9f5f-65ab-4ccb-94a6-c530d924962f
md"""
The map $F$ is affine and its Fréchet derivative $DF(x) \in \mathscr{B}(\ell^1_\mathbb{N})$, at some $x \in \ell^1_\mathbb{N}$, reads

```math
[DF(x)]h = y * h, \qquad \forall h \in \ell^1_\mathbb{N}.
```

In other words, $DF(x) = \mathcal{M}_y$, where $\mathcal{M}_y$ denotes the multiplication operator associated with $y$.
As a bounded linear operator acting on infinite sequences in $\ell^1_\mathbb{N}$, $\mathcal{M}_y$ can be visualized as a matrix with an infinite number of rows and columns:

```math
\mathcal{M}_y =
\begin{pmatrix}
y_0 & 0   & 0 & \cdots \\
y_1 & y_0 & 0 & \cdots \\
y_2 & y_1 & y_0 & \ddots \\
\vdots & \vdots & \vdots & \ddots
\end{pmatrix}.
```

Our goal is to find a numerical approximation $\bar{x}$ of the zero of $F$ as an element of the truncated space $\Pi_N \ell^1_\mathbb{N}$.
The general idea is to apply Newton's method on $\Pi_N \circ F \circ \Pi_N$, whose iterates are given by

```math
x \in \Pi_N\ell^1_\mathbb{N} \mapsto x - (\Pi_N DF(\Pi_N x) \Pi_N)^{-1} \Pi_N F(\Pi_N x) \in \Pi_N\ell^1_\mathbb{N}.
```

For our map $F$, these reduce to solving the following linear system in $\ell^1_\mathbb{N}$

```math
\Pi_N \mathcal{M}_y \Pi_N x = 1,
```

so that we seek $\bar{x}$ such that

```math
\Pi_N \mathcal{M}_y \Pi_N \bar{x} \approx 1.
```

If the truncation dimension $N$ is large enough we may hope that also $\Pi_{> N} F (\bar{x}) \approx 0$, since we expect the coefficients of $\bar{x}$ to decrease for moderately large $n$.
"""

# ╔═╡ 82c9f85d-f74f-4106-a5b8-35043e2d6119
y = Sequence(Taylor(1), [2.0, 1.0])

# ╔═╡ 1b058c77-7241-41ac-9a6e-84160b8d79da
N = 3

# ╔═╡ 82143c8d-7009-4537-975f-c832330caa9e
M = project(Multiplication(y), Taylor(N), Taylor(N))

# ╔═╡ a3124717-6741-4769-8fb0-2b489f13dc94
x_bar = M \ Sequence(Taylor(N), [1.0 ; zeros(N)])

# ╔═╡ d2d76ed8-911d-4247-9460-8c7df7669912
begin
	plot(LinRange(-1, 1, 101), t -> x_bar(t); linewidth = 3, label = "x_bar")
	xlims!(-1, 1)
	xlabel!("t")
end

# ╔═╡ 2346fe00-da58-4854-b476-6e3570941da2
md"""
We can see how close $F(\bar{x})$ is to zero.
"""

# ╔═╡ 1c61be31-8515-4b32-9d21-9a8910693f5a
F(x_bar, y)

# ╔═╡ 16e5b184-8d8c-43b9-a226-cdad95830b4f
md"""
Of course, in this example, we chose very few Taylor coefficients, so the residue is not as small as one would ideally aim for.
"""

# ╔═╡ 1d1277bf-2aa3-4c89-b34d-b2a2e5fc6fc5
md"## Constructing $A$"

# ╔═╡ 8a4513a5-9049-48df-ae9b-ea8c5d5f1076
md"""
For $T$ to be a contraction, we want $A$ to be a good approximate inverse of $\mathcal{M}_y$.
Recall that our numerical zero $\bar{x} \in \Pi_N \ell^1_\mathbb{N}$ satisfies $\bar{x} * y \approx 1$.
Thus, by taking

```math
A \overset{\text{def}}{=} \mathcal{M}_{\bar{x}},
```

we have that $A DF(\bar{x}) \approx I$ as desired:

```math
A [DF(\bar{x})] h = \bar{x} * y * h \approx 1 * h = h.
```

!!! note "Remark"
    In general, for an invertible function $w(t) = \sum_{n \ge 0} z_n t^n$, with its sequence of coefficients $z \in \ell^1_\mathbb{N}$, we have that $\mathcal{M}_z^{-1} = \mathcal{M}_{z^{-1}$}.
"""

# ╔═╡ 7edef822-8f17-4a17-8dc7-f476fd610a29
md"## Verifying the contraction"

# ╔═╡ 63e93962-935b-4c4a-8d10-04b43411b609
md"""
We now verify that $T(x) \overset{\text{def}}{=} x - A F(x)$ is a contraction around $\bar{x} \in \Pi_N \ell^1_\mathbb{N}$.
To rigorously compute the bounds $Y, Z_1, Z_2$, we enclose all relevant data into intervals.
"""

# ╔═╡ 5c57ed8e-7bb1-4947-8ebd-f0c59b7ae8d4
iy = interval(y)

# ╔═╡ bec0f423-5b0a-4dee-bedf-9e3b7fcf5d3d
ix_bar = interval(x_bar)

# ╔═╡ f1fa7997-684d-443b-a5ed-5377ebcd3aee
md"##### Computing $Y$"

# ╔═╡ e1dd0737-2db8-4cf3-9cec-bea25648201b
md"""
!!! lemma "Lemma"
	If $x \in \Pi_N \ell^1_\mathbb{N}$ and $y \in \Pi_{N'} \ell^1_\mathbb{N}$, then $x * y \in \Pi_{N+N'} \ell^1_\mathbb{N}$.
"""

# ╔═╡ 25df704a-3e01-45aa-a3de-ea2cb3f18504
md"""
```math
\| A F(\bar{x}) \|_1 = \| \bar{x} * (\bar{x} * y - 1) \|_1 \le Y.
```

Since $y \in \Pi_1 \ell^1_\mathbb{N}, \bar{x} \in \Pi_N \ell^1_\mathbb{N}$, we have that $A F(\bar{x}) \in \Pi_{2N + 1} \ell^1_\mathbb{N}$ and the above quantity is computable.
"""

# ╔═╡ d27d66a6-4717-4d6d-895f-cdebd649cebd
Y = norm(ix_bar * F(ix_bar, iy), 1)

# ╔═╡ 39ed4568-b971-4cf3-accc-1ad32a44ce3d
md"##### Computing $Z_1$"

# ╔═╡ d5b89e7d-01bf-433a-a6e0-302951343bcd
md"""
!!! lemma "Lemma"
	For all $x \in \ell^1_\mathbb{N}$ we have that $\| \mathcal{M}_x \|_{\mathscr{B}(\ell^1_\mathbb{N})} = \| x \|_1$.
"""

# ╔═╡ 5718bc5d-a86f-4907-8d5f-1707c8b51f58
md"""
```math
\| A DF(\bar{x}) - I \|_{\mathscr{B}(\ell^1_\mathbb{N})} = \| \bar{x} * y - 1 \|_1 = \| F(\bar{x}) \|_1 \le Z_1.
```

The above quantity is computable since $F(\bar{x}) \in \Pi_{N + 1} \ell^1_\mathbb{N}$.
"""

# ╔═╡ 9d468890-c155-442d-bcbf-6a33bc991f72
Z₁ = norm(F(ix_bar, iy), 1)

# ╔═╡ 3d9b3db1-ac82-434c-8521-ce5ce89f570d
md"##### Computing $Z_2$"

# ╔═╡ 2d273b10-b485-49aa-9d06-8485ac166dc8
md"""
```math
\|A(DF(x) - DF(\bar{x}))\|_{\mathscr{B}(\ell^1_\mathbb{N})} = 0, \qquad \forall x \in \ell^1_\mathbb{N}.
```

So $Z_2 = 0$ is a constant.
In other words, we have a Lipschitz control for $DF$ on the whole space $\mathbb{R}^3$, and we can freely choose $R = \infty$.
"""

# ╔═╡ cc5576ac-2894-48ef-a7a9-c4c16b589173
R = Inf

# ╔═╡ d02c07ee-77f2-4317-9692-558323e37b54
Z₂ = interval(0)

# ╔═╡ 6b63d5e1-048a-4a34-98a0-960d28a51a0b
md"##### Finishing the CAP"

# ╔═╡ 17cc5d87-c680-4ccd-a63d-dd354059bdc5
md"""
To determine the values of $r$ that satisfy the radii polynomial inequalities, we use the `interval_of_existence` function from RadiiPolynomial.
This function returns an interval of existence, within which all contained values satisfy the contraction conditions.
"""

# ╔═╡ e87f4cbe-e511-4c78-882b-f86961007932
ie = interval_of_existence(Y, Z₁, Z₂, R)

# ╔═╡ 7427f6e4-337c-4af8-ad59-dfd296b0a453
r = inf(ie)

# ╔═╡ c498e384-b1a2-40a5-bd72-edabdcd1fc4d
md"""
To guarantee that the fixed-point of $T$ is a zero of $F$, we must check that $A$ is injective.
This property comes as a by-product of the contraction argument.
Indeed, $\|I - A DF(\bar{x})\| \le Z_1 < 1$ ensures that $A DF(\bar{x})$ is invertible, which implies that $A$ is surjective.
Since $A [DF(\bar{x})] h = \bar{x} * y * h = y * \bar{x} * h = [DF(\bar{x})] A h$, for all $h \in \ell^1_\mathbb{N}$, it follows that $A$ must also be injective.

Therefore, we have proved the existence of a zero $\tilde{x}$ of $F$ such that

```math
\tilde{x} = \bar{x} + \gamma,
```

where $\gamma \in \ell^1_\mathbb{N}$ is not known explicitly, but satisfies $\| \gamma \|_1 \le r$.
Our choice of norm also gives us a $C^0$-error bound for the function

```math
\begin{aligned}
\sup_{t \in [-1,1]} | v(t)^{-1} - \sum_{n=0}^N \bar{x}_n t^n |
&= \sup_{t \in [-1,1]} | \left(\sum_{n \ge 0} y_n t^n \right)^{-1} - \sum_{n=0}^N \bar{x}_n t^n | \\
&= \sup_{t \in [-1,1]} | \sum_{n \ge 0} (\tilde{x}_n - \bar{x}_n) t^n | \\
&\le \| \tilde{x} - \bar{x} \|_1 \\
&\le r.
\end{aligned}
```
"""

# ╔═╡ 04417cef-c3d9-4adb-a449-7b26bbe20056
begin
	rigorous_eval(x, r, t) = interval(x(t), r; format = :midpoint)

	ts = LinRange(-1, 1, 101)

	plot([[interval(ts[i], ts[i+1]), rigorous_eval(ix_bar, r, interval(ts[i], ts[i+1]))] for i = 1:length(ts)-1];
		label = "", color = :palegreen3)

	plot!(ts, t -> inv(2+t);
		label = "1/(2+t)", color = :salmon, linewidth = 3)

	plot!(ts, t -> x_bar(t);
		label = "x_bar", color = :royalblue1, linewidth = 3, linestyle = :dash)

	xlims!(-1, 1)
	xlabel!("t")
end

# ╔═╡ d61da3a0-ac88-41a4-ae50-86a02ea43cbf
md"""
Lastly, this example is simple enough that we can compute by hand both the first 5 terms of $(2+t)^{-1}$

```math
\frac{1}{2+t} = \frac{1}{2} - \frac{t}{4} + \frac{t^2}{8} - \frac{t^3}{16} + \frac{t^4}{32} - \frac{t^5}{64} + O(t^6),
```

and compare with the result of our CAP: the approximate finite coefficients $\bar{x}$ and the error $r$.
"""

# ╔═╡ Cell order:
# ╟─ec5584b5-8b4d-4149-97c7-54fe0ee7fb8b
# ╠═51bbeeac-7ab4-4afc-9c2d-3f32b65c9b6d
# ╟─2af9b21a-2dcf-497f-bd7c-9ebbf80b9428
# ╟─11810156-e8c9-4145-af2a-93e0757a8cc9
# ╟─44fc05a5-a0c1-4428-adee-edb5b8b869e3
# ╟─9162d799-c8a8-44cd-8921-b6fa4c6f291f
# ╟─37827087-01ad-4fe7-bcb6-df2239ea0bca
# ╟─d794367f-d346-416f-b8f1-1f442f876e20
# ╟─f6da61f7-64a3-4c33-b9ba-7917ffb67b2a
# ╟─65cebd36-0474-42ee-b983-70c5f9c1614f
# ╟─2dc02de6-f250-44db-b347-7d4d65c61d89
# ╟─f9d764b9-024a-4845-98ad-c3e79a4a5cd9
# ╟─56bf8063-cc76-47ee-9074-f6a59466ec60
# ╟─3cb18cba-18e5-4e6f-8b16-eccd1dcef9fd
# ╟─7390523a-ecf1-4b7e-a3e6-cce50d43083b
# ╟─2a4448c5-e14a-4791-926f-172f3cb6e1f0
# ╟─faec7172-c732-470c-88a0-8de5d92075b0
# ╟─4290872b-923b-45dc-81b0-3e280a5a77c1
# ╟─8356c027-ead5-4073-9d0e-b25f41bf9543
# ╟─c11ce6e3-4e23-46bc-8bbd-217908d966af
# ╟─b6655a7c-32d0-473a-ac1d-ddb96e31de81
# ╟─266d8eb2-fc1f-42ca-802e-90a4a358b313
# ╟─0dab4ba7-4ad0-4ea8-a4c9-5892427969e7
# ╠═4f37b759-371e-4042-82d0-720523718021
# ╟─17013b31-0c79-4dd0-a6e3-5946528d69f6
# ╟─151f9f5f-65ab-4ccb-94a6-c530d924962f
# ╠═82c9f85d-f74f-4106-a5b8-35043e2d6119
# ╠═1b058c77-7241-41ac-9a6e-84160b8d79da
# ╠═82143c8d-7009-4537-975f-c832330caa9e
# ╠═a3124717-6741-4769-8fb0-2b489f13dc94
# ╟─d2d76ed8-911d-4247-9460-8c7df7669912
# ╟─2346fe00-da58-4854-b476-6e3570941da2
# ╠═1c61be31-8515-4b32-9d21-9a8910693f5a
# ╟─16e5b184-8d8c-43b9-a226-cdad95830b4f
# ╟─1d1277bf-2aa3-4c89-b34d-b2a2e5fc6fc5
# ╟─8a4513a5-9049-48df-ae9b-ea8c5d5f1076
# ╟─7edef822-8f17-4a17-8dc7-f476fd610a29
# ╟─63e93962-935b-4c4a-8d10-04b43411b609
# ╠═5c57ed8e-7bb1-4947-8ebd-f0c59b7ae8d4
# ╠═bec0f423-5b0a-4dee-bedf-9e3b7fcf5d3d
# ╟─f1fa7997-684d-443b-a5ed-5377ebcd3aee
# ╟─e1dd0737-2db8-4cf3-9cec-bea25648201b
# ╟─25df704a-3e01-45aa-a3de-ea2cb3f18504
# ╠═d27d66a6-4717-4d6d-895f-cdebd649cebd
# ╟─39ed4568-b971-4cf3-accc-1ad32a44ce3d
# ╟─d5b89e7d-01bf-433a-a6e0-302951343bcd
# ╟─5718bc5d-a86f-4907-8d5f-1707c8b51f58
# ╠═9d468890-c155-442d-bcbf-6a33bc991f72
# ╟─3d9b3db1-ac82-434c-8521-ce5ce89f570d
# ╟─2d273b10-b485-49aa-9d06-8485ac166dc8
# ╠═cc5576ac-2894-48ef-a7a9-c4c16b589173
# ╠═d02c07ee-77f2-4317-9692-558323e37b54
# ╟─6b63d5e1-048a-4a34-98a0-960d28a51a0b
# ╟─17cc5d87-c680-4ccd-a63d-dd354059bdc5
# ╠═e87f4cbe-e511-4c78-882b-f86961007932
# ╠═7427f6e4-337c-4af8-ad59-dfd296b0a453
# ╟─c498e384-b1a2-40a5-bd72-edabdcd1fc4d
# ╟─04417cef-c3d9-4adb-a449-7b26bbe20056
# ╟─d61da3a0-ac88-41a4-ae50-86a02ea43cbf
