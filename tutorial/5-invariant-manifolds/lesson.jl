### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Parameterization method"

using Markdown
using InteractiveUtils

# ╔═╡ c13e544f-5d8f-490c-ac94-50b6e3eba79b
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ 9e0f949b-ea89-4bc7-acd4-d9ad1575f38e
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
    using PlutoTeachingTools, PlutoUI, RadiiPolynomial, LinearAlgebra
end

# ╔═╡ be25c236-6e65-47c7-a36c-c0a6eb594f57
html"""
<style>
main {
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ 5482d7a2-df89-4ca9-b3f7-5fac7967f6e1
TableOfContents(title = "Table of Contents"; indent = true, depth = 4, aside = true)

# ╔═╡ e507fbe8-e81c-4386-b33b-b0482815f8c1
md"""
# Local stable manifolds of equilibria
"""

# ╔═╡ 5c7f76e1-ba35-4cb7-9fdc-cfbf05a9ab13
md"""
This part of the tutorial is concerned with the computation of local stable manifolds (the computation of local unstable manifolds is done similarly). Our focus is the development of efficient high order approximation methods --with validated error bounds--
for these manifolds.

To be more specific, consider a differential equation $\dot{x} = f(x)$ with associated flow $\varphi\colon \mathbb{R} \times \mathbb{R}^n \to \mathbb{R}^n$ and assume that $\tilde{x}$ is an equilibrium.
The *stable set* for $\tilde{x}$ is defined to be

$\begin{align}
W^s(\tilde{x}) = \{x\in\mathbb{R}^n : \lim_{t\to \infty} \varphi(t,x) = \tilde{x} \}.
\end{align}$

At this level of generality there is very little that can be said about $W^s(\tilde{x})$, however a natural step is to simplify the problem by restricting the dynamics to a neighborhood of $\tilde{x}$. Hence, let us define a more localized notion of a stable set. Given a neighborhood $U$ of $\tilde{x}$ the associated *local stable manifold* is given by

$\begin{align}
W^s_{\text{loc}}(\tilde{x}) = W^s(\tilde{x},U) := \{x \in W^s(\tilde{x}) \mid \varphi([0,\infty),x) \subset U \}.
\end{align}$"""

# ╔═╡ 51622c96-c674-4079-ad9d-3554b6562ed4
md"""
## The Parameterization Method
"""

# ╔═╡ a707c66a-8e9f-454f-a0ad-bcc4d628d82d
md"""
### Assumptions
"""

# ╔═╡ 1e145c27-90cd-43f1-b8b2-4ae19b8b9ec3
md"""

**For sake of simplicity of the presentation**, we make the following assumptions:

-  $\tilde{x}$ is an equilibrium solution of $\dot{x}=f(x)$
-  there are exactly $m$ stable eigenvalues $\{\lambda_1, \ldots, \lambda_m \}$, they are real and there exist $m$ linearly independent associated eigenvector $\xi_i \in \mathbb{R}^n$

Denote
$\begin{align}
\Lambda =
\begin{pmatrix}
 \lambda_1 &  \ldots & 0  \\
 \vdots & \ddots  & \vdots  \\
 0  & \ldots  & \lambda_m
\end{pmatrix} \in M_m(\mathbb{R})
\quad \text{and} \quad A_0 = [\xi_1 | \ldots | \xi_m]
\end{align}$
an $n \times m$ matrix whose columns are the associated (stable) eigenvectors.

The linearized equation restricted to the stable subspace is $\dot{\theta} = \Lambda \theta$ ($\theta \in \mathbb{R}^m$) with associated flow $e^{\Lambda t}$.
"""

# ╔═╡ aa2f4387-e3d9-48b1-bee1-af8c80dade1f
md"""
### The conjugacy relation
"""

# ╔═╡ c60b9061-3143-46b4-997f-c5c072b81da5
md"""
Denote by $B^m := \{ \theta = (\theta_1,\dots,\theta_m) : |\theta_j| \le 1 \}$ the unit closed $\ell^\infty$ ball in $\mathbb{R}^m$.

To solve for the parameterization $P:B^m \to \mathbb{R}^n$ of the local stable manifold, we impose

(1) $P(0) = \tilde{x}$ (equilibria are mapped to one another)

(2) $DP(0) = A_0$ ($P$ is tangent to the stable eigenspace $\mathbb{E}^s$ at $0$})

(3) $\varphi (t, P(\theta)) = P( e^{\Lambda t} \theta)$, for all $\theta \in B^m$ (**conjugacy relation**)
"""

# ╔═╡ e4163286-ceda-4277-8da9-5a9394b15e54
Foldable("""Why such parameterization provides a local stable manifold?""", md"""
We get that $P(B^m) \subset W^s_{\text{loc}}(\tilde{x})$, since for $\theta \in B^m$,

$\begin{align}
\lim_{t \to \infty} \varphi (t, P(\theta) ) = \lim_{t \to \infty}
	P( e^{\Lambda t} \theta)
			= P\left(\lim_{t \to \infty} e^{\Lambda t} \theta \right)
		    = P(0)
		    = \tilde{x}.
\end{align}$
""")

# ╔═╡ 6d00a202-4d6b-497e-b423-76e3323c4f24
md"""
### The parameterization lemma
"""

# ╔═╡ 07dcb89c-9926-432e-94d4-44cb98c2c8a1
md"""

*Remark.* The utility of (3) above is limited by the appearance of the flow $\varphi$, which is only known implicitly. Let us introduce a more practical infinitesimal version of the conjugacy relation, which requires the following important result."""

# ╔═╡ 68c49689-e740-4ddb-bc67-6f811e8b891d
Markdown.MD(Markdown.Admonition("tip", "Parameterization Lemma",[md"""
Let $P \colon B^m  \to \mathbb{R}^n$ be a function satisfying (1) and (2). Then $P$ satisfies (3) $\Longleftrightarrow$ $P$ solves the *partial* differential equation

$\begin{align}
\text{(4)} \quad DP(\theta) \Lambda \theta = f(P(\theta)), \quad \theta \in B^m.
\end{align}$"""]))

# ╔═╡ 4d46da29-6ec9-4d0b-b28e-843c1774b2c4
Foldable("""Proof""", md"""
($\Longleftarrow$) Assume that the assumption (4) above holds and let $\theta \in B^m$. Let $\gamma(t) := P(e^{\Lambda t} \theta)$. Then, $\gamma(0) = P(\theta)$ and

$\begin{align}
\gamma'(t) = \frac{d}{dt}  P\left( e^{\Lambda t} \theta \right)
   	= D P\left( e^{\Lambda t} \theta \right) \Lambda e^{\Lambda t} \theta
 	= f\left( P\left( e^{\Lambda t} \theta \right)\right)
 	= f ( \gamma(t)),
\end{align}$

that is $\gamma$ solves the IVP $x'=f(x), x(0) = P(\theta)$. By uniqueness, $\varphi( t, P(\theta)) = \gamma(t) = P(e^{\Lambda t} \theta)$.

($\Longrightarrow$) Suppose that $P$ solves (3). Differentiate on both sides with respect to $t$

$\begin{align}
f(\varphi(t,P(\theta))) = DP( e^{\Lambda t} \theta) \Lambda e^{\Lambda t} \theta, \quad \text{and take } t=0. \quad \square
\end{align}$
""")

# ╔═╡ 4e9479fe-2ddd-4ea5-a8dd-93ef889aef33
md"""
From the Parameterization Lemma, to solve for the parameterization $P:B^m \to \mathbb{R}^n$, we impose

(1) $P(0) = \tilde{x}$ (equilibria are mapped to one another)

(2) $DP(0) = A_0$ ($P$ is tangent to the stable eigenspace $\mathbb{E}^s$ at $0$)

(4) $DP(\theta) \Lambda \theta = f(P(\theta))$, $\forall~\theta \in B^m$ (invariance equation)

The invariance equation can be written more explicitly, for $\theta = (\theta_1, \ldots, \theta_m) \in B^m$, as

$\begin{align}
(4) \quad \lambda_1 \theta_1 \frac{\partial}{\partial \theta_1} P(\theta_1, \ldots, \theta_m)
+ \ldots + \lambda_m \theta_m \frac{\partial}{\partial \theta_m} P(\theta_1, \ldots, \theta_m)
= f(P(\theta_1, \ldots, \theta_m))
\end{align}$
"""

# ╔═╡ a85db154-e650-40bd-8191-175a28a5de02
md"""
### Taylor series representation
"""

# ╔═╡ ac318edc-cca6-456a-b4ac-7ece084089c5
md"""
We use a Taylor series representation of the form

$\begin{align}
P(\theta) = \sum_{|\alpha|=0}^\infty \begin{pmatrix} (a_1)_\alpha \\ (a_2)_\alpha \\ \vdots \\ (a_n)_\alpha \end{pmatrix}
\theta^\alpha
\end{align}$

where $\alpha = (\alpha_1,\dots,\alpha_m) \in \mathbb{N}^m$, $|\alpha| = \alpha_1 + \dots + \alpha_m$, $\theta=(\theta_1,\dots,\theta_m)$ and $\theta^\alpha = \theta_1^{\alpha_1} \cdots \theta_m^{\alpha_m}$.

Denote $a = (a_1,\dots,a_n)$, where each $a_i = ((a_i)_\alpha)_{|\alpha| \ge 0}$. Plugging the above Taylor series in equations (1), (2) and (4):

(1*) $\begin{pmatrix} (a_1)_{(0,\dots,0)} \\ \vdots \\ (a_n)_{(0,\dots,0)} \end{pmatrix} = \tilde{x}$, (for $|\alpha|=0$: the equilibrium solution)

(2*) $\begin{pmatrix} (a_1)_{e_j} \\ \vdots \\ (a_n)_{e_j} \end{pmatrix} = \xi_j$, $j=1,\dots,m$
 (for $|\alpha|=1$: the $m$ stable eigenvectors)

(4*) $(\alpha \cdot \lambda) (a_i)_{\alpha} = (\phi_i(a))_{\alpha}$, $i=1,\dots,n$ (for $|\alpha| \ge 2$)

where

$\begin{align}
\alpha \cdot \lambda = \alpha_1 \lambda_1 +\cdots + \alpha_m \lambda_m \quad \text{and} \quad f_i(P(\theta))= \sum_{|\alpha| \ge 0} (\phi_i(a))_\alpha \theta^\alpha
\end{align}$
"""

# ╔═╡ da95dea7-9ec5-46f4-beea-9a22cf5ca6ee
md"""
### The Zero Finding Problem, Banach Space and Norm
"""

# ╔═╡ 510a3921-abe3-414f-ad17-c165415d8e40
md"""
For each $i=1,\dots,n$, define $F_i(a)=\{ (F_i(a))_\alpha \}_{|\alpha| \ge 0}$ component-wise by

$\begin{align}
(F_i(a))_\alpha := \begin{cases}
(a_i)_{(0,\dots,0)}  - \tilde{x}_i,&  |\alpha| = 0, \\
(a_i)_{e_j}  - (\xi_{j})_i, & |\alpha| = 1, {\rm~with~} \alpha=e_j {\rm~for~}j=1,\dots,m,  \\
(a_i)_{\alpha} - \frac{1}{\alpha \cdot \lambda}(\phi_i(a))_{\alpha}, & |\alpha| \ge 2.
\end{cases}
\end{align}$

For each $i=1,\dots,n$, let $c_i = \{ (c_i)_\alpha\}_{|\alpha| \ge 0}$ component-wise as

$\begin{align}
(c_i)_\alpha =
\begin{cases}
\tilde{x}_i, &  {\rm~if~} |\alpha| = 0, \\
\xi_{j}, &  {\rm~for~} |\alpha| = 1 {\rm~with~} \alpha=e_j {\rm~for~} j=1,\dots,m,  \\
0, & {\rm~if~} |\alpha| \ge 2
\end{cases}
\end{align}$

and let $c := (c_1,\dots,c_n)$. Define the diagonal operator $K$ acting on a sequence $b = \{ b_\alpha\}_{|\alpha| \ge 0}$ as

$\begin{align}
(K b)_\alpha =
\begin{cases}
0, &   {\rm~if~} |\alpha| = 0, \\
0, &   {\rm~if~} |\alpha| = 1,  \\
\frac{1}{\alpha \cdot \lambda} b_\alpha, &  {\rm~if~} |\alpha| \ge 2.
\end{cases}
\end{align}$

Define the operator $\mathcal{K}$ acting on $a=(a_1,\dots,a_n)$ as
$\begin{align}
\mathcal{K} a = (K a_1,K a_2,\dots,K a_n).
\end{align}$

With the above notation, for each $i=1,\dots,n$,

$\begin{align}
F_i(a) = a_i - K \phi_i(a) - c_i.
\end{align}$

Letting $F(a) = (F_1(a),\dots,F_n(a))$ and $\phi(a) = (\phi_1(a),\dots,\phi_n(a))$, one can densely write

$\begin{align}
F(a) = a - \mathcal{K} \phi(a) - c
\end{align}$

Define the Banach space $\mathcal{X} = (\ell^1)^n$, where

$\begin{align}
\ell^1 := \left\{ c = (c_\alpha)_{|\alpha| \ge 0} : c_\alpha \in \mathbb{R}, \quad \|c\|_1 := \sum_{|\alpha|=0}^\infty |c_\alpha| < \infty  \right\}.
\end{align}$

The norm on $\mathcal{X}$ is given by

$\begin{align}
\|a\|_{\mathcal{X}} = \sum_{i=1}^n \| a_i \|_1
\end{align}$
"""

# ╔═╡ 32179842-b1fc-4cbb-bb98-a39b6ebce9ee
md"""
###### Code that implements the operator $K$
"""

# ╔═╡ f1bb1b2a-c2b5-48a3-897a-bcee907333e3
function K(λ, N)
    op = zeros(typeof(λ), Taylor(N), Taylor(N))
    for n = 2:N
        op[n,n] = inv(λ * n)
    end
    return op
end

# ╔═╡ ccb1e680-934b-46b0-ad8f-4bd9516010c6
K(-2.0,5)

# ╔═╡ 3cfbcab7-c3cc-45ce-99cc-b9e800bec98b
Markdown.MD(Markdown.Admonition("note", "Exercise (Compact operator)",[md"""
Show that the operator $\mathcal{K}:\mathcal{X} \to \mathcal{X}$ is compact.
"""]))

# ╔═╡ b9ad918e-6b48-4c14-a2cc-ead87c667070
md"""
We now formalize how a solution to the above zero finding problem leads to local stable manifolds of equilibria of our differential equation.
"""

# ╔═╡ 0548b47b-cc62-4730-8293-26b2afe3f793
Markdown.MD(Markdown.Admonition("lemma", "Lemma (A solution to the zero finding problem yields a parameterization of the local stable manifold)",[md""" Assume that
-  $\tilde{x}$ is an equilibrium solution of $\dot{x}=f(x)$;
-  there are exactly $m$ stable eigenvalues $\{\lambda_1, \ldots, \lambda_m \}$, they are real and there exist $m$ linearly independent associated eigenvector $\xi_i \in \mathbb{R}^n$.

If $a \in \mathcal{X}$ is a solution to $F(a)=0$, then the corresponding Taylor series expansion

$\begin{align}
P(\theta) = \sum_{|\alpha|=0}^\infty \begin{pmatrix} (a_1)_\alpha \\ (a_2)_\alpha \\ \vdots \\ (a_n)_\alpha \end{pmatrix}
\theta^\alpha
\end{align}$

converges on $B^m := \{ \theta = (\theta_1,\dots,\theta_m) : |\theta_j| \le 1 \}$ and the function $P:B^m \to \mathbb{R}^n$ provides a parameterization of the local stable manifold $W^s_{\text{loc}}(\tilde{x})$, that is

$\begin{align}
P(B^m) = W^s_{\text{loc}}(\tilde{x}).
\end{align}$
"""]))

# ╔═╡ ce8ea7a5-dea3-4cda-9dea-eab93c54b24d
md"""
Rather than presenting the bounds for the Newton-Kantorovich Theorem in general, let us rather focus on an example. Specifically, we compute a rigorous parameterization of a one-dimensional local stable manifold of an equilibrium solution in the Lorenz system.
"""

# ╔═╡ Cell order:
# ╟─c13e544f-5d8f-490c-ac94-50b6e3eba79b
# ╠═9e0f949b-ea89-4bc7-acd4-d9ad1575f38e
# ╟─be25c236-6e65-47c7-a36c-c0a6eb594f57
# ╟─5482d7a2-df89-4ca9-b3f7-5fac7967f6e1
# ╟─e507fbe8-e81c-4386-b33b-b0482815f8c1
# ╟─5c7f76e1-ba35-4cb7-9fdc-cfbf05a9ab13
# ╟─51622c96-c674-4079-ad9d-3554b6562ed4
# ╟─a707c66a-8e9f-454f-a0ad-bcc4d628d82d
# ╟─1e145c27-90cd-43f1-b8b2-4ae19b8b9ec3
# ╟─aa2f4387-e3d9-48b1-bee1-af8c80dade1f
# ╟─c60b9061-3143-46b4-997f-c5c072b81da5
# ╟─e4163286-ceda-4277-8da9-5a9394b15e54
# ╟─6d00a202-4d6b-497e-b423-76e3323c4f24
# ╟─07dcb89c-9926-432e-94d4-44cb98c2c8a1
# ╟─68c49689-e740-4ddb-bc67-6f811e8b891d
# ╟─4d46da29-6ec9-4d0b-b28e-843c1774b2c4
# ╟─4e9479fe-2ddd-4ea5-a8dd-93ef889aef33
# ╟─a85db154-e650-40bd-8191-175a28a5de02
# ╟─ac318edc-cca6-456a-b4ac-7ece084089c5
# ╟─da95dea7-9ec5-46f4-beea-9a22cf5ca6ee
# ╟─510a3921-abe3-414f-ad17-c165415d8e40
# ╟─32179842-b1fc-4cbb-bb98-a39b6ebce9ee
# ╠═f1bb1b2a-c2b5-48a3-897a-bcee907333e3
# ╠═ccb1e680-934b-46b0-ad8f-4bd9516010c6
# ╟─3cfbcab7-c3cc-45ce-99cc-b9e800bec98b
# ╟─b9ad918e-6b48-4c14-a2cc-ead87c667070
# ╟─0548b47b-cc62-4730-8293-26b2afe3f793
# ╟─ce8ea7a5-dea3-4cda-9dea-eab93c54b24d
