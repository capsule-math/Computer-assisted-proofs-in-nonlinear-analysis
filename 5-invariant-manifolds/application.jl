### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Lorenz system"

using Markdown
using InteractiveUtils

# ╔═╡ 97c19551-3248-47f6-a1cb-bc4fb732e6df
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ 97df61b8-00b4-4c51-bbc0-b09df0e4cbe1
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
    using PlutoTeachingTools, PlutoUI, RadiiPolynomial, LinearAlgebra
end

# ╔═╡ 7fc40507-eda3-474d-a454-04e9173a7adb
html"""
<style>
main {
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ 8f59e349-74a6-4a6f-9a47-6ca26d187aeb
html"""
<style>
main {
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ 4412bd33-8b3f-41de-9e1a-a8f1187c793f
TableOfContents(title = "Table of Contents"; indent = true, depth = 4, aside = true)

# ╔═╡ d9604e75-0e9e-44b6-bf2e-885f01a96bdc
md"""
## A one-dimensional local stable manifold in the Lorenz system
"""

# ╔═╡ 62719aa9-bfda-4913-8b74-f3694bfc33c0
md"""
The Lorenz system is given by

$\begin{align}
\dot x = f(x) := \begin{pmatrix}
\sigma (x_2 - x_1) \\
\rho x_1 - x_2 - x_1 x_3 \\
- \beta x_3 + x_1 x_2
\end{pmatrix}
\end{align}$

where $\sigma,\rho, \beta>0$ are parameters.

It is straightforward (**Exercise**) to verify that if $\rho>1$, there are two nontrivial equilibria given by
$\left\{
\begin{pmatrix} -\sqrt{\beta(\rho-1)} \\ -\sqrt{\beta(\rho-1)} \\ \rho-1 \end{pmatrix},
\quad
\begin{pmatrix} \sqrt{\beta(\rho-1)} \\ \sqrt{\beta(\rho-1)} \\ \rho-1 \end{pmatrix} \right\}$.
"""

# ╔═╡ fc381b29-6671-4364-a9df-b35db07af45b
md"""
We make the following **assumptions**:
-  Assume that $\rho>1$, and let us choose the equilibrium solution to be
$\begin{align}
\tilde x := \begin{pmatrix} \sqrt{\beta(\rho-1)} \\ \sqrt{\beta(\rho-1)} \\ \rho-1 \end{pmatrix}.
\end{align}$
-  Assume that there is a unique stable real eigenvalue $\lambda<0$ with associated eigenvector $\xi \in \mathbb{R}^3$. Hence, the dimension of $W^s_{\text{loc}}(\tilde{x})$ is one, that is $m=1$.
"""

# ╔═╡ 194ee83b-7442-42f6-862a-2f5bd1c8b710
Markdown.MD(Markdown.Admonition("warning", "Remark",[md""" For $\rho>24.75$, $\sigma=10$ and $\beta=8/3$, $Df(\tilde{x})$ has one real negative eigenvalues and two complex conjugate eigenvalues with positive real parts (**Exercise**). In this case, the local stable manifold $W^s_{\text{loc}}(\tilde{x})$ is one-dimensional.
"""]))

# ╔═╡ 9829987c-2bf5-4b33-8ae1-6b80d43ed978
md"""
Under the above assumptions, $m=1$, and the map $F$ is given by

$\begin{align}
F(a) = a - \mathcal{K} \phi(a) - c, \quad \quad \quad {\rm where } \quad \phi(a):=\begin{pmatrix}
\sigma (a_2 - a_1) \\
\rho a_1 - a_2 - a_1*a_3 \\
- \beta a_3 + a_1*a_2
\end{pmatrix}
\end{align}$

where $a_i*a_j$ is the **Cauchy product** defined by

$\begin{align}
(a_i*a_j)_n = \sum_{k=0}^n (a_i)_k (a_j)_{n-k}
\end{align}$

The following result will be useful when performing the nonlinear analysis.
"""

# ╔═╡ 1f8fba92-d023-4444-b4e3-18c751f678c5
Markdown.MD(Markdown.Admonition("note", "Exercise (Banach algebra)",[md"""
Prove that for all $a_1,a_2 \in \ell^1$, $\|a_1*a_2\|_1 \le \|a_1\|_1 \|a_2\|_1$.
"""]))

# ╔═╡ 11aa1d55-1a2c-48cf-a7f5-30ff16e2e92a
md"""
Note that in the case of a one-dimensional manifold, the sequences are indexed over a single index running over $\mathbb{N}$, as opposed to the more general multi-index case of an $m$-dimensional stable manifold. In the 1D case,

$\begin{align}
\ell^1 = \left\{ c = (c_n)_{n \ge 0} : c_n \in \mathbb{R}, \quad \|c\|_1 := \sum_{n=0}^\infty |c_n| < \infty  \right\}
\end{align}$
"""

# ╔═╡ 668058ef-e66b-4fd4-ad67-2a8389589e86
Markdown.MD(Markdown.Admonition("note", "Exercise",[md""" Show that for all $a = (a_1,a_2,a_3) \in \mathcal{X} := (\ell^1)^3$, $F(a) \in \mathcal{X}$. That shows that $F:\mathcal{X} \to \mathcal{X}$.
"""]))

# ╔═╡ 0562b884-0653-4d5b-8a66-e35467e6693e
md"""
#### Evaluating the operator norm
"""

# ╔═╡ a47a4e2a-371e-4244-8574-493f00d456df
Markdown.MD(Markdown.Admonition("note", "Exercise (characterization of operator norms on product of sequence spaces)",[md"""

Consider the Banach space $\mathcal{X}=(\ell^1)^n$ with norm

$\begin{align}
\|a\|_{\mathcal{X}} = \sum_{i=1}^n \| a_i \|_1.
\end{align}$

Let $\Gamma \in B(\mathcal{X})$ be a bounded linear operator acting on $\mathcal{X} = (\ell^1)^n$ and represented by

$\begin{align}
\Gamma = \begin{pmatrix}
\Gamma_{1,1} & \cdots & \Gamma_{1,n} \\
\vdots & \ddots & \vdots \\
\Gamma_{n,1} & \cdots & \Gamma_{n,n}
\end{pmatrix}
\end{align}$

where each $\Gamma_{i,j}$ is a bounded linear operator on $\ell^1$. More explicitly, $\Gamma$ acts on an element $a = (a_1,\dots,a_n) \in \mathcal{X}$ as

$\begin{align}
\Gamma a = \left( (\Gamma a)_1,\dots,(\Gamma a)_n \right) = \left( \sum_{j=1}^n \Gamma_{1,j} a_j , \dots, \sum_{j=1}^n \Gamma_{n,j} a_j \right),
\end{align}$

where$(\Gamma_{i,j} a_j)_{n_1} = \sum_{n_2 \ge 0 } (\Gamma_{i,j})_{n_1,n_2} (a_j)_{n_2}$. Then, show that

$\begin{align}
\boxed{
  \|\Gamma\|_{B(\mathcal{X})} = \max_{j=1,\dots,n} \sum_{i=1}^n \| \Gamma_{i,j} \|_{B(\ell^1)}
}
\end{align}$

where

$\begin{align}
  \|\Gamma_{i,j} \|_{B(\ell^1)} = \sup_{n_2 \ge 0} \sum_{n_1 \ge 0} |(\Gamma_{i,j})_{n_1,n_2}| .
\end{align}$
"""]))

# ╔═╡ beca7238-ca9b-409a-8d5a-240ca720f58a
md"""
Here is an example of how to compute the operator norm of $\pi^{\le 5}\mathcal{K}$
"""

# ╔═╡ c8c8856f-fe43-4b77-aad9-57723b6dcf9c
function K(λ, N)
    op = zeros(typeof(λ), Taylor(N), Taylor(N))
    for n = 2:N
        op[n,n] = inv(λ * n)
    end
    return op
end

# ╔═╡ 5e6482db-9b76-4bf0-bc58-0a6f779c4556
begin
	mat = K(-2.0, 5)
	M = zeros(Taylor(5)^3, Taylor(5)^3) # The operator K
	component(M, 1, 1) .= mat
	component(M, 2, 2) .= mat
	component(M, 3, 3) .= mat
	ell1 = Ell1() #no entry means that the weight is 1
	X = NormedCartesianSpace((ell1, ell1, ell1), ell1) # each space has l1 and the outer norm is l1
	opnorm(M, X) # this is equivalent to opnorm(M, 1)
end

# ╔═╡ 8c0ac2a2-b3db-455a-8a6d-3de6728b61d3
Markdown.MD(Markdown.Admonition("note", "Exercise",[md"""
Show theoretically that $\|\mathcal{K}\|_{B(\mathcal{X})} \le \frac{1}{2|\lambda|}$
"""]))

# ╔═╡ 7cfbe32b-d677-4611-97d9-ed775c60e48e
md"""
#### Finite dimensional projection
"""

# ╔═╡ 13f3b4e4-0409-49a2-a544-0c34ff9289db
md"""
Computing an approximate solution to $F=0$ requires first considering a finite dimensional projection. We need the projection operator on a finite number of Taylor coefficients ($N+1$ coefficients) and its complement:

$\begin{align}
(\pi^{\leq N} b)_n := \begin{cases}
b_n, & n=0,\dots,N, \\
0, & n > N,
\end{cases}
& \qquad\qquad
(\pi^{>N} b)_n := \begin{cases}
0, & n=0,\dots,N, \\
b_n, & n > N.
\end{cases}
\end{align}$

Now, given three sequences $a_1,a_2,a_3$ and for $a=(a_1,a_2,a_3)$, consider the finite dimensional projection on a finite number of Taylor coefficients ($3(N+1)$ coefficients) and its complement:

$\begin{align}
\pi^{\leq N} a = \left( \pi^{\leq N} a_1,\pi^{\leq N} a_2,\pi^{\leq N} a_3 \right) \qquad {\rm and} \qquad
\pi^{> N} a = \left( \pi^{> N} a_1,\pi^{> N} a_2,\pi^{> N} a_3 \right)
\end{align}$

so that $a=\pi^{\leq N}a+\pi^{>N}a$ and $\pi^{\leq N}\pi^{>N}=0=\pi^{>N}\pi^{\leq N}$.

The range $\pi^{\leq N} \mathcal{X}$ is finite dimensional, and the restriction of $\pi^{\leq N} F$ to $\pi^{\leq N} \mathcal{X}$, which we denote by $F^{\leq N}$, is what we will work with in the computer. An implementation is given below, as well as an implementation of its derivative.
"""

# ╔═╡ 40156798-59d5-4afe-a054-29c753a51e46
md"""
## Numerics
"""

# ╔═╡ 42ad70a3-114e-4cd3-8dd1-2429ac80bfc2
md"""
##### Implementation code for the finite dimensional projection $F^{\leq N}$
"""

# ╔═╡ cdbbb2e1-56bc-4beb-93d3-2db1a2aa71d5
function F(a, x, λ, ξ, par, N)
    a₁, a₂, a₃ = eachcomponent(a)

    σ, ρ, β = par

    F₁ = project(a₁ - K(λ, N) * (σ*(a₂ - a₁)), Taylor(N))
    F₁[0] -= x[1]
    F₁[1] -= ξ[1]

    F₂ = project(a₂ - K(λ, N) * (a₁*(ρ - a₃) - a₂), Taylor(N))
    F₂[0] -= x[2]
    F₂[1] -= ξ[2]

    F₃ = project(a₃ - K(λ, N) * (a₁*a₂ - β*a₃), Taylor(N))
    F₃[0] -= x[3]
    F₃[1] -= ξ[3]

    return Sequence(Taylor(N)^3, [coefficients(F₁) ; coefficients(F₂) ; coefficients(F₃)])
end

# ╔═╡ aaf278f6-b8ee-4094-9173-ffde2cb0b461
md"""
##### Implementation code for the Jacobian matrix $DF^{\leq N}$
"""

# ╔═╡ d77a192e-4be0-4987-8023-8d91b22d8d9f
function DF(a, λ, par, N)
    a₁, a₂, a₃ = eachcomponent(a)

    σ, ρ, β = par

    D₁F₁ = I - K(λ, N) * (-σ)
    D₁F₁[0,:] .= 0; D₁F₁[1,:] .= 0; D₁F₁[0,0] = D₁F₁[1,1] = 1
    D₂F₁ =   - K(λ, N) * σ
    D₃F₁ = zeros(Taylor(N), Taylor(N))

    D₁F₂ =   - K(λ, N) * project(Multiplication(ρ - a₃), Taylor(N), Taylor(N))
    D₂F₂ = I - K(λ, N) * (-1)
    D₂F₂[0,:] .= 0; D₂F₂[1,:] .= 0; D₂F₂[0,0] = D₂F₂[1,1] = 1
    D₃F₂ =   - K(λ, N) * project(Multiplication(-a₁),    Taylor(N), Taylor(N))

    D₁F₃ =   - K(λ, N) * project(Multiplication(a₂), Taylor(N), Taylor(N))
    D₂F₃ =   - K(λ, N) * project(Multiplication(a₁), Taylor(N), Taylor(N))
    D₃F₃ = I - K(λ, N) * (-β)
    D₃F₃[0,:] .= 0; D₃F₃[1,:] .= 0; D₃F₃[0,0] = D₃F₃[1,1] = 1

    return LinearOperator(Taylor(N)^3, Taylor(N)^3,
        [coefficients(D₁F₁) coefficients(D₂F₁) coefficients(D₃F₁);
         coefficients(D₁F₂) coefficients(D₂F₂) coefficients(D₃F₂);
         coefficients(D₁F₃) coefficients(D₂F₃) coefficients(D₃F₃)])
end

# ╔═╡ 7410e377-715c-4464-a83d-d7cd48c4d1ae
md"""
From now on, we fix the parameters in the Lorenz system to be the *classical* ones, namely
$\begin{align}
\rho=28, \quad \sigma=10 \quad \text{and} \quad \beta=8/3.
\end{align}$
"""

# ╔═╡ fbeac0c7-8fcb-4e4b-9662-9f1d48b23fb3
begin
σ, ρ, β = 10.0, 28.0, 8/3
par = [σ, ρ, β]
tσ, tρ, tβ = interval(10.0), interval(28.0), interval(8)/interval(3)
tpar = [tσ, tρ, tβ]
end

# ╔═╡ ae88c826-f486-4786-9447-affcc466d1b5
md"""
For these parameters, there is a unique real negative eigenvalue $\lambda$ of $Df(\tilde{x})$ with associated eigenvector $\xi \in \mathbb{R}^3$ (see the above Remark).
"""

# ╔═╡ 76aa085c-d7ab-4e2f-9648-89b3eebe11cd
md"""
### Rigorous computation of the stable eigenvalue $\lambda$ and eigenvector $\xi \in \mathbb{R}^3$
"""

# ╔═╡ 585b99f6-fad0-44b1-b1ca-f7aa6f8f52f8
md"""
While we could work hard and compute the eigenpair $(\lambda,\xi)$ by hand, we resort to a computer-assisted proof to obtain a rigorous enclosure. Here is a code that performs this task.
"""

# ╔═╡ 400fe5e0-0dfd-4d82-baf9-b66d931133a8
md"""
##### Numerical approximation and rigorous enclosure of the equilibrium solution $\tilde{x}$
"""

# ╔═╡ e0f75d0e-555d-4df7-8077-9871a6d7b197
begin
bx = [sqrt(β*(ρ-1)), sqrt(β*(ρ-1)), ρ-1] #Numerical approximation
tx = [sqrt(tβ*(tρ-interval(1))), sqrt(tβ*(tρ-interval(1))), tρ-interval(1)] #Rigorous enclosure
end

# ╔═╡ deb48c67-19f1-481e-b2de-1765e449e6ae
md"""
##### Jacobian matrix for the Lorenz system
"""

# ╔═╡ c106822b-2c0c-4759-84fb-9ddba91fe295
function Df(u, par)
    σ, ρ, β = par
    u₁, u₂, u₃ = u
    return [-σ     σ         zero(u₃)
             ρ-u₃ -one(u₂)  -u₁
             u₂    u₁       -β]
end

# ╔═╡ f75c26b4-8e6f-4520-b1a8-f30c20b9a194
md"""
##### Computation of an approximate eigenpair
"""

# ╔═╡ 05a4dd64-b5c9-4014-a18c-f26938e02ebd
begin
approx_eigenvalues, approx_eigenvectors = eigen(Df(bx, par))
bλ, bξ = real(approx_eigenvalues[1]), real(approx_eigenvectors[:,1])
end

# ╔═╡ 3cb601bf-43ea-430c-950d-402b0929e515
md"""
##### Zero finding problem for an eigenpair (see Part I of the tutorial)
"""

# ╔═╡ 2a76f46f-c203-4028-80ba-e7bc1f481076
function F_eig(X, M, ū)
	λ = X[1]
	u = X[2:end]
	return Sequence(
		[sum(u .* conj(ū)) - 1;
		 M*u - λ*u])
end

# ╔═╡ 8a802d5f-e2f5-4e03-b0d9-b7cd0c081d2c
md"""
##### The Jacobian of the zero finding map for an eigenpair
"""

# ╔═╡ 5ef06ffe-d2f3-4cd2-946a-df8fd78f863f
function DF_eig(X, M, ū)
	λ = X[1]
	u = X[2:end]
	n = length(X) - 1
	temp = zeros(eltype(X), n+1, n+1)
	temp[2:end,1] = -u
	temp[1,2:end] = ū'
	temp[2:end,2:end] = M - λ*I
	return LinearOperator(temp)
end

# ╔═╡ e4a015b0-e7e6-4841-b97f-4a845f462196
md"""
##### A Newton-Kantorovich Theorem to obtain a computer-assisted proof for an eigenpair
"""

# ╔═╡ c3d0d2ae-5711-451b-946d-6cdca682fe06
function validate_eigenpair(λ, u, M)
	X = [λ ; u]
	ū = u
	iA = interval.(inv(DF_eig(X, mid.(M), ū)))
	iX = interval.(X)
	Y = norm(iA * F_eig(iX, M, ū), 1)
	Z₁ = opnorm(I - iA * DF_eig(iX, M, ū), 1)
	Z₂ = opnorm(iA, 1)
	return interval_of_existence(Y, Z₁, Z₂, Inf)
end

# ╔═╡ 4ff1b050-3036-4e0b-b9db-b07b4d17327c
md"""
###### The rigorous enclosure of the eigenpair
"""

# ╔═╡ ec3cd93f-877d-419e-a608-cef748306927
begin
err_eig = inf(validate_eigenpair(bλ, bξ, Df(tx, tpar)))
tλ = interval(bλ, err_eig; format = :midpoint)
tξ = interval.(bξ, err_eig; format = :midpoint)
end

# ╔═╡ 6e715b86-4604-49c3-8df1-929921e38f9a
md"""
##### We introduce a rescaling for the eigenvector (will impact the decay of the Taylor coefficients)
"""

# ╔═╡ ffe5218b-88ca-42e0-8a0d-074ab2e6d898
begin
eigenvector_scaling = 60;
bξ_scaled = eigenvector_scaling*bξ
tξ_scaled = eigenvector_scaling*tξ
end

# ╔═╡ 98c9ae61-0a84-48c0-8773-0fe9bbd260ac
md"""
#### Compute an numerical approximation of the Taylor coefficients with Newton's method
"""

# ╔═╡ 5838bc82-6896-49a7-bc2c-b700fb3a3eb8
begin

N = 150 #Choose a Taylor order

initialguess = zeros(Taylor(N)^3) # The initial guess consists of the linear approximation of the manifold
for i = 1:3
    component(initialguess, i)[0] = bx[1]
    component(initialguess, i)[1] = bξ[1]
end

a, _ = newton(a -> (F(a, bx, bλ, bξ_scaled, par, N), DF(a, bλ, par, N)), initialguess)

end

# ╔═╡ 0bff34a8-2ad3-43d8-9f25-60a926f5b5c6
# md"""
# ##### Plot the image of the parameterization
# """

# ╔═╡ 3c23bf6e-fd79-4770-83c3-385832eddac6
# begin
# 	using OrdinaryDiffEq
# 	using GLMakie
# 	function lorenz(du, u, p, t)
# 	    du[1] = 10.0(u[2] - u[1])
# 	    du[2] = u[1] * (28.0 - u[3]) - u[2]
# 	    du[3] = u[1] * u[2] - (8 / 3) * u[3]
# 	end
# 	u0 = [0.0; 10.0; 10.0]
# 	tspan = (0.0, 50.0)
# 	prob = ODEProblem(lorenz, u0, tspan)
# 	sol = solve(prob, Tsit5())
# 	fig = Figure()
# 	ax = Axis3(fig[1,1]; aspect = :data, azimuth = 3π/5)
# 	lines!(ax, [Point3f(sol(t)) for t = LinRange(tspan..., 10_001)]; color = :royalblue1)
# 	lines!(ax, [Point3f(a(t)) for t = LinRange(-1, 1, 101)]; color = :salmon, linewidth = 3)
# 	meshscatter!(ax, [Point3f(bx), Point3f(-bx[1], -bx[2], bx[3]), Point3f(0,0,0)]; markersize = 1, color = :palegreen3)
# 	fig
# end

# ╔═╡ c4363fdf-6500-4830-bc42-9777d616ac34
md"""
### Choosing the approximate inverse $A$ for the Fréchet derivative $DF(\bar{a})$
"""

# ╔═╡ e72d43a1-b804-455a-bff0-8b1b715ddcee
md"""
The derivative of $F(a)=a-\mathcal{K}\phi(a)-c$ can now be written compactly as

$\begin{align}
DF(a) = I - \mathcal{K} D\phi(a)
\end{align}$

acting on an element $b=(b_1,b_2,b_3) \in \mathcal{X}$ as

$\begin{align}
DF(a)b = b  - \mathcal{K} D\phi(a) b.
\end{align}$
"""

# ╔═╡ 24e1a5af-569e-4ae8-a3d7-19b7962db48b
md"""
We use the projection to split $A$ into a finite part and a tail part:

$\begin{align}A := A^{\le N}\pi^{\le N} + I^{>N}\pi^{>N},\end{align}$

where $A^{\le N}$ is a linear map (to be chosen below) on $\pi^{\le N} \mathcal{X}$ and $I^{>N}$ is the identity on $\pi^{> N} \mathcal{X}$. Hence by construction one may also write $A$ as:

$\begin{align}A = \pi^{\le N}A^{\le N}\pi^{\le N} + \pi^{>N}I^{>N}\pi^{>N}.\end{align}$
"""

# ╔═╡ 922f5fa4-f9db-4eae-bd81-dc3dd525d212
md"""
We compute a numerical inverse of the Jacobian $DF^{\le N}(\bar{a})$ and we denote this inverse matrix by $A^{\le N}$.
Moreover, we compute the operator norm of the operator $A$.
"""

# ╔═╡ bb5ea822-f29b-499c-bebb-53208cf42eb1
begin
A_finite = interval.( inv(DF(a, bλ, par, N)) )
opnorm_A_finite = opnorm(A_finite, 1)
opnormA = max(opnorm_A_finite, 1)
end

# ╔═╡ 0d75d077-53c3-4f03-99d3-af796d973d66
md"""
## The Newton-Kantorovich proof
"""

# ╔═╡ 9433b6f1-3460-48c4-b587-9009823327c0
md"""
We are now in the Newton-Kantorovich context and thus need to derive the appropriate bounds.
"""

# ╔═╡ 4741e38b-0c29-4f02-a540-c543f7b3f234
md"""
The following code is used to compute the operator $\pi^{> N}$ acting on a single sequence or a list of sequences
"""

# ╔═╡ 656ae53a-1dbc-4abc-8837-234762db80c0
function tail(x, N)
    # get the tail of a Taylor series
    y = zeros(eltype(x), space(x))
    for n = N+1:order(x)
        y[n] = x[n]
    end
    return y
end

# ╔═╡ 9fe5120a-824a-4395-afe8-325a36ab9a07
md"""
### The bound $Y$
"""

# ╔═╡ eb2ae1fd-5062-436b-ac50-eb3b169bca42
md"""
Recall that that bound $Y$ satisfies $\| A F(\bar{a}) \|_{\mathcal{X}} \le Y$. Note that $A F(\bar{a}) = \pi^{\le 2N} A F(\bar{a}) + \stackrel{=0}{\overbrace{\pi^{> 2N} A F(\bar{a})}}$, since $(\bar{a}_1*\bar{a}_2)_n=(\bar{a}_1*\bar{a}_3)_n=0$ for $n>2N$. Hence,

$\begin{align}
A F(\bar{a}) &= \pi^{\le 2N} A F(\bar{a})  \\
& = \pi^{\le N} A F(\bar{a})+ (\pi^{\le 2N} - \pi^{\le N}) A F(\bar{a}) \\
& = \pi^{\le N} A^{\le N} F^{\le N}(\bar{a}) + (\pi^{\le 2 N}-\pi^{\le N}) (\bar{a} - \mathcal{K} \phi(\bar{a})-c) \\
& = \pi^{\le N} A^{\le N} F^{\le N}(\bar{a}) - (\pi^{\le 2 N}-\pi^{\le N}) \mathcal{K} \phi(\bar{a})
\end{align}$

The code just below computes with interval arithmetic the bound $Y$ satisfying

$\begin{align}
\| A F(\bar{a}) \|_{\mathcal{X}}
\le \| \pi^{\le N} A^{\le N} F^{\le N}(\bar{a}) \|_{\mathcal{X}} + \frac{1}{|\lambda| (N+1)} \| (\pi^{\le 2 N}-\pi^{\le N}) \phi(\bar{a}) \|_{\mathcal{X}} \le Y.
\end{align}$
"""

# ╔═╡ 3021e7c0-508a-4a95-9c16-098fa35b65a0
begin
full_F = F(interval.(a), tx, tλ, tξ_scaled, par, 2N)
Y = norm(A_finite * full_F, 1) + sum(i -> norm(tail(component(full_F, i), N), 1), 1:3) / (abs(tλ)*(N+1))
end

# ╔═╡ fc18e3ff-4a39-4aba-a2d2-94521035013e
md"""
### The bound $Z_1$
"""

# ╔═╡ 40a4ef18-8cf6-46bf-bf84-a7c89d898f51
md"""
We begin with three Exercises that will be useful in deriving the $Z_1$ bound.
"""

# ╔═╡ 454efc82-4918-4c60-8b9d-fb630302beec
Markdown.MD(Markdown.Admonition("note", "Exercise 1 for Z_1",[md""" Show that $\| \pi^{>N} \mathcal{K} \|_{B(\mathcal{X})} \le \frac{1}{|\lambda|(N+1)}$.
"""]))

# ╔═╡ c52ec039-8ef2-41fc-a9cd-34fcef45f735
Markdown.MD(Markdown.Admonition("note", "Exercise 2 for Z_1",[md""" Show that
$\begin{align} \| D\phi(\bar{a}) \|_{B(\mathcal{X})}
= \max_{i=1,2,3} \sum_{j=1}^3 \| \partial_i \phi_j(\bar{a})\|_1 \end{align}
= \max \left\{
\sigma + \|\rho-a_3\|_1 + \|a_2\|_1,
\sigma + 1 + \|a_1\|_1,
\|a_1\|_1 + \beta
\right\}.$
"""]))

# ╔═╡ 7229c43e-dcda-436c-86fc-1e626fe43daf
Foldable("""Hint""", md"""
Note that

$\begin{align}
D\phi(a) = \begin{pmatrix}
-\sigma I & \sigma I & 0 \\
\rho I - \mathcal{M}_{a_3} & -I & - \mathcal{M}_{a_1} \\
\mathcal{M}_{a_2} & \mathcal{M}_{a_1} & - \beta I
\end{pmatrix}
\end{align},$

where
$\mathcal{M}_{a_i} \in B(\ell^1)$ is the multiplication operator defined by $\mathcal{M}_{a_i} b = a_i*b$ with norm $\| \mathcal{M}_{a_i} \|_{B(\ell^1)} = \| a_i\|_1$.
""")

# ╔═╡ bc448153-fe09-4a05-b857-5b9fd51dcdae
Markdown.MD(Markdown.Admonition("note", "Exercise 3 for Z_1",[md""" Show that
$\pi^{\le N} D\phi(\bar{a})\pi^{> N}=0$.
"""]))

# ╔═╡ 1d3d35a7-c7a8-419e-b751-a144b86bbe1b
Foldable("""Hint""", md"""
Note that each none blocks of $D\phi(\bar{a})$ is lower triangular since of the Cauchy products.
""")

# ╔═╡ 002c7c0a-2b38-4a31-9e08-a424fa5679b8
md"""
Recall that the bound $Z_1$ satisfies
$\begin{align}
\| I - A DF(\bar{a}) \|_{B(\mathcal{X})} \le Z_1.
\end{align}$

We use a similar splitting as in the case of periodic orbits using Fourier series (see Part 3 of the tutorial), and observe that

$\begin{align}
\|\Gamma\|_{B(\mathcal{X})} = \max \Big( \|\Gamma \pi^{\le N}\|_{B(\mathcal{X})}, \|\Gamma \pi^{> N}\|_{B(\mathcal{X})} \Big).
\end{align}$

We compute separately two bounds $M_1,M_2$ satisfying $\|(I - A DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})} \le M_1$ and $\| (I - A DF(\bar{a})) \pi^{> N}\|_{B(\mathcal{X})} \le M_2$.

Recall that $A = A^{\le N}\pi^{\le N} + I^{>N}\pi^{>N}$.
"""

# ╔═╡ 9f022768-29a7-4ce9-bf9a-ceb04ea0adf2
md"""
###### Computing $M_1$
"""

# ╔═╡ 090c74fc-9bf7-47c0-af5f-887d6bf43ccd
md"""
Note first that

$\begin{align}
\|(I - A DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})} &\le
\|(I - A^{\le N}\pi^{\le N} DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})}
 + \| \pi^{>N} DF(\bar{a}) \pi^{\le N}\|_{B(\mathcal{X})},
\end{align}$

where $\pi^{>N} DF(\bar{a}) \pi^{\le N} = \pi^{>N} (I-\mathcal{K} D\phi(\bar{a})) \pi^{\le N}
= \pi^{>N} \pi^{\le N} - \pi^{>N} \mathcal{K} D\phi(\bar{a}) \pi^{\le N}
= - \pi^{>N} \mathcal{K} D\phi(\bar{a}) \pi^{\le N}$, and hence by Exercises 1 and 2:

$\begin{align}
\| \pi^{>N} DF(\bar{a}) \pi^{\le N}\|_{B(\mathcal{X})} & \le
\| \pi^{>N} \mathcal{K} \|_{B(\mathcal{X})} \| \pi^{>N} D\phi(\bar{a}) \|_{B(\mathcal{X})} \\
& \le \frac{1}{|\lambda|(N+1)} \| D\phi(\bar{a}) \|_{B(\mathcal{X})} \\
& \le \frac{1}{|\lambda|(N+1)} \max \left\{
\sigma + \|\rho-a_3\|_1 + \|a_2\|_1,
\sigma + 1 + \|a_1\|_1,
\|a_1\|_1 + \beta
\right\}
\end{align}$

Using interval arithmetic, one can compute an upper bound $M_1$ such that

$\begin{align}
\|(I - A DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})} &\le
\|(I - A^{\le N}\pi^{\le N} DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})} +
\frac{1}{|\lambda|(N+1)} \max \left\{
\sigma + \|\rho-a_3\|_1 + \|a_2\|_1,
\sigma + 1 + \|a_1\|_1,
\|a_1\|_1 + \beta
\right\}\le M_1
\end{align}$
"""

# ╔═╡ 9e8d66d9-ad6e-42b5-a24f-6d0a9060f688
md"""
###### Computing $M_2$
"""

# ╔═╡ 6297e819-a105-4451-bbb7-c2241ff77f51
md"""
Now, for $M_2$:

$\begin{align}
\|(I - A DF(\bar{a})) \pi^{> N}\|_{B(\mathcal{X})} & =
\|(I - (A^{\le N}\pi^{\le N} + I^{>N}\pi^{>N}) DF(\bar{a})) \pi^{> N}\|_{B(\mathcal{X})} \\
&\le
\| A^{\le N}\stackrel{=0}{\overbrace{\pi^{\le N} DF(\bar{a})\pi^{> N}}}\|_{B(\mathcal{X})}
 + \| \pi^{> N}-\pi^{>N} DF(\bar{a}) \pi^{> N}\|_{B(\mathcal{X})},
\end{align}$

since $\pi^{\le N} DF(\bar{a})\pi^{> N} = \pi^{\le N} (I- \mathcal{K} D\phi(\bar{a}))\pi^{> N}
=\pi^{\le N}\pi^{> N} - \pi^{\le N} D\phi(\bar{a}))\pi^{> N} = - \pi^{\le N} D\phi(\bar{a}))\pi^{> N}=0$ by Exercise 3.

Now, using a similar idea as for $M_1$, one can show that

$\begin{align}
\| \pi^{> N}-\pi^{>N} DF(\bar{a}) \pi^{> N}\|_{B(\mathcal{X})} \le \frac{1}{|\lambda|(N+1)} \max \left\{
\sigma + \|\rho-a_3\|_1 + \|a_2\|_1,
\sigma + 1 + \|a_1\|_1,
\|a_1\|_1 + \beta
\right\}\end{align}$.
"""

# ╔═╡ fb8527a5-4fa9-40d3-b0af-2cce8cacac3d
md"""
###### Computing $Z_1$
"""

# ╔═╡ 7e60324b-f65e-4d08-b642-13542b3637f3
md"""
Since clearly $M_2 \le M_1$, then we let $Z_1=M_1$. Then, using interval arithmetic, one can compute an upper bound $Z_1$ such that

$\begin{align}
\| I - A DF(\bar{a}) \|_{B(\mathcal{X})}
 &\le
\|(I - A^{\le N}\pi^{\le N} DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})} +
\frac{1}{|\lambda|(N+1)} \max \left\{
\sigma + \|\rho-a_3\|_1 + \|a_2\|_1,
\sigma + 1 + \|a_1\|_1,
\|a_1\|_1 + \beta
\right\}\le Z_1
\end{align}$
which we do below.
"""

# ╔═╡ c5d63b4c-7823-4122-a81f-9e3a18a634ff
begin
a₁, a₂, a₃ = eachcomponent(interval.(a))

opnorm_Dϕ = max(
    tσ + norm(ρ - a₃, 1) + norm(a₂, 1),
    tσ + 1               + norm(a₁, 1),
         norm(a₁, 1)     + β)

Z₁ = opnorm(I - A_finite * DF(interval.(a), tλ, tpar, N), 1) + opnorm_Dϕ / (abs(tλ)*(N+1))
end

# ╔═╡ c0ce2c85-eaae-4b08-8ef7-2a9297b10e05
md"""
### The bound $Z_2$
"""

# ╔═╡ 6f346c67-879a-4427-afd4-6acdeb17af21
md"""
This concerns the bound on the operator $A[DF(a)-DF(\bar{a})]$. Since this is a Lipschitz bound on the first derivative, one approach is to bound the second derivative. However, here we just do a direct computation to see that
$A[DF(a)-DF(\bar{a})]b = - A \mathcal{K} [D\phi(a)-D\phi(\bar{a})]b$.

Note that

$\begin{align}
(D\phi(a) - D\phi(\bar{a}))b = \begin{pmatrix}
0 & 0 & 0 \\
0 - (\mathcal{M}_{a_3}-\mathcal{M}_{\bar{a}_3}) & 0 & - (\mathcal{M}_{a_1}-\mathcal{M}_{\bar{a}_1}) \\
(\mathcal{M}_{a_2}-\mathcal{M}_{\bar{a}_2}) & (\mathcal{M}_{a_1}-\mathcal{M}_{\bar{a}_1}) & 0
\end{pmatrix} \begin{pmatrix} b_1 \\ b_2 \\ b_3 \end{pmatrix}
= \begin{pmatrix} 0 \\ -(a_3-\bar{a}_3)*b_1 -(a_1-\bar{a}_1)*b_3 \\ (a_2-\bar{a}_2)*b_1 + (a_1-\bar{a}_1)*b_2\end{pmatrix}
\end{align},$

hence by the Banach algebra property and the definition of the operator norm we arrive at

$\begin{align}
\|A[DF(a)-DF(\bar{a})]b\|_{\mathcal{X}} &\le \|A\|_{B(\mathcal{X})} \| \mathcal{K} \|_{B(\mathcal{X})}
(
\| a_3-\bar{a}_3 \|_1 \|b_1\|_1 + \| a_1-\bar{a}_1\|_1 \|b_3\|_1 + \|a_2-\bar{a}_2\|_1 \|b_1\|_1 + \|a_1-\bar{a}_1\|_1 \|b_2\|_1
)
\\
& \le \|A\|_{B(\mathcal{X})} \| \mathcal{K} \|_{B(\mathcal{X})} \|a-\bar{a}\|_{\mathcal{X}} \|b\|_{\mathcal{X}}\\
& = \frac{1}{2|\lambda|} \|A\|_{B(\mathcal{X})} \|a-\bar{a}\|_{\mathcal{X}} \|b\|_{\mathcal{X}}.
\end{align}$

Computing $Z_2$ boils down to compute (or bound) $\frac{1}{2|\lambda|} \|A\|_{B(\mathcal{X})}$.
"""

# ╔═╡ 6688fe5d-0019-485d-9059-bb5bd2be7d6e
Z₂ = (2/abs(tλ)) * opnormA

# ╔═╡ 5d0b6b50-0170-4f8b-ad5a-0dda212a00e1
md"""
## Finishing the CAP
"""

# ╔═╡ a9f597fe-4058-47a8-8a39-8d6b36594a96
md"""
We now evaluate the radii polynomial to finish the Newton-Kantorovich proof.
"""

# ╔═╡ 1082f423-448c-4d74-8baa-76b2e6508084
begin
r_star = Inf
r = interval_of_existence(Y, Z₁, Z₂, r_star)
end

# ╔═╡ Cell order:
# ╟─97c19551-3248-47f6-a1cb-bc4fb732e6df
# ╠═97df61b8-00b4-4c51-bbc0-b09df0e4cbe1
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─8f59e349-74a6-4a6f-9a47-6ca26d187aeb
# ╟─4412bd33-8b3f-41de-9e1a-a8f1187c793f
# ╟─d9604e75-0e9e-44b6-bf2e-885f01a96bdc
# ╟─62719aa9-bfda-4913-8b74-f3694bfc33c0
# ╟─fc381b29-6671-4364-a9df-b35db07af45b
# ╟─194ee83b-7442-42f6-862a-2f5bd1c8b710
# ╟─9829987c-2bf5-4b33-8ae1-6b80d43ed978
# ╟─1f8fba92-d023-4444-b4e3-18c751f678c5
# ╟─11aa1d55-1a2c-48cf-a7f5-30ff16e2e92a
# ╟─668058ef-e66b-4fd4-ad67-2a8389589e86
# ╟─0562b884-0653-4d5b-8a66-e35467e6693e
# ╟─a47a4e2a-371e-4244-8574-493f00d456df
# ╟─beca7238-ca9b-409a-8d5a-240ca720f58a
# ╠═c8c8856f-fe43-4b77-aad9-57723b6dcf9c
# ╠═5e6482db-9b76-4bf0-bc58-0a6f779c4556
# ╟─8c0ac2a2-b3db-455a-8a6d-3de6728b61d3
# ╟─7cfbe32b-d677-4611-97d9-ed775c60e48e
# ╟─13f3b4e4-0409-49a2-a544-0c34ff9289db
# ╟─40156798-59d5-4afe-a054-29c753a51e46
# ╟─42ad70a3-114e-4cd3-8dd1-2429ac80bfc2
# ╠═cdbbb2e1-56bc-4beb-93d3-2db1a2aa71d5
# ╟─aaf278f6-b8ee-4094-9173-ffde2cb0b461
# ╠═d77a192e-4be0-4987-8023-8d91b22d8d9f
# ╟─7410e377-715c-4464-a83d-d7cd48c4d1ae
# ╠═fbeac0c7-8fcb-4e4b-9662-9f1d48b23fb3
# ╟─ae88c826-f486-4786-9447-affcc466d1b5
# ╟─76aa085c-d7ab-4e2f-9648-89b3eebe11cd
# ╟─585b99f6-fad0-44b1-b1ca-f7aa6f8f52f8
# ╟─400fe5e0-0dfd-4d82-baf9-b66d931133a8
# ╠═e0f75d0e-555d-4df7-8077-9871a6d7b197
# ╟─deb48c67-19f1-481e-b2de-1765e449e6ae
# ╠═c106822b-2c0c-4759-84fb-9ddba91fe295
# ╟─f75c26b4-8e6f-4520-b1a8-f30c20b9a194
# ╠═05a4dd64-b5c9-4014-a18c-f26938e02ebd
# ╟─3cb601bf-43ea-430c-950d-402b0929e515
# ╠═2a76f46f-c203-4028-80ba-e7bc1f481076
# ╟─8a802d5f-e2f5-4e03-b0d9-b7cd0c081d2c
# ╠═5ef06ffe-d2f3-4cd2-946a-df8fd78f863f
# ╟─e4a015b0-e7e6-4841-b97f-4a845f462196
# ╠═c3d0d2ae-5711-451b-946d-6cdca682fe06
# ╟─4ff1b050-3036-4e0b-b9db-b07b4d17327c
# ╠═ec3cd93f-877d-419e-a608-cef748306927
# ╟─6e715b86-4604-49c3-8df1-929921e38f9a
# ╠═ffe5218b-88ca-42e0-8a0d-074ab2e6d898
# ╟─98c9ae61-0a84-48c0-8773-0fe9bbd260ac
# ╠═5838bc82-6896-49a7-bc2c-b700fb3a3eb8
# ╟─0bff34a8-2ad3-43d8-9f25-60a926f5b5c6
# ╟─3c23bf6e-fd79-4770-83c3-385832eddac6
# ╟─c4363fdf-6500-4830-bc42-9777d616ac34
# ╟─e72d43a1-b804-455a-bff0-8b1b715ddcee
# ╟─24e1a5af-569e-4ae8-a3d7-19b7962db48b
# ╟─922f5fa4-f9db-4eae-bd81-dc3dd525d212
# ╠═bb5ea822-f29b-499c-bebb-53208cf42eb1
# ╟─0d75d077-53c3-4f03-99d3-af796d973d66
# ╟─9433b6f1-3460-48c4-b587-9009823327c0
# ╟─4741e38b-0c29-4f02-a540-c543f7b3f234
# ╠═656ae53a-1dbc-4abc-8837-234762db80c0
# ╟─9fe5120a-824a-4395-afe8-325a36ab9a07
# ╟─eb2ae1fd-5062-436b-ac50-eb3b169bca42
# ╠═3021e7c0-508a-4a95-9c16-098fa35b65a0
# ╟─fc18e3ff-4a39-4aba-a2d2-94521035013e
# ╟─40a4ef18-8cf6-46bf-bf84-a7c89d898f51
# ╟─454efc82-4918-4c60-8b9d-fb630302beec
# ╟─c52ec039-8ef2-41fc-a9cd-34fcef45f735
# ╟─7229c43e-dcda-436c-86fc-1e626fe43daf
# ╟─bc448153-fe09-4a05-b857-5b9fd51dcdae
# ╟─1d3d35a7-c7a8-419e-b751-a144b86bbe1b
# ╟─002c7c0a-2b38-4a31-9e08-a424fa5679b8
# ╟─9f022768-29a7-4ce9-bf9a-ceb04ea0adf2
# ╟─090c74fc-9bf7-47c0-af5f-887d6bf43ccd
# ╟─9e8d66d9-ad6e-42b5-a24f-6d0a9060f688
# ╟─6297e819-a105-4451-bbb7-c2241ff77f51
# ╟─fb8527a5-4fa9-40d3-b0af-2cce8cacac3d
# ╟─7e60324b-f65e-4d08-b642-13542b3637f3
# ╠═c5d63b4c-7823-4122-a81f-9e3a18a634ff
# ╟─c0ce2c85-eaae-4b08-8ef7-2a9297b10e05
# ╟─6f346c67-879a-4427-afd4-6acdeb17af21
# ╠═6688fe5d-0019-485d-9059-bb5bd2be7d6e
# ╟─5d0b6b50-0170-4f8b-ad5a-0dda212a00e1
# ╟─a9f597fe-4058-47a8-8a39-8d6b36594a96
# ╠═1082f423-448c-4d74-8baa-76b2e6508084
