### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Rigorous computation of an eigenpair"

using Markdown
using InteractiveUtils

# ╔═╡ c8978e2a-0308-4549-b00c-ae24ef25c208
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ 7117324e-6419-4fb1-b44f-1bccf47657b1
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
    using PlutoTeachingTools, LinearAlgebra
end

# ╔═╡ 7fc40507-eda3-474d-a454-04e9173a7adb
html"""
<style>
main {
    max-width: 1000px;
    margin-left: auto;
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ c0a3bcb6-33b5-40a9-9696-7e37a2c9c432
md"""
**1.** Consider a matrix $M$, and an approximate eigenpair $(\bar{\lambda}, \bar{v})$ of $M$.
The goal of this exercise is to validate $(\bar{\lambda}, \bar{v})$, i.e. to prove that there exists an exact eigenpair $(\tilde{\lambda}, \tilde{v})$ nearby.
Assuming that the corresponding exact eigenvalue $\tilde{\lambda}$ is simple, define a suitable $F = 0$ problem, and derive the bounds needed to apply the Newton-Kantorovich theorem in that context.
"""

# ╔═╡ 7748e568-afc9-43cc-b2bd-5a231d86f455
Foldable("Hint",
md"""
The "natural" zero-finding problem is $G(\lambda, v) \overset{\text{def}}{=} (M - \lambda I)v$, but it has one too many unknowns.
This is consistent with the fact that zeros of $G$ are not isolated (one can always rescale the eigenvector).
Therefore, a suitable zero-finding problem needs to incorporate a normalization condition, for instance:

```math
F(\lambda, v) \overset{\text{def}}{=}
\begin{pmatrix}
\langle v, \bar{v} \rangle - 1 \\
(M - \lambda I)v
\end{pmatrix}.
```
"""
)

# ╔═╡ cab728f1-9ff5-4bdd-8101-5c39718c4d53
md"""
**2.** For any positive integer $N$, the Wilkinson matrix $W_{N}$ is the following $(2N+1)\times(2N+1)$ tridiagonal matrix:

```math
W_{N} =
\begin{pmatrix}
N & 1 & & & & & \\
1 & N-1 & 1 & & & & \\
 & 1 & \ddots & \ddots & & & \\
 & & \ddots & 0 & \ddots & & \\
 & & & \ddots & \ddots & 1 & \\
 & & & & 1 & N-1 & 1 \\
 & & & & & 1 & N
\end{pmatrix}
```

We provide below approximate eigenvalues and eigenvectors of $W_3$.
Rigorously enclose all eigenpairs of $W_3$.
"""

# ╔═╡ d61514c3-3b0e-4658-8b31-de9f9514a9c3
function wilkinson(N)
	M = zeros(2N+1, 2N+1)
	for i = 1:2N+1
		M[i,i] = abs(N - i + 1)
		if i+1 ≤ 2N+1
			M[i,i+1] = 1
		end
		if i-1 ≥ 1
			M[i,i-1] = 1
		end
	end
	return M
end

# ╔═╡ 1bba510b-be86-44b0-a3c9-419b3b6ada37
N = 3

# ╔═╡ fe0054f0-4fd5-489f-9fcb-3af086876699
W = wilkinson(N)

# ╔═╡ 3509fe96-4a83-461f-8fed-23343d74dc8c
eigenvalues, eigenvectors = eigen(W)

# ╔═╡ Cell order:
# ╟─c8978e2a-0308-4549-b00c-ae24ef25c208
# ╠═7117324e-6419-4fb1-b44f-1bccf47657b1
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─c0a3bcb6-33b5-40a9-9696-7e37a2c9c432
# ╟─7748e568-afc9-43cc-b2bd-5a231d86f455
# ╟─cab728f1-9ff5-4bdd-8101-5c39718c4d53
# ╠═d61514c3-3b0e-4658-8b31-de9f9514a9c3
# ╠═1bba510b-be86-44b0-a3c9-419b3b6ada37
# ╠═fe0054f0-4fd5-489f-9fcb-3af086876699
# ╠═3509fe96-4a83-461f-8fed-23343d74dc8c
