### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Rigorous control of the entire spectrum"

using Markdown
using InteractiveUtils

# ╔═╡ 7287283f-16f1-4336-bbf2-b8b06b5177ad
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ bf81d459-ab70-4b6b-abcb-5c1a4fb67f11
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
    using PlutoTeachingTools
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

# ╔═╡ 45e10c91-f161-43a5-9c9e-5b4dca6a8e53
md"""
Exercise 3 can be used to enclose eigenvalues of a given matrix one by one.
We now present an alternate strategy to enclose the entire spectrum at once (but not the corresponding eigenvectors), which can sometimes also be adapted in infinite dimension.
This strategy relies on the [Gershgorin circle theorem](https://en.wikipedia.org/wiki/Gershgorin_circle_theorem), which we recall below.
"""

# ╔═╡ 916f5d37-2bde-4a35-923a-2ff9e606371d
md"""
!!! theorem "Theorem (Gershgorin)"
	For any matrix $M = \left(M_{i,j}\right)_{1 \le i, j \le d}$, its spectrum $\sigma(M)$ satisfies

	```math
	\begin{align}
	\sigma(M) \subset \bigcup_{i=1}^d D\left( M_{i,i},\, \sum_{j \ne i} | M_{i,j} | \right),
	\end{align}
	```

	where $D(z,r)$ denotes the closed disk of center $z$ and radius $r$ in the complex plane.
	Moreover, if $\mathscr{J} \subset \{1, \dots, d\}$ is such that $\bigcup_{i \in \mathscr{J}} D\left(M_{i,i},\, \sum_{j \ne i} | M_{i,j} |\right)$ is disjoint from $\bigcup_{i \notin \mathscr{J}} D\left(M_{i,i},\, \sum_{j \ne i} | M_{i,j} |\right)$, then $\bigcup_{i \in \mathscr{J}} D\left(M_{i,i},\, \sum_{j \ne i} | M_{i,j} |\right)$ contains exactly $| \mathscr{J} |$ eigenvalues.
"""

# ╔═╡ d01c5817-c4b6-4502-81f6-51ae5715117b
md"""
**1.** Using the Gershgorin circle theorem, get as tight as possible rigorous enclosures of all eigenvalues of $W_{3}$ (defined in Exercise 3).
"""

# ╔═╡ 02939955-c9aa-4152-8e08-f31e1e0c0e9c
Foldable("Hint",
md"""
You may first compute numerically a matrix $P$ of approximate eigenvectors of $W_3$, then rigorously compute $\tilde W_{3} = P^{-1}W_3 P$, and finally apply the Gershgoring circle theorem to $\tilde W_{3}$.
"""
)

# ╔═╡ bc052916-db0e-43a2-bbb0-76b917d6f638
md"""
**2.** Prove that $W_{1000}$ has exactly one eigenvalue with negative real part.
"""

# ╔═╡ 3686e5b2-daac-470d-b9e3-4bd5706e5894
Foldable("Hint",
md"""
You do not need to numerically diagonalize all of $W_{1000}$: for most rows, the corresponding Gershgorin disk already lies in the left half of the complex plane.
""")

# ╔═╡ Cell order:
# ╟─7287283f-16f1-4336-bbf2-b8b06b5177ad
# ╠═bf81d459-ab70-4b6b-abcb-5c1a4fb67f11
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─45e10c91-f161-43a5-9c9e-5b4dca6a8e53
# ╟─916f5d37-2bde-4a35-923a-2ff9e606371d
# ╟─d01c5817-c4b6-4502-81f6-51ae5715117b
# ╟─02939955-c9aa-4152-8e08-f31e1e0c0e9c
# ╟─bc052916-db0e-43a2-bbb0-76b917d6f638
# ╟─3686e5b2-daac-470d-b9e3-4bd5706e5894
