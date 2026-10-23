### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Rigorous inverse of a matrix"

using Markdown
using InteractiveUtils

# ╔═╡ 41c7c8f6-8876-4406-b326-9b4639c22ab3
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ 0af0057a-893d-4f5c-8157-61a55cc52a5d
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

# ╔═╡ 81ed6d83-622f-46ac-b0c3-0fae0c8ef378
md"""
In this exercise, we consider a matrix $M$, and a numerically computed approximate inverse $\bar{X}$ of $M$.
Our goal will be to guarantee a posteriori that $M$ is indeed invertible, and to provide a computable error bound between $M^{-1}$ and $\bar{X}$, for any submultiplicative matrix norm (i.e. $\| A B \| \le \| A \| \| B \|$).
"""

# ╔═╡ 6ad9d1f6-b2b3-49f7-b117-f2db2b7228fd
md"""
**1.** Denoting $\delta = \Vert I - M \bar{X} \|$, and assuming $\delta < 1$, show that $M$ is invertible and

```math
\| M^{-1} - \bar{X} \| \le \frac{\delta}{1-\delta} \| \bar{X} \|.
```
"""

# ╔═╡ 05ad3f54-e92c-4ab5-b276-aa7763ba36b3
Foldable("Hint",
md"""
Write $\left(M\bar{X}\right)^{-1}$ as $\left(I + M\bar{X} - I\right)^{-1}$ and mutliply to the left by $\bar{X}$ in order to get a power series expansion of $M^{-1}$.
"""
)

# ╔═╡ 82a5527f-0661-4a32-b758-5708bb184968
md"""
**2.** Try to obtain a similar estimate using an appropriate zero-finding problem and the Newton-Kantorovich approach.
"""

# ╔═╡ e34560b8-93af-4f5b-8e4f-f5c48ef29d3c
Foldable("Hint",
md"""
You may consider $F(X) = MX-I$ (or $F(X) = XM-I$).
""")

# ╔═╡ Cell order:
# ╟─41c7c8f6-8876-4406-b326-9b4639c22ab3
# ╠═0af0057a-893d-4f5c-8157-61a55cc52a5d
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─81ed6d83-622f-46ac-b0c3-0fae0c8ef378
# ╟─6ad9d1f6-b2b3-49f7-b117-f2db2b7228fd
# ╟─05ad3f54-e92c-4ab5-b276-aa7763ba36b3
# ╟─82a5527f-0661-4a32-b758-5708bb184968
# ╟─e34560b8-93af-4f5b-8e4f-f5c48ef29d3c
