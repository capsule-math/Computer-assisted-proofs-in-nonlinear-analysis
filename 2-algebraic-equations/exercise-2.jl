### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Inverse function for infinite sequences"

using Markdown
using InteractiveUtils

# ╔═╡ 4a72be52-115f-4fcb-be90-486f627eeca4
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ d50a2f5f-2e01-48d1-9fbe-360ee62320a9
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

# ╔═╡ f88bd126-8272-4c63-ac31-ce6ba72700e5
md"""
**1.** Let $v(t) = \sum_{n \ge 0} y_n t^n$ such that $\Pi_N y$ is known, while its tail satisfies the estimate $\| (I - \Pi_N) y \|_1 \le \epsilon$.
For the zero-finding problem $F(x) \overset{\text{def}}{=} x * y - 1$, derive explicit formulas for the bounds $Y, Z_1, Z_2$ in terms of $\Pi_N y$ and $\epsilon$.
"""

# ╔═╡ 800c404b-d415-4a93-a5c4-e2296068a534
md"""
**2.** Find the inverse of $v(t) = e^t$ as a Taylor series for all $t \in [-1, 1]$.
"""

# ╔═╡ d874d58b-3124-427e-b45a-6926173622fc
Foldable("Hint",
md"""
Recall that $v(t) = e^t = \sum_{n \ge 0} y_n t^n$ where $y_n = 1/n!$.
Using the Taylor remainder, obtain a bound for $\| (I - \Pi_N) y \|_1$.
""")

# ╔═╡ Cell order:
# ╟─4a72be52-115f-4fcb-be90-486f627eeca4
# ╠═d50a2f5f-2e01-48d1-9fbe-360ee62320a9
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─f88bd126-8272-4c63-ac31-ce6ba72700e5
# ╟─800c404b-d415-4a93-a5c4-e2296068a534
# ╟─d874d58b-3124-427e-b45a-6926173622fc
