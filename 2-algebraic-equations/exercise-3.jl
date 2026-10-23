### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Radius of convergence"

using Markdown
using InteractiveUtils

# ╔═╡ 7756e7d4-a00c-4e9d-9bf9-866ab34907b9
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ b186bc19-ff2d-4914-914a-b423dcfea0f1
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

# ╔═╡ 1a0e2f9c-dd27-4141-97fc-a4ebddb61cb8
md"""
For $\nu \ge 1$, consider the weighted sequence space

```math
\ell^1_{\nu, \mathbb{N}} \overset{\text{def}}{=} \left\{ a \in \mathbb{R}^\mathbb{N} \, : \, \| a \|_\nu \overset{\text{def}}{=} \sum_{n \ge 0} |a_n| \nu^n < \infty \right\}.
```
"""

# ╔═╡ 15be0e8a-408b-4db2-af7e-15261f54238e
md"""
**1.** Show that $\ell^1_{\nu, \mathbb{N}}$ equipped with the Cauchy product is a Banach algebra.
"""

# ╔═╡ 4ea373ca-e44b-49cf-9e6c-74a5043dee79
md"""
**2.** Show that $\|M\|_{\mathscr{B}(\ell^1_{\nu, \mathbb{N}})} = \max_{l \ge 0} \nu^{-l} \sum_{n \ge 0} |M_{n,l}| \nu^n$ for all $M \in \mathscr{B}(\ell^1_{\nu, \mathbb{N}})$.
In particular, for all $x \in \ell^1_{\nu, \mathbb{N}}$, show that $\| \mathcal{M}_x \|_{\mathscr{B}(\ell^1_{\nu, \mathbb{N}})} = \| x \|_\nu$, where $\mathcal{M}_x$ is the multiplication operator associated with $x$.
"""

# ╔═╡ c3cf2193-105f-46da-8bbc-9a95082a805d
md"""
**3.** Find the inverse of $v(t) = 2 + t$ for all $t \in [-\nu, \nu]$ as a Taylor series.
How large can you take $\nu$?
"""

# ╔═╡ 0746bcdf-dd75-4726-92dc-696ec0a510b1
md"""
**4.** Find the inverse of $v(t) = e^t$ for all $t \in [-\nu, \nu]$ as a Taylor series.
How large can you take $\nu$?
"""

# ╔═╡ Cell order:
# ╟─7756e7d4-a00c-4e9d-9bf9-866ab34907b9
# ╠═b186bc19-ff2d-4914-914a-b423dcfea0f1
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─1a0e2f9c-dd27-4141-97fc-a4ebddb61cb8
# ╟─15be0e8a-408b-4db2-af7e-15261f54238e
# ╟─4ea373ca-e44b-49cf-9e6c-74a5043dee79
# ╟─c3cf2193-105f-46da-8bbc-9a95082a805d
# ╟─0746bcdf-dd75-4726-92dc-696ec0a510b1
