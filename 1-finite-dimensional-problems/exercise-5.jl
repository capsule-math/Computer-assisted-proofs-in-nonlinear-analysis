### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Toward rigorous continuation"

using Markdown
using InteractiveUtils

# ╔═╡ 45f01879-f2e0-4732-8792-a3c7223b310b
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ 1b522e51-248f-4b20-9fb3-d6d9f67bc69a
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
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

# ╔═╡ afef60df-375e-439b-a973-ee36439731f2
md"""
Prove that there exists a period-3 orbit to $x_{n+1} = \mu x_n (1-x_n)$ for all $\mu$ in an interval containing $3.9$.
Make that interval as large as you can.
"""

# ╔═╡ Cell order:
# ╟─45f01879-f2e0-4732-8792-a3c7223b310b
# ╠═1b522e51-248f-4b20-9fb3-d6d9f67bc69a
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─afef60df-375e-439b-a973-ee36439731f2
