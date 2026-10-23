### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "2D manifold"

using Markdown
using InteractiveUtils

# ╔═╡ b6c2d706-a283-47a7-827d-a65cae7b8e88
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ 35d4db2e-4cd3-4bcb-89a5-1c9ec6b18425
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
end

# ╔═╡ a3982eb8-7843-11ef-0e0b-6be567b62442
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

# ╔═╡ 7e454145-e9b8-4348-85d7-56b79edea6ea
md"""
Compute the 2D stable manifold of the origin for Lorenz.
"""

# ╔═╡ Cell order:
# ╟─b6c2d706-a283-47a7-827d-a65cae7b8e88
# ╠═35d4db2e-4cd3-4bcb-89a5-1c9ec6b18425
# ╟─a3982eb8-7843-11ef-0e0b-6be567b62442
# ╟─7e454145-e9b8-4348-85d7-56b79edea6ea
