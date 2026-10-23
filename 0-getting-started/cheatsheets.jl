### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Cheatsheets"

using Markdown
using InteractiveUtils

# ╔═╡ 5a4bebd1-2c7d-4c5c-a3f2-0976e3fd44ee
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ b2fa25d8-9577-47dd-acf5-8cdcc052f36d
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))   # the course's shared environment
    Pkg.instantiate()
    using IntervalArithmetic, RadiiPolynomial, Plots
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

# ╔═╡ 334ef43a-796d-471a-a82d-d608b92cdc14
md"### Cheatsheets"

# ╔═╡ 8407f06b-2392-49de-a734-b0339f7376fa
md"""
- [Fast track to Julia](https://juliadocs.github.io/Julia-Cheat-Sheet/)
- [MATLAB-Julia-Python comparative cheatsheet](https://cheatsheets.quantecon.org/)
- [Plots.jl cheatsheet](https://github.com/sswatson/cheatsheets/blob/master/plotsjl-cheatsheet.pdf)
"""

# ╔═╡ 55933fb6-6db8-4fad-9753-c317916660e4
md"### IntervalArithmetic"

# ╔═╡ 24ba79d7-b92c-4f42-a940-24907321104e
interval([1.0, 2.0])

# ╔═╡ 3744963a-1f6b-4dbc-b796-83dad65e87f5
interval([1.0 2.0 ; 3.0 4.0])

# ╔═╡ 1672e4dd-5f52-4be5-8d5f-2b7a77a6af99
md"### RadiiPolynomial"

# ╔═╡ 86b71145-bf35-4310-9914-78925141a5d1
x = Sequence(Taylor(1), [1.0, 1.0]) # 1 + t

# ╔═╡ 960765e6-2c84-4e91-ace1-191126fc6589
x * x # (1 + t)^2 = 1 + 2t + t^2

# ╔═╡ f2acf06e-4be5-4a30-b036-e554c0d2ffa4
Derivative(1)

# ╔═╡ 822cc3b2-dcec-46c8-a5c6-726de36f90e7
Derivative(1) * x

# ╔═╡ a9bf21bb-1a92-4942-b76a-aef1ffc98887
Multiplication(x)

# ╔═╡ 0cdf0577-d1b3-4c42-8090-6b3bc198b044
project(Multiplication(x), Taylor(5), Taylor(5))

# ╔═╡ 1547e260-583c-41ca-bbbe-29b79707708f
y = Sequence(Fourier(1, 1.0), [im*0.5, 0.0, -im*0.5]) # sin(t)

# ╔═╡ bf2d955d-9d9c-4ad3-85d5-d36696699b9b
plot(LinRange(-π, π, 101), t -> real(y(t)))

# ╔═╡ 6c89b8af-6c6c-46da-a441-52f75930b377
interval(x)

# ╔═╡ cdb5640a-2778-4285-a0cc-6c669bed839f
norm(x, 1)

# ╔═╡ 2583c7b4-f78d-47c8-a925-99cbe8ed6103
ν = interval(1.1)

# ╔═╡ 57a51438-53d3-4669-b410-c6338165043c
norm(interval(x), Ell1(GeometricWeight(ν)))

# ╔═╡ Cell order:
# ╟─5a4bebd1-2c7d-4c5c-a3f2-0976e3fd44ee
# ╠═b2fa25d8-9577-47dd-acf5-8cdcc052f36d
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─334ef43a-796d-471a-a82d-d608b92cdc14
# ╟─8407f06b-2392-49de-a734-b0339f7376fa
# ╟─55933fb6-6db8-4fad-9753-c317916660e4
# ╠═24ba79d7-b92c-4f42-a940-24907321104e
# ╠═3744963a-1f6b-4dbc-b796-83dad65e87f5
# ╟─1672e4dd-5f52-4be5-8d5f-2b7a77a6af99
# ╠═86b71145-bf35-4310-9914-78925141a5d1
# ╠═960765e6-2c84-4e91-ace1-191126fc6589
# ╠═f2acf06e-4be5-4a30-b036-e554c0d2ffa4
# ╠═822cc3b2-dcec-46c8-a5c6-726de36f90e7
# ╠═a9bf21bb-1a92-4942-b76a-aef1ffc98887
# ╠═0cdf0577-d1b3-4c42-8090-6b3bc198b044
# ╠═1547e260-583c-41ca-bbbe-29b79707708f
# ╠═bf2d955d-9d9c-4ad3-85d5-d36696699b9b
# ╠═6c89b8af-6c6c-46da-a441-52f75930b377
# ╠═cdb5640a-2778-4285-a0cc-6c669bed839f
# ╠═2583c7b4-f78d-47c8-a925-99cbe8ed6103
# ╠═57a51438-53d3-4669-b410-c6338165043c
