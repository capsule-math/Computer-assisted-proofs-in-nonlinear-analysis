### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Sequence space properties"

using Markdown
using InteractiveUtils

# ╔═╡ 0b67c7ba-b589-4673-9e0d-4ba2b732eebb
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ 05452044-9b3f-4202-97e8-180c99020184
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

# ╔═╡ f88bd126-8272-4c63-ac31-ce6ba72700e5
md"""
**1.** Show that the sequence space

```math
\ell^1_\mathbb{N} \overset{\text{def}}{=} \left\{ x \in \mathbb{R}^\mathbb{N} \, : \, \| x \|_1 \overset{\text{def}}{=} \sum_{n \ge 0} |x_n| < \infty \right\}
```

is a Banach space.
"""

# ╔═╡ 800c404b-d415-4a93-a5c4-e2296068a534
md"""
**2.** Show that $\ell^1_\mathbb{N}$, equipped with the Cauchy product given by

```math
(x * y)_n \overset{\text{def}}{=} \sum_{l=0}^n x_{n-l} y_l, \qquad n \ge 0,
```

satisfies the Banach algebra property $\| x * y \|_1 \le \| x \|_1 \| y \|_1$.
"""

# ╔═╡ fb9aaf0f-5266-4b15-bf12-1bf537345544
md"""
**3.** Let $x \in \ell^1_\mathbb{N}$ and $\mathcal{M}_x : \ell^1_\mathbb{N} \to \ell^1_\mathbb{N}$ the multiplication operator associated with $x$, i.e. $\mathcal{M}_x h \overset{\text{def}}{=} x * h$ for all $h \in \ell^1_\mathbb{N}$.
Show that $\| \mathcal{M}_x \|_{\mathscr{B}(\ell^1_\mathbb{N})} = \| x \|_1$.
"""

# ╔═╡ Cell order:
# ╟─0b67c7ba-b589-4673-9e0d-4ba2b732eebb
# ╠═05452044-9b3f-4202-97e8-180c99020184
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─f88bd126-8272-4c63-ac31-ce6ba72700e5
# ╟─800c404b-d415-4a93-a5c4-e2296068a534
# ╟─fb9aaf0f-5266-4b15-bf12-1bf537345544
