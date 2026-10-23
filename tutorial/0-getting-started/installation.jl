### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> title = "Installation"

using Markdown
using InteractiveUtils

# ╔═╡ 3f1c6a2e-8b4d-4e2a-9c71-0a5d2b7e1f01
md"""[← Computer-assisted proofs in nonlinear analysis](../index.html)"""

# ╔═╡ 3f1c6a2e-8b4d-4e2a-9c71-0a5d2b7e1f02
md"""
# Installation

The notebooks of this course run in [Pluto](https://plutojl.org) and share one environment (the files `Project.toml` and `Manifest.toml` at the root of the repository), so every notebook uses the same package versions.

**1. Install Julia** with [juliaup](https://julialang.org/downloads) (Julia 1.11 or later).

**2. Install Pluto**, once:
```
julia -e 'using Pkg; Pkg.add("Pluto")'
```

**3. Get the course and its packages:**
```
git clone https://github.com/capsule-math/Computer-assisted-proofs-in-nonlinear-analysis
cd Computer-assisted-proofs-in-nonlinear-analysis
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
(Without git: download the repository as a zip from GitHub and unpack it.)

**4. Start Pluto**
```
julia -e 'using Pluto; Pluto.run()'
```
and open a notebook from the course folder, e.g. `1-finite-dimensional-problems/lesson.jl`.

The first cell of every notebook activates the course environment; the first run may take a few minutes while packages compile.
"""

# ╔═╡ Cell order:
# ╟─3f1c6a2e-8b4d-4e2a-9c71-0a5d2b7e1f01
# ╟─3f1c6a2e-8b4d-4e2a-9c71-0a5d2b7e1f02
