# Computer-assisted proofs in nonlinear analysis

J. B. van den Berg, M. Breden, O. Hénot, J.-P. Lessard, J. D. Mireles James.

Lessons, applications and exercises as [Pluto](https://plutojl.org) notebooks,
using [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
and [RadiiPolynomial.jl](https://github.com/OlivierHnt/RadiiPolynomial.jl).

Read online: https://capsule-math.github.io/Computer-assisted-proofs-in-nonlinear-analysis/

## Installation

1. Install Julia with [juliaup](https://julialang.org/downloads) (Julia 1.11 or later).
2. Install Pluto:
   ```
   julia -e 'using Pkg; Pkg.add("Pluto")'
   ```
3. Get the course and its packages:
   ```
   git clone https://github.com/capsule-math/Computer-assisted-proofs-in-nonlinear-analysis
   cd Computer-assisted-proofs-in-nonlinear-analysis
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```
4. Start Pluto and open any notebook, e.g. `1-finite-dimensional-problems/lesson.jl`:
   ```
   julia -e 'using Pluto; Pluto.run()'
   ```

Each notebook starts by activating the course environment (`Project.toml`,
`Manifest.toml` at the root), so all notebooks use the same package versions.

## Layout

```
Project.toml, Manifest.toml   the shared environment
N-topic/                      one folder per module, in order
  README.md                   its title (first heading)
  lesson.jl                   theory
  application.jl              worked examples
  exercise-K.jl               exercises
build.jl                      writes the course page _site/index.html
```

To add a notebook, put it in a module folder with one of these names and give
it a title (Pluto: frontmatter, or `#> title = "…"` in the header). Its first
cell should be the environment cell of the existing notebooks. To add a module,
create `N-topic/README.md`. The course page picks both up on the next push.

## Site

`.github/workflows/site.yml` runs every notebook, exports them to HTML, builds
the course page and publishes it with GitHub Pages (on push to `main`; also
monthly against the newest allowed package versions).

To preview locally:
```
julia -e 'using Pkg; Pkg.activate(temp=true); Pkg.add("PlutoSliderServer"); using PlutoSliderServer; export_directory("."; Export_output_dir="_site", Export_create_index=false)'
julia build.jl
```
then open `_site/index.html`.
