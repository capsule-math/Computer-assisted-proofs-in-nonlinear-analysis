#!/usr/bin/env julia
# Write _site/index.html: one section per module folder, one link per notebook.
#
#   module folders   N-name/            (sorted by N; title = first heading of its README.md)
#   notebooks        lesson.jl, application.jl, exercise-K.jl, anything else
#                    (link text = the notebook's frontmatter title)
#
# The notebook pages themselves (_site/N-name/*.html) are exported by the workflow.
# Run from the repository root:  julia build.jl
using Dates

esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

modules = sort(filter(d -> isdir(d) && occursin(r"^\d+-", d), readdir(".")),
               by = d -> parse(Int, match(r"^\d+", d).match))

function module_title(dir)
    f = joinpath(dir, "README.md")
    if isfile(f)
        m = match(r"^#\s+(.+)$"m, read(f, String))
        m === nothing || return strip(m[1])
    end
    return uppercasefirst(replace(replace(dir, r"^\d+-" => ""), "-" => " "))
end

function notebook_title(path)
    m = match(r"^#> title = \"(.*)\"$"m, read(path, String))
    return m === nothing ? splitext(basename(path))[1] : m[1]
end

# installation, lesson, application, exercise-1, exercise-2, …, then anything else
function kind(name)
    name == "installation" && return (-1, 0, "")
    name == "lesson"      && return (0, 0, "Lesson")
    name == "application" && return (1, 0, "Application")
    m = match(r"^exercise-(\d+)$", name)
    m === nothing || return (2, parse(Int, m[1]), "Exercise $(m[1])")
    return (3, 0, "")
end

sections = String[]
for dir in modules
    nbs = filter(f -> endswith(f, ".jl") && startswith(read(joinpath(dir, f), String), "### A Pluto.jl notebook ###"),
                 readdir(dir))
    items = sort([(kind(splitext(f)[1])..., f) for f in nbs], by = x -> (x[1], x[2], x[4]))
    lines = map(items) do (_, _, label, f)
        title = esc(notebook_title(joinpath(dir, f)))
        href  = "$dir/$(splitext(f)[1]).html"
        isempty(label) ? """<li><a href="$href">$title</a></li>""" :
                         """<li>$label: <a href="$href">$title</a></li>"""
    end
    number = match(r"^\d+", dir).match
    push!(sections, """<h2>$number · $(esc(module_title(dir)))</h2>
<ul>
$(join(lines, "\n"))
</ul>""")
end

page = """<!doctype html>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Computer-assisted proofs in nonlinear analysis</title>
<style>
  html{background:#fff;color:#000}
  body{font:14px/1.45 ui-monospace,Menlo,Consolas,"Liberation Mono",monospace;max-width:76rem;margin:3rem auto;padding:0 1.5rem}
  h1{font-size:1rem;font-weight:700;letter-spacing:.08em;text-transform:uppercase;margin:0 0 .5rem}
  h2{font-size:.85rem;font-weight:700;letter-spacing:.06em;text-transform:uppercase;margin:2rem 0 .5rem;padding-top:.6rem;border-top:2px solid #000}
  p{margin:0 0 2rem;max-width:48rem}
  ul{list-style:none;margin:0;padding:0}
  li{padding:.25rem 0;border-bottom:1px solid #000}
  a{color:#000}
  .muted{color:#666}
  footer{margin-top:2rem;font-size:.75rem;color:#666}
</style>
<h1>Computer-assisted proofs in nonlinear analysis</h1>
<p>J. B. van den Berg, M. Breden, O. Hénot, J.-P. Lessard, J. D. Mireles James.<br>
<span class="muted">Lessons, applications and exercises, with <a href="https://github.com/JuliaIntervals/IntervalArithmetic.jl">IntervalArithmetic.jl</a> and <a href="https://github.com/OlivierHnt/RadiiPolynomial.jl">RadiiPolynomial.jl</a>.<br><br>
<a href="https://github.com/capsule-math/Computer-assisted-proofs-in-nonlinear-analysis">source</a> ·
<a href="https://capsule-math.github.io/index/">proof packages</a></span></p>
$(join(sections, "\n"))
<footer>generated $(Dates.format(today(), "yyyy-mm-dd"))</footer>
"""

mkpath("_site")
write(joinpath("_site", "index.html"), page)
println(length(modules), " modules")
