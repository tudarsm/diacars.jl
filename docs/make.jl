using DiaCARS
using Documenter

DocMeta.setdocmeta!(DiaCARS, :DocTestSetup, :(using DiaCARS); recursive=true)

makedocs(;
    modules=[DiaCARS],
    authors="Max Greifenstein, TU Darmstadt RSM",
    repo="https://github.com/tudarsm/diacars.jl/blob/{commit}{path}#{line}",
    sitename="DiaCARS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tudarsm.github.io/DiaCARS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Function Reference" => "functionreference.md",
        # "Development" => "development.md",
        "Examples" => "examples.md",
    ],
)

deploydocs(; repo = "github.com/tudarsm/diacars.jl")
