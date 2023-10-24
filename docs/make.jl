using MyOscar
using Documenter

DocMeta.setdocmeta!(MyOscar, :DocTestSetup, :(using MyOscar); recursive=true)

makedocs(;
    modules=[MyOscar],
    authors="LucaRemke <luca.remke@idsr.uni-stuttgart.de> and contributors",
    repo="https://github.com/LucaRemke/MyOscar.jl/blob/{commit}{path}#{line}",
    sitename="MyOscar.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://LucaRemke.github.io/MyOscar.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LucaRemke/MyOscar.jl",
    devbranch="main",
)
