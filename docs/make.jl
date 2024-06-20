using StaircaseShenanigans
using Documenter

DocMeta.setdocmeta!(StaircaseShenanigans, :DocTestSetup, :(using StaircaseShenanigans); recursive=true)

makedocs(;
    modules=[StaircaseShenanigans],
    authors="Josef Bisits <jbisits@gmail.com>",
    sitename="StaircaseShenanigans.jl",
    format=Documenter.HTML(;
        canonical="https://jbisits.github.io/StaircaseShenanigans.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jbisits/StaircaseShenanigans.jl",
    devbranch="main",
)
