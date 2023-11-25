using Quaternions
using Documenter
using HalfIntegers

DocMeta.setdocmeta!(Quaternions, :DocTestSetup, :(using Quaternions); recursive=true)

makedocs(;
    modules=[Quaternions],
    repo="https://github.com/JuliaGeometry/Quaternions.jl/blob/{commit}{path}#{line}",
    sitename="Quaternions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGeometry.github.io/Quaternions.jl",
        assets = ["assets/custom.css", "assets/favicon.ico"],
        repolink = "https://github.com/JuliaGeometry/Quaternions.jl"
    ),
    pages=[
        "Home" => "index.md",
        "APIs" => "api.md",
        "Examples" => [
            "examples/basics.md",
            "examples/type_parameter.md",
            "examples/rotations.md",
            "examples/dual_quaternions.md"
        ],
    ],
    warnonly=true
)

deploydocs(;
    repo="github.com/JuliaGeometry/Quaternions.jl", devbranch="main",
)
