using Documenter
using FLOYao
using DocThemeIndigo

indigo = DocThemeIndigo.install(FLOYao)

makedocs(;
    modules = [FLOYao],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://QuantumBFS.github.io/FLOYao.jl",
        assets=String[indigo],
    ),
    pages = [
        "Home" => "index.md",
    ],
    repo = "https://github.com/QuantumBFS/FLOYao.jl",
    sitename = "FLOYao.jl",
)

deploydocs(; repo = "github.com/QuantumBFS/FLOYao.jl")
