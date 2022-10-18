using Documenter
using FLOYao
using YaoAPI
using DocThemeIndigo
using YaoArrayRegister

indigo = DocThemeIndigo.install(FLOYao)

const PAGES = [
    "Home" => "index.md",
    "Example: VQE for the TFIM" => "vqe_example.md",
    "Features" => "features/features.md",
    "Supported Gates" => "features/supported_gates.md",
    "Mathematical Background" => "background.md",
    "Adding custom gates" => "adding_gates.md",
    "Known restrictions" => "known_restrictions.md",
]

makedocs(;
    modules = [FLOYao,YaoAPI],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://QuantumBFS.github.io/FLOYao.jl",
        assets=String[indigo],
    ),
    pages = PAGES,
    repo = "https://github.com/QuantumBFS/FLOYao.jl",
    sitename = "FLOYao.jl",
)

deploydocs(; repo = "github.com/QuantumBFS/FLOYao.jl")
