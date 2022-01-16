using QuantumMaterials
using Documenter

DocMeta.setdocmeta!(QuantumMaterials, :DocTestSetup, :(using QuantumMaterials); recursive = true)

makedocs(;
    modules = [QuantumMaterials],
    authors = "Kun Chen",
    repo = "https://github.com/numericaleft/QuantumMaterials.jl/blob/{commit}{path}#{line}",
    sitename = "QuantumMaterials.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://numericaleft.github.io/QuantumMaterials.jl",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Background" => [
            "background/silicon_empirical_pseudopotential.md",
        ],
    ]
)

deploydocs(;
    repo = "github.com/numericalEFT/QuantumMaterials.jl",
    devbranch = "main"
)
