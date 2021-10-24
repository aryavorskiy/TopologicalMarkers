using Documenter
using TopologicalMarkers

push!(LOAD_PATH,"../src/")

pages = [
    "Home" => "index.md",
    "Calculate" => [
        "Linear operators" => "operators.md",
        "Unitary evolution" => "evolution.md"
    ],
    "Visualize" => "visual.md",
    "Library" => "scope.md"
]

makedocs(
    modules = [TopologicalMarkers],
    sitename = "TopologicalMarkers",
    pages = pages
)

deploydocs(
    repo = "github.com/aryavorskiy/TopologicalMarkers.git",
)