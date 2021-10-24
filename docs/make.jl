using Documenter
using TopologicalMarkers

push!(LOAD_PATH,"../src/")

pages = [
    "Home" => "index.md",
    "Scope" => "scope.md"
]

makedocs(
    modules = [TopologicalMarkers],
    sitename = "TopologicalMarkers",
    pages = pages
)