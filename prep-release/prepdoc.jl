using BrainWave, Lexicon
include("extract_docstrings.jl")

extract_docstrings(["../src/BrainWave.jl"], "../docs/API.md")
#include("extract_docs.jl")
#Lexicon.save("../docs/API.md", BrainWave)

cd("../")
run(`mkdocs build`)
cd("prep-release")

