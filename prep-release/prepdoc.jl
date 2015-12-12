using BrainWave, Lexicon

Lexicon.save("../docs/API.md", BrainWave)
cd("../")
run(`mkdocs build`)
cd("prep-release")

