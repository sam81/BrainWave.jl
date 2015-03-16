using ElectroJulia, Lexicon


Lexicon.save("../docs/API.md", ElectroJulia)
cd("../")
run(`mkdocs build`)
cd("prep-release")
