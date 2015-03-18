using Weave

tangle("tutorial.md", out_path=:pwd, informat="markdown")
reload("tutorial.jl")
