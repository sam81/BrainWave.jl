using Weave

cd("src/")
weave("index.Rmd", informat="markdown",
      out_path = ".", doctype = "github") 
cd("../")
