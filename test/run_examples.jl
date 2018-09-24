using BrainWave

fToProcess = ["BrainWave.jl"]
for fName in fToProcess
#fName = fToProcess[1]
    fIn = open(string("../src/", fName), "r")
    fOut = open(string("_run_examples_", fName), "w")
    lns = readlines(fIn)
    idxStart = (Int)[]
    idxStop = (Int)[]
    for i=1:length(lns)
        if lns[i] == "```julia"
            push!(idxStart, i+1)
        elseif lns[i] == "```"
            push!(idxStop, i-1)
        end
    end
    for i=1:length(idxStart)
        theseLines = lns[idxStart[i]:idxStop[i]]
        for currLine in theseLines
            write(fOut, string(currLine, "\n"))
        end
    end
    close(fIn); close(fOut)
end

for fName in fToProcess
    fNameTest = string("_run_examples_", fName)
    include(fNameTest)
end

