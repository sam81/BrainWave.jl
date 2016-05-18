fPaths = ["../src/BrainWave.jl"]

fHandle = open(fPaths[1], "r")
lns = readlines(fHandle)
close(fHandle)

startPnts = (Int)[]
stopPnts = (Int)[]
lnsToWrite = (ByteString)[]
for i=1:length(lns)
    if strip(lns[i]) == "@doc doc\"\"\""
        push!(startPnts, i)
    elseif strip(lns[i]) == "\"\"\"->"
        push!(stopPnts, i)
    end
end

for s=1:length(startPnts)
    for nl=1:stopPnts[s]-startPnts[s]-1
        push!(lnsToWrite, lns[startPnts[s]+nl])
    end
    push!(lnsToWrite, "\n")
end


fOut = open("../docs/API.md", "w")
write(fOut, lnsToWrite)
close(fOut)
