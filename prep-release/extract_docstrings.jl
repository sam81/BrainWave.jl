
function extract_docstrings(fPaths, outFile::ASCIIString)
    lnsToWrite = (ByteString)[]
    for fn=1:length(fPaths)
        fPath = fPaths[fn]
        fHandle = open(fPath, "r")
        lns = readlines(fHandle)
        close(fHandle)

        startPnts = (Int)[]
        stopPnts = (Int)[]

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
    end

    fOut = open(outFile, "w")
    write(fOut, lnsToWrite)
    close(fOut)
end
