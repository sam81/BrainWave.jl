
using BrainWave, Compat, MAT, Winston

chanTab = readdlm("eeglab_chan32.locs", '\t')
chanLabels = chanTab[:,4]
chanLabels[2] = "EOG1"; chanLabels[6] = "EOG2"
[chanLabels[i] = strip(chanLabels[i]) for i=1:length(chanLabels)]

fileIn = matopen("eeglab_data.set")
dset = read(fileIn, "EEG")
close(fileIn)

nEvents = length(dset["event"]["type"])
code = zeros(Int, nEvents)
idx = zeros(Int, nEvents)
sampRate = int(dset["srate"])
for i=1:nEvents
    idx[i] = int(dset["event"]["latency"][i])
    code[i] = ifelse(dset["event"]["type"][i] == "square", 1, 2)
end

data = dset["data"]
evtTab = @compat Dict{String,Any}("code" => code,
                                  "idx" => idx)

#rerefCnt!(data, 14)
idxEOG = [find(chanLabels .== "EOG1")[1], find(chanLabels .== "EOG2")[1]]
data = data .- mean(deleteSlice2D(data, idxEOG, 1), 1)

filterContinuous!(data, sampRate, "bandpass", 512, [1, 30], transitionWidth=0.2)

epochStart = -1
epochEnd = 2
segs, nRaw = segment(data, evtTab, epochStart, epochEnd,
                     sampRate, eventsList=[1])

baselineCorrect!(segs, epochStart, abs(epochStart), sampRate)

thresh = 75
toRemoveDict = findArtefactThresh(segs, thresh)
removeEpochs!(segs, toRemoveDict)

ave, nSegs = averageEpochs(segs)

tArr = ([0:size(ave["1"],2)]/sampRate) + epochStart
p = FramedPlot()
c1 = Curve(tArr, ave["1"][find(chanLabels .== "POz")[1],:])
hLine1 = LineY(0, "type", "dot", "color", "black")
vLine1 = LineX(0, "type", "dot", "color", "black")
setattr(p, "xlabel", "Time (s)")	  
setattr(p, "ylabel", "Amplitude (\\mu V)")
add(p, c1, hLine1,vLine1)
savefig(p, "plot1.png")
