# Brainwave.jl Tutorial

We first load the necessary  modules. We use the `MAT.jl` module to 
read in the EEGLab dataset that is saved in MATLAB v5 format.

```julia
using BrainWave, Compat, MAT, Winston
```
we the read in the channel labels:

```julia
chanTab = readdlm("eeglab_chan32.locs", '\t')
chanLabels = chanTab[:,4]
chanLabels[2] = "EOG1"; chanLabels[6] = "EOG2"
[chanLabels[i] = strip(chanLabels[i]) for i=1:length(chanLabels)]
```
then we read in the EEGLab dataset:

```julia
fileIn = matopen("eeglab_data.set")
dset = read(fileIn, "EEG")
close(fileIn)
```

we generate the event table from the one stored with the EEGLab dataset:

```julia
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
```

we rereference the data to an average reference, excluding the EOG1 and EOG2 channels

```julia
#rerefCnt!(data, 14)
idxEOG = [find(chanLabels .== "EOG1")[1], find(chanLabels .== "EOG2")[1]]
data = data .- mean(deleteSlice2D(data, idxEOG, 1), 1)
```

we filter the recording with a bandpass filter between 1 and 30 Hz:

```julia
filterContinuous!(data, sampRate, "bandpass", 512, [1, 30], transitionWidth=0.2)
```

we epoch the recording:

```julia
epochStart = -1
epochEnd = 2
segs, nRaw = segment(data, evtTab, epochStart, epochEnd,
                     sampRate, eventsList=[1])
```

baseline correct the recording:

```julia
baselineCorrect!(segs, epochStart, abs(epochStart), sampRate)
```

detect artefacts and remove them:

```julia
thresh = 75
toRemoveDict = findArtefactThresh(segs, thresh)
removeEpochs!(segs, toRemoveDict)
```

average the epochs to obtain the ERP:

```julia
ave, nSegs = averageEpochs(segs)
```

and plot the ERP:

```julia
tArr = ([0:size(ave["1"],2)]/sampRate) + epochStart
p = FramedPlot()
c1 = Curve(tArr, ave["1"][find(chanLabels .== "POz")[1],:])
hLine1 = LineY(0, "type", "dot", "color", "black")
vLine1 = LineX(0, "type", "dot", "color", "black")
setattr(p, "xlabel", "Time (s)")	  
setattr(p, "ylabel", "Amplitude (\\mu V)")
add(p, c1, hLine1,vLine1)
savefig(p, "plot1.png")
```

![ERP](plot1.png)
