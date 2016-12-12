using Base.Test
@everywhere using DistributedArrays
@everywhere using BrainWave
srand(1234)
sampRate = 8192; nTaps=512
nChans=8
rec, evtTab = simulateRecording(nChans=nChans, dur=500, sampRate=sampRate)
rec2 = deepcopy(rec)
chans = collect(1:nChans)

filterContinuous!(rec, sampRate, "highpass", nTaps, [30], channels=chans, transitionWidth=0.2)

nProcsToUse = min(length(workers()), size(rec, 1))
workersToUse = workers()[1:nProcsToUse]
println(workersToUse)
dRec = distribute(rec2, procs=workersToUse, dist=[nProcsToUse, 1])
filterContinuous!(dRec, sampRate, "highpass", nTaps, [30], workersToUse, channels=chans, transitionWidth=0.2)

rec2[:,:] = convert(Array, dRec)
@test unique(rec2 .== rec)[1] .== true
