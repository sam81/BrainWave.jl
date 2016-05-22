using Base.Test
@everywhere using BrainWave

sampRate = 8192; nTaps=512
nChans=8
rec, evtTab = simulateRecording(nChans=nChans, dur=500, sampRate=sampRate)
rec2 = deepcopy(rec)
chans = collect(1:nChans)


filterContinuous!(rec, sampRate, "highpass", nTaps, [30], channels=chans, transitionWidth=0.2);


sharedRec = convert(SharedArray, rec2);
filterContinuous!(sharedRec, sampRate, "highpass", nTaps, [30], channels=chans, transitionWidth=0.2);
rec2 = convert(Array, sharedRec)

@test unique(rec2 .== rec)[1] .== true
