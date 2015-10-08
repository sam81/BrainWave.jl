using BrainWave, Compat
## ### simulate data of an EEG recording
## sampRate = 256
## nSec = 120
## # 16 channels, 120 seconds
## rec = rand(16, sampRate*nSec)*150
## evt = repeat([1, 2], outer=[12])
## shuffle!(evt)
## startPoints = [257:256*5:256*119]
## evtTab = @compat Dict{String,Any}("code" => evt,
##                                   "idx" => startPoints)


sampRate = 256
rec, evtTab = simulateRecording(sampRate=256, minVolt=-85, maxVolt=85)

filterContinuous!(rec, sampRate, "bandpass", 512, [1, 30], transitionWidth=0.2)

epochStart = -0.2
epochEnd = 0.5
segs, nRaw = segment(rec, evtTab, epochStart, epochEnd,
                     sampRate, eventsList=[1, 2])

baselineCorrect!(segs, epochStart, abs(epochStart), sampRate)

thresh = 75
toRemoveDict = findArtefactThresh(segs, thresh)
removeEpochs!(segs, toRemoveDict)

ave, nSegs = averageEpochs(segs)


