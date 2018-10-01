using BrainWave, Statistics, Test

sampRate = 1000
bsDur = 0.1
sigDur = 1-0.1

bsAmp = 0.5
sigAmp = 1

nSampBs = round(Int, bsDur*sampRate)
nSampSig = round(Int, sigDur*sampRate)
nSampTot = nSampBs+nSampSig
bs = zeros(nSampBs) .+ bsAmp
sig = zeros(nSampSig) .+ sigAmp

recVals = vcat(bs, sig)'

epochDur=sigDur; preDur=bsDur; events=[1,2]; 
rec, evtTab = simulateRecording(nChans=1, dur=10, epochDur=epochDur, preDur=preDur, events=events, sampRate=sampRate)
segs, nSegs = segment(rec, evtTab, -preDur, epochDur, sampRate, eventsList=[1, 2], eventsLabelsList=["cnd1", "cnd2"])


nSegsCnd1 = size(segs["cnd1"])[3]
nSegsCnd2 = size(segs["cnd2"])[3]

for i=1:nSegsCnd1
    segs["cnd1"][:, 1:nSampBs, i] .= bsAmp
    segs["cnd1"][:, nSampBs+1:nSampTot, i] .= sigAmp
end

for i=1:nSegsCnd2
    segs["cnd2"][:, 1:nSampBs, i] .= bsAmp
    segs["cnd2"][:, nSampBs+1:nSampTot, i] .= sigAmp
end

baselineCorrect!(segs, -bsDur, bsDur, sampRate)

@test segs["cnd1"][:,nSampBs+1,1][1] .== 0.5

for i=1:nSegsCnd1
    @test mean(segs["cnd1"][:,nSampBs+1:nSampTot,i]) == 0.5
end

for i=1:nSegsCnd2
    @test mean(segs["cnd2"][:,nSampBs+1:nSampTot,i]) == 0.5
end


#test case with `bsStart` and `bsEnd`

for i=1:nSegsCnd1
    segs["cnd1"][:, 1:50, i] .= bsAmp-0.5
    segs["cnd1"][:, 51:101, i] .= bsAmp
    segs["cnd1"][:, 102:nSampTot, i] .= sigAmp
end

for i=1:nSegsCnd2
    segs["cnd2"][:, 1:50, i] .= bsAmp-0.5
    segs["cnd2"][:, 51:101, i] .= bsAmp
    segs["cnd2"][:, 102:nSampTot, i] .= sigAmp
end

baselineCorrect!(segs, -bsDur/2, 0, preDur, sampRate)

@test segs["cnd1"][:,102,1][1] .== 0.5

for i=1:nSegsCnd1
    @test mean(segs["cnd1"][:,102:nSampTot,i]) == 0.5
end

for i=1:nSegsCnd2
    @test mean(segs["cnd2"][:,102:nSampTot,i]) == 0.5
end

@test_throws(ErrorException, baselineCorrect!(segs, 0, 0, preDur, sampRate))
@test_throws(ErrorException, baselineCorrect!(segs, -0.05, -0.06, preDur, sampRate))
