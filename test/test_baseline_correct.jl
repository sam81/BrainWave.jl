using BrainWave, Statistics, Test

sampRate = 1000
bsDur = 0.01
sigDur = 1-0.01

bsAmp = 0.5
sigAmp = 1

nSampBs = round(Int, bsDur*sampRate)
nSampSig = round(Int, sigDur*sampRate)
bs = zeros(nSampBs) .+ bsAmp
sig = zeros(nSampSig) .+ sigAmp

recVals = vcat(bs, sig)'

epochDur=sigDur; preDur=bsDur; events=[1,2]; 
rec, evtTab = simulateRecording(nChans=1, dur=10, epochDur=epochDur, preDur=preDur, events=events, sampRate=sampRate)
segs, nSegs = segment(rec, evtTab, -preDur, epochDur, sampRate, eventsList=[1, 2], eventsLabelsList=["cnd1", "cnd2"])


nSegsCnd1 = size(segs["cnd1"])[3]
nSegsCnd2 = size(segs["cnd2"])[3]

for i=1:nSegsCnd1
    segs["cnd1"][:, 1:100, i] .= bsAmp
    segs["cnd1"][:, 101:1000, i] .= sigAmp
end

for i=1:nSegsCnd2
    segs["cnd2"][:, 1:100, i] .= bsAmp
    segs["cnd2"][:, 101:1000, i] .= sigAmp
end

baselineCorrect!(segs, -bsDur, bsDur, sampRate)

@test segs["cnd1"][:,101,1][1] .== 0.5

for i=1:nSegsCnd1
    @test mean(segs["cnd1"][:,101:1000,i]) == 0.5
end

for i=1:nSegsCnd2
    @test mean(segs["cnd2"][:,101:1000,i]) == 0.5
end
