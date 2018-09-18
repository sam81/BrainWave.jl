module BrainWave

export averageAverages, averageEpochs, averageEpochsIterativeWeighted, baselineCorrect!,
deleteSlice2D, deleteSlice3D, detrendEEG!, epochVariance, filterContinuous!, #_centered,
fftconvolve, findArtefactThresh, findExtremum, findInflections, findPeaks, findTroughs, FMP, FMPIterativeWeighted, getACF, getACF2, getAutocorrelogram, getAutocorrelogram2,
getSNR, getSNR2, getPhaseSpectrum, getSpectrogram, getSpectrum,
iterativeWeightedAverage,
meanERPAmplitude, mergeEventTableCodes!, nextPowTwo,
removeEpochs!, removeSpuriousTriggers!, rerefCnt!, RMS,
segment, simulateRecording,
toRawEEG, RawEEG,
input, plotRawEEG
#filterContinuousParallel!

#getNoiseSidebands, #chainSegments,#getFRatios,


using DataFrames, Distributed, Distributions, FFTW, DSP, PyCall
using DistributedArrays, SharedArrays
#import Compat.String
import PyPlot; const plt = PyPlot
using DocStringExtensions 

#pyinitialize("python3")

#@pyimport scipy.signal as scisig
const scisig = PyNULL()

function __init__()
    copy!(scisig, pyimport_conda("scipy.signal", "scipy"))
end

#include("findABRPeaks.jl")
include("AbstractionLayer.jl")


"""
Perform a weighted average of a list of averages. The weight of
each average in the list is determined by the number of segments
from which it was obtained.

$(SIGNATURES)

##### Arguments

* `aveListArray::Array{Dict{String, Array{Real, 2}}}`: The list of averages for each experimental condition.
* `nSegments::Array{Dict{String, Integer}}`: The number of epochs on which each average is based.

##### Returns

* `weightedAve::Dict{String,Array{Real,2}}`: The weighted average of the averages in the list.

##### Examples

```julia
    epochDur=0.5; preDur=0.15; events=[1,2]; sampRate=256;
    rec, evtTab = simulateRecording(dur=120, epochDur=epochDur, preDur=preDur, events=events)
    segs, nRaw = segment(rec, evtTab, -preDur, epochDur, sampRate)
    ave, nSegs = averageEpochs(segs)

    rec2, evtTab2 = simulateRecording(dur=120, epochDur=epochDur, preDur=preDur, events=events)
    segs2, nRaw2 = segment(rec2, evtTab2, -preDur, epochDur, sampRate)
    ave2, nSegs2 = averageEpochs(segs2)

    aveList = [ave, ave2]; nCleanByBlock = [nRaw, nRaw2]
    aveAll, nSegsAll = averageAverages(aveList, nCleanByBlock)
```

"""
function averageAverages(aveList::Array{Dict{String, Array{T, 2}}}, nSegments::Array{Dict{String, P}}) where {T<:Real, P<:Integer}
    eventList = collect(keys(aveList[1]))
    weightedAve = Dict{String,Array{eltype(aveList[1][eventList[1]]),2}}() #(String => Array{eltype(aveList[1][eventList[1]]),2})[]
    nSegsSum = Dict{String,Int}()#(String => Int)[]
    for i=1:length(eventList)
        event = eventList[i]
        nSegsSum[event] = 0
        for j=1:length(aveList)
            nSegsSum[event] = nSegsSum[event] + nSegments[j][event]
        end
    end

    for i=1:length(eventList)
        event = eventList[i]
        weightedAve[event] = zeros(eltype(aveList[1][eventList[1]]), size(aveList[1][event]))
        for j=1:length(aveList)
            weightedAve[event] = weightedAve[event] + aveList[j][event] * (nSegments[j][event]/nSegsSum[event])
        end
    end

    return weightedAve, nSegsSum
end

"""
Average the epochs of a segmented recording.

$(SIGNATURES)

##### Arguments

* `rec::Dict{String,Array{T,3}}`: Dictionary containing the segmented recordings for each condition.
        The segmented recordings consist of 3-dimensional arrays (n_channels x n_samples x n_epochs).

##### Returns

* `ave::Dict{String,Array{Real,2}}`: The averaged epochs for each condition.
* `n_segs::Dict{String,Integer}`: The number of epochs averaged for each condition.

##### Examples

```julia
    epochDur=0.5; preDur=0.15; events=[1,2]; sampRate=256;
    rec, evtTab = simulateRecording(dur=120, epochDur=epochDur, preDur=preDur, events=events)
    segs, nRaw = segment(rec, evtTab, -preDur, epochDur, sampRate)
    ave, nSegs = averageEpochs(segs)
```

##### averageEpochs{T<:Real}(rec::Array{T,3})
"""
function averageEpochs(rec::Array{T,3}) where {T<:Real}

    nSegs = size(rec)[3]
    ave = dropdims(mean(rec, dims=3), dims=3)

    return ave, nSegs
end

function averageEpochs(rec::Dict{String,Array{T,3}}) where {T<:Real}


    eventList = collect(keys(rec))
    ave = Dict{String,Array{eltype(rec[eventList[1]]),2}}()
    nSegs = Dict{String,Int}()
    for i=1:length(eventList)
        #nSegs[eventList[i]] = size(rec[eventList[i]])[3]
        #ave[eventList[i]] = mean(rec[eventList[i]], 3)[:,:,1]
        ave[eventList[i]], nSegs[eventList[i]] = averageEpochs(rec[eventList[i]])
    end
    return ave, nSegs
end

"""
Average the epochs of a segmented recording using iterative weighted averaging algorithm (see References below).

$(SIGNATURES)

##### Arguments

* `rec::Dict{String,Array{T,3}}`: Dictionary containing the segmented recordings for each condition.
        The segmented recordings consist of 3-dimensional arrays (n_channels x n_samples x n_epochs).
* `noiseWinStart`::Real: Time in seconds at which the noise estimate should start relative to the start of the epoch.
* `noiseWinStop`::Real: Time in seconds at which the noise estimate should stop relative to the start of the epoch.
* `preDur::Real`: Duration of recording before the experimental event, in seconds.
* `sampRate::Real`: The samplig rate of the EEG recording.
* `noiseEstimate::String`: if `global` the estimate of the noise used to weight individual
segments is derived from all the channels. If `byChannel` an noise estimate is derived for each
channel and segments are weighted differently for each channel depending on the channel noise
estimate.

##### Returns

* `ave::Dict{String,Array{Real,2}}`: The averaged epochs for each condition.

##### Examples

```julia
    epochDur=0.5; preDur=0.15; events=[1,2]; sampRate=256;
    rec, evtTab = simulateRecording(dur=120, epochDur=epochDur, preDur=preDur, events=events)
    segs, nRaw = segment(rec, evtTab, -preDur, epochDur, sampRate)
    ave = averageEpochsIterativeWeighted(segs)
```

##### References

* Riedel, H., Granzow, M., & Kollmeier, B. (2001). Single-sweep-based methods to improve the quality of auditory brain stem responses Part II: Averaging methods. Z Audiol, 40(2), 62–85.

"""
function averageEpochsIterativeWeighted(rec::Dict{String,Array{T,3}}; noiseEstimate::String="global") where {T<:Real}
    eventList = collect(keys(rec))
    ave = Dict{String,Array{eltype(rec[eventList[1]]),2}}()
    for i=1:length(eventList)
        ave[eventList[i]] = iterativeWeightedAverage(rec[eventList[i]], noiseEstimate=noiseEstimate)
    end
    return ave
end


function averageEpochsIterativeWeighted(rec::Dict{String,Array{T,3}}, noiseWinStart::Real, noiseWinStop::Real, preDur::Real, sampRate::Real; noiseEstimate::String="global") where {T<:Real}
    eventList = collect(keys(rec))
    ave = Dict{String,Array{eltype(rec[eventList[1]]),2}}()
    for i=1:length(eventList)
        ave[eventList[i]] = iterativeWeightedAverage(rec[eventList[i]], noiseWinStart, noiseWinStop, preDur, sampRate, noiseEstimate=noiseEstimate)
    end
    return ave
end

"""
Perform baseline correction by subtracting the average pre-event
voltage from each channel of a segmented recording.

$(SIGNATURES)

##### Arguments

* `rec::Dict{String,Array{T,3}}`: The segmented recording.
* `baselineStart::Real`: Start time of the baseline window relative to the event onset, in seconds.
                          The absolute value of `baselineStart` cannot be greater than `preDur`.
                          In practice `baselineStart` allows you to define a baseline window shorter
                          than the time window before the experimental event (`preDur`).
* `preDur::Real`: Duration of recording before the experimental event, in seconds.
* `sampRate::Integer`: The samplig rate of the EEG recording.

##### Examples

```julia
    epochDur=0.5; preDur=0.15; events=[1,2]; sampRate=256;
    rec, evtTab = simulateRecording(dur=120, epochDur=epochDur, preDur=preDur, events=events)
    segs, nRaw = segment(rec, evtTab, -preDur, epochDur, sampRate)

    #baseline window has the same duration of preDur
    baselineCorrect!(segs, -preDur, preDur, sampRate)
    #now with a baseline shorter than preDur
    baselineCorrect!(segs, -0.15, preDur, sampRate)
```

"""
function baselineCorrect!(rec::Dict{String,Array{T,3}}, baselineStart::Real, preDur::Real, sampRate::Integer) where {T<:Real}
    eventList = collect(keys(rec))
    epochStartSample = round(Int, preDur*sampRate)
    baselineStartSample = round(Int, (epochStartSample+1) - abs(round(baselineStart*sampRate)))

    for i=1:length(eventList) #for each event
        for j=1:size(rec[eventList[i]])[3] #for each epoch
            for k=1: size(rec[eventList[i]])[1] #for each electrode
                thisBaseline = mean(rec[eventList[i]][k,baselineStartSample:epochStartSample,j])
                rec[eventList[i]][k,:,j] = rec[eventList[i]][k,:,j] .- thisBaseline
            end
        end
    end


end

function baselineCorrect!(rec::Dict{String,Array{T,3}}, baselineStart::Real, baselineEnd::Real, preDur::Real, sampRate::Integer) where {T<:Real}
    eventList = collect(keys(rec))
    epochStartSample = round(Int, preDur*sampRate)
    baselineStartSample = round(Int, (epochStartSample+1) - abs(round(baselineStart*sampRate)))
    baselineEndSample = round(Int, (epochStartSample+1) - abs(round(baselineEnd*sampRate)))

    for i=1:length(eventList) #for each event
        for j=1:size(rec[eventList[i]])[3] #for each epoch
            for k=1: size(rec[eventList[i]])[1] #for each electrode
                thisBaseline = mean(rec[eventList[i]][k,baselineStartSample:baselineEndSample,j])
                rec[eventList[i]][k,:,j] = rec[eventList[i]][k,:,j] .- thisBaseline
            end
        end
    end


end


"""
Delete a row or a column from a 2-dimensional array.

$(SIGNATURES)

##### Arguments

* `x::AbstractMatrix{T}`: the 2-dimensional array from which rows of columns should be deleted.
* `toRemove::Union{Integer, AbstractVector{Integer}}`: an integer or a list of integers indicating the rows/columns to remove.
* `axis::Integer`: an integer indicating whether rows or columns should be removed. 1 corresponds to rows, 2 corresponds to columns.

##### Returns

* `x::AbstractMatrix{T}`: a 2-dimensional array with the selected row or columns removed.

##### Examples

```julia
    x = [1 2 3 4;
         5 6 7 8;
         9 10 11 12;
         13 14 15 16
        ]

    #remove first row
    isequal(deleteSlice2D(x, 1, 1), x[2:end,:])
    #remove rows 1 and 4
    isequal(deleteSlice2D(x, [1,4], 1), x[2:3,:])
    # remove columns 1 and 4
    isequal(deleteSlice2D(x, [1,4], 2), x[:,2:3])
```

"""
function deleteSlice2D(x::AbstractMatrix{T}, toRemove::Union{P, AbstractVector{P}}, axis::Integer) where {T<:Any, P<:Integer}
    if in(axis, [1,2]) == false
        error("axis must be either 1 (rows), or 2 (columns)")
    end
    if length(toRemove) > 1
        toRemove = sort(toRemove)
    end
    if axis == 1
        for i=1:length(toRemove)
            #x = x[[setdiff(1:size(x, 1), toRemove[i]-i+1)],:]
            x = x[collect(setdiff(1:size(x, 1), toRemove[i]-i+1)),:]
        end
    elseif axis == 2
        for i=1:length(toRemove)
            #x = x[:, [setdiff(1:size(x, 2), toRemove[i]-i+1)]]
            x = x[:, collect(setdiff(1:size(x, 2), toRemove[i]-i+1))]
        end
    end
    return(x)
end

"""
Delete a slice from a 3-dimensional array.

$(SIGNATURES)

##### Arguments

* `x::Array{T, 3}`: the 3-dimensional array from which rows of columns should be deleted.
* `toRemove::Union{Integer, AbstractVector{Integer}}`: an integer or a list of integers indicating the slices to remove.
* `axis::Integer`: an integer indicating the axis along which slices should be removed.

##### Returns

* `x::Array{T, 3}`: a 3-dimensional array with the selected row or columns removed.

##### Examples

```julia
    x = reshape(collect(1:27), 3,3,3)
    deleteSlice3D(x, 2, 1)
    deleteSlice3D(x, [2,3], 3)
    isequal(deleteSlice3D(x, [2,3], 3), x[:,:, [1]])
```

"""
function deleteSlice3D(x::Array{T,3}, toRemove::Union{P, AbstractVector{P}}, axis::Integer) where {T<:Any, P<:Integer}

    if in(axis, [1,2,3]) == false
        error("axis must be either 1, 2, or 3")
    end
    if length(toRemove) > 1
        toRemove = sort(toRemove)
    end
    if axis == 1
        for i=1:length(toRemove)
            x = x[collect(setdiff(1:size(x, 1), toRemove[i]-i+1)),:,:]
        end
    elseif axis == 2
        for i=1:length(toRemove)
            x = x[:, collect(setdiff(1:size(x, 2), toRemove[i]-i+1)),:]
        end
    elseif axis == 3
        for i=1:length(toRemove)
            x = x[:, :, collect(setdiff(1:size(x, 3), toRemove[i]-i+1))]
        end
    end

    return x
end

"""
Remove the mean value from each channel of an EEG recording.

$(SIGNATURES)

##### Arguments

* `rec::AbstractMatrix{T}`: The EEG recording.

##### Examples

```julia
    x = [1 2 3; 4 5 6]
    detrendEEG!(x)
```

"""
function detrendEEG!(rec::AbstractMatrix{T}) where {T<:Real}

    nChannels = size(rec)[1]
    for i=1:nChannels
        rec[i,:] = rec[i,:] .- mean(rec[i,:])
    end

end

"""
$(SIGNATURES)
"""
function epochVariance(rec::Dict{String,Array{T,3}}, winStart::Real, winStop::Real, preDur::Real, sampRate::Real) where {T<:Real}
    eventList = collect(keys(rec))
    epochVar = Dict{String, Real}()
    for i=1:length(eventList)
        epochVar[eventList[i]] = epochVariance(rec[eventList[i]], winStart, winStop, preDur, sampRate)
    end
    return epochVar
end

function epochVariance(sweeps::Array{T,3}, winStart::Real, winStop::Real, preDur::Real, sampRate::Real) where {T<:Real}
    nSweeps = size(sweeps)[3]
    nSamp = size(sweeps)[2]
    nChans = size(sweeps)[1]

    winStartSample = round(Int, (winStart+preDur)*sampRate)+1
    winStopSample = round(Int, (winStop+preDur)*sampRate)
    if winStartSample < 1
        winStartSample = 1
    end
    if winStopSample > nSamp
        winStopSample = nSamp
    end
    ## println(winStartSample)
    ## println(winStopSample)
    epochVar = mean(var(sweeps[:,winStartSample:winStopSample,:], 3))
end

"""
Filter a continuous EEG recording.

$(SIGNATURES)

##### Arguments

* `rec::AbstractMatrix{Real}`: The nChannelsXnSamples array with the EEG recording.
* `sampRate::Integer`: The EEG recording sampling rate.
* `filterType::String`:  The filter type., one of "lowpass", "highpass", or "bandpass".
* `nTaps::Integer`: The number of filter taps.
* `cutoffs::Union{Real, AbstractVector{Real}}`:: The filter cutoffs. If "filterType" is "lowpass" or "highpass"
        the "cutoffs" array should contain a single value. If "filterType"
        is bandpass the "cutoffs" array should contain the lower and
        the upper cutoffs in increasing order.
* `channels::Union{Q, AbstractVector{Q}}`: The channel number, or the list of channels numbers that should be filtered.
* `transitionWidth::Real`: The width of the filter transition region, normalized between 0-1.
        For a lower cutoff the nominal transition region will go from
        `(1-transitionWidth)*cutoff` to `cutoff`. For a higher cutoff
        the nominal transition region will go from cutoff to
        `(1+transitionWidth)*cutoff`.

##### Examples

```julia
    sampRate = 2048; nTaps=512
    rec, evtTab = simulateRecording(nChans=4, dur=120, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "highpass", nTaps, [30], channels=[1,2,3,4], transitionWidth=0.2)

    #using the parallel processing version of the function with a SharedArray
    @everywhere using BrainWave
    rec, evtTab = simulateRecording(nChans=4, dur=120, sampRate=sampRate);
    #to exploit parallel processing you need to call julia with -p `n` parameter, or use addprocs(`n`)
    #where n is the number of processes you want to use
    #convert recording to SharedArray
    sharedRec = convert(SharedArray, rec)
    #call the filtering function
    filterContinuous!(rec, sampRate, "highpass", nTaps, [30], channels=[1,2,3,4], transitionWidth=0.2)
    #convert back to Array type
    rec = convert(Array, sharedRec)

    #using the parallel processing version of the function with a DistributedArray
    @everywhere using BrainWave, DistributedArrays
    rec, evtTab = simulateRecording(nChans=4, dur=120, sampRate=sampRate);
    #to exploit parallel processing you need to call julia with -p `n` parameter, or use addprocs(`n`)
    #where n is the number of processes you want to use
    nProcsToUse = min(length(workers()), size(rec, 1))
    workersToUse = workers()[1:nProcsToUse]
    #convert recording to DistributedArray
    dRec = distribute(rec, procs=workersToUse, dist=[nProcsToUse, 1])
    #call the filtering function, note that you need to pass the ids of the workers to use
    filterContinuous!(dRec, sampRate, "highpass", nTaps, [30], workersToUse, channels=[1,2,3,4], transitionWidth=0.2)
    #convert back to Array type
    rec[:,:] = convert(Array, dRec)
```

"""
function filterContinuous!(rec::AbstractMatrix{T},
                           sampRate::Integer,
                           filterType::String,
                           nTaps::Integer,
                           cutoffs::Union{P, AbstractVector{P}};
                           channels::Union{Q, AbstractVector{Q}}=collect(1:size(rec,1)),
                           transitionWidth::Real=0.2) where {T<:Real, P<:Real, Q<:Integer}
    if filterType == "lowpass"
        f3 = cutoffs[1]
        f4 = cutoffs[1] * (1+transitionWidth)
        f3 = (f3*2) / sampRate
        f4 = (f4*2) / sampRate
        f = [0, f3, f4, 1]
        m = [1, 1, 0.00003, 0]
    elseif filterType == "highpass"
        f1 = cutoffs[1] * (1-transitionWidth)
        f2 = cutoffs[1]
        f1 = (f1*2) / sampRate
        f2 = (f2*2) / sampRate
        f = [0, f1, f2, 0.999999, 1] #high pass
        m = [0, 0.00003, 1, 1, 0]
    elseif filterType == "bandpass"
        f1 = cutoffs[1] * (1-transitionWidth)
        f2 = cutoffs[1]
        f3 = cutoffs[2]
        f4 = cutoffs[2] * (1+transitionWidth)
        f1 = (f1*2) / sampRate
        f2 = (f2*2) / sampRate
        f3 = (f3*2) / sampRate
        f4 = (f4*2) / sampRate
        f = [0, f1, f2, ((f2+f3)/2), f3, f4, 1]
        m = [0, 0.00003, 1, 1, 1, 0.00003, 0]
    end

    b = convert(Array{eltype(rec),1}, scisig[:firwin2](nTaps,f,m))
    nChannels = size(rec)[1]

    for i=1:nChannels
        if in(i, channels) == true
            if VERSION < v"0.5-"
                rec[i,:] = fftconvolve(reshape(rec[i,:], size(rec[i,:], 2)), b, "same")
                rec[i,:] = reverse(fftconvolve(reverse(reshape(rec[i,:], size(rec[i,:], 2)), dims=1), b, "same"), dims=1)
            else
                rec[i,:] = fftconvolve(rec[i,:], b, "same")
                rec[i,:] = reverse(fftconvolve(reverse(rec[i,:], dims=1), b, "same"), dims=1)
            end
           
        end
    end
    #return rec
end

#####################################
# `filterContinuous!` for SharedArray
#####################################
function filterContinuous!(rec::SharedArray{T, 2},
                           sampRate::Integer,
                           filterType::String,
                           nTaps::Integer,
                           cutoffs::Union{P, AbstractVector{P}};
                           channels::Union{Q, AbstractVector{Q}}=collect(1:size(rec,1)),
                           transitionWidth::Real=0.2) where {T<:Real, P<:Real, Q<:Integer}
    if filterType == "lowpass"
        f3 = cutoffs[1]
        f4 = cutoffs[1] * (1+transitionWidth)
        f3 = (f3*2) / sampRate
        f4 = (f4*2) / sampRate
        f = [0, f3, f4, 1]
        m = [1, 1, 0.00003, 0]
    elseif filterType == "highpass"
        f1 = cutoffs[1] * (1-transitionWidth)
        f2 = cutoffs[1]
        f1 = (f1*2) / sampRate
        f2 = (f2*2) / sampRate
        f = [0, f1, f2, 0.999999, 1] #high pass
        m = [0, 0.00003, 1, 1, 0]
    elseif filterType == "bandpass"
        f1 = cutoffs[1] * (1-transitionWidth)
        f2 = cutoffs[1]
        f3 = cutoffs[2]
        f4 = cutoffs[2] * (1+transitionWidth)
        f1 = (f1*2) / sampRate
        f2 = (f2*2) / sampRate
        f3 = (f3*2) / sampRate
        f4 = (f4*2) / sampRate
        f = [0, f1, f2, ((f2+f3)/2), f3, f4, 1]
        m = [0, 0.00003, 1, 1, 1, 0.00003, 0]
    end

    b = convert(Array{eltype(rec),1}, scisig[:firwin2](nTaps,f,m))
    #b = ones(Float32,nTaps)
    nChannels = size(rec)[1]
    ## if channels == nothing
    ##     channels = [1:nChannels]
    ## end

    @sync @distributed for i=1:nChannels
        if in(i, channels) == true
            if VERSION < v"0.5-"
                rec[i,:] = fftconvolve(reshape(rec[i,:], size(rec[i,:], 2)), b, "same")
                rec[i,:] = reverse(fftconvolve(reverse(reshape(rec[i,:], size(rec[i,:], 2)), dims=1), b, "same"), dims=1)
            else
                rec[i,:] = fftconvolve(rec[i,:], b, "same")
                rec[i,:] = reverse(fftconvolve(reverse(rec[i,:], dims=1), b, "same"), dims=1)
            end
           
        end
    end
    #return rec
end

#################################
# `filterContinuous!` for DArray
#################################
function filterContinuous!(rec::DArray{T},
                           sampRate::Integer,
                           filterType::String,
                           nTaps::Integer,
                           cutoffs::Union{P, AbstractVector{P}},
                           workersToUse::AbstractVector{W};
                           channels::Union{Q, AbstractVector{Q}}=collect(1:size(rec,1)),
                           transitionWidth::Real=0.2) where {T<:Real, P<:Real, Q<:Integer, W<:Integer}
    if filterType == "lowpass"
        f3 = cutoffs[1]
        f4 = cutoffs[1] * (1+transitionWidth)
        f3 = (f3*2) / sampRate
        f4 = (f4*2) / sampRate
        f = [0, f3, f4, 1]
        m = [1, 1, 0.00003, 0]
    elseif filterType == "highpass"
        f1 = cutoffs[1] * (1-transitionWidth)
        f2 = cutoffs[1]
        f1 = (f1*2) / sampRate
        f2 = (f2*2) / sampRate
        f = [0, f1, f2, 0.999999, 1] #high pass
        m = [0, 0.00003, 1, 1, 0]
    elseif filterType == "bandpass"
        f1 = cutoffs[1] * (1-transitionWidth)
        f2 = cutoffs[1]
        f3 = cutoffs[2]
        f4 = cutoffs[2] * (1+transitionWidth)
        f1 = (f1*2) / sampRate
        f2 = (f2*2) / sampRate
        f3 = (f3*2) / sampRate
        f4 = (f4*2) / sampRate
        f = [0, f1, f2, ((f2+f3)/2), f3, f4, 1]
        m = [0, 0.00003, 1, 1, 1, 0.00003, 0]
    end

    b = convert(Array{eltype(rec),1}, scisig[:firwin2](nTaps,f,m))
    nChannels = size(rec)[1]

    ## nProcsToUse = min(length(workers()), size(rec, 1))
    ## dRec = distribute(rec, procs=workers()[1:nProcsToUse], dist=[nProcsToUse,1])
    #println(workersToUse)
    @sync for w in workersToUse 
        @spawnat w doWorkFilterContinuousParallel!(rec, b, channels=channels)
    end
    #rec[:,:] = convert(Array, dRec)
    #close(dRec)
    #darray_closeall()
end

## #for DistributedArray
## function filterContinuousParallel!{T<:Real, P<:Real, Q<:Integer}(rec::AbstractMatrix{T},
##                                                          sampRate::Integer,
##                                                          filterType::String,
##                                                          nTaps::Integer,
##                                                          cutoffs::Union{P, AbstractVector{P}};
##                                                          channels::Union{Q, AbstractVector{Q}}=collect(1:size(rec,1)),
##                                                          transitionWidth::Real=0.2)
##     if filterType == "lowpass"
##         f3 = cutoffs[1]
##         f4 = cutoffs[1] * (1+transitionWidth)
##         f3 = (f3*2) / sampRate
##         f4 = (f4*2) / sampRate
##         f = [0, f3, f4, 1]
##         m = [1, 1, 0.00003, 0]
##     elseif filterType == "highpass"
##         f1 = cutoffs[1] * (1-transitionWidth)
##         f2 = cutoffs[1]
##         f1 = (f1*2) / sampRate
##         f2 = (f2*2) / sampRate
##         f = [0, f1, f2, 0.999999, 1] #high pass
##         m = [0, 0.00003, 1, 1, 0]
##     elseif filterType == "bandpass"
##         f1 = cutoffs[1] * (1-transitionWidth)
##         f2 = cutoffs[1]
##         f3 = cutoffs[2]
##         f4 = cutoffs[2] * (1+transitionWidth)
##         f1 = (f1*2) / sampRate
##         f2 = (f2*2) / sampRate
##         f3 = (f3*2) / sampRate
##         f4 = (f4*2) / sampRate
##         f = [0, f1, f2, ((f2+f3)/2), f3, f4, 1]
##         m = [0, 0.00003, 1, 1, 1, 0.00003, 0]
##     end

##     b = convert(Array{eltype(rec),1}, scisig[:firwin2](nTaps,f,m))
##     #b = ones(Float32,nTaps)
##     nChannels = size(rec)[1]
##     ## if channels == nothing
##     ##     channels = [1:nChannels]
##     ## end

##     nProcsToUse = min(length(workers()), size(rec, 1))
##     dRec = distribute(rec, procs=workers()[1:nProcsToUse], dist=[nProcsToUse,1])
##     @sync for w=1:length(workers()[1:nProcsToUse])
##         @spawnat workers()[w] doWorkFilterContinuousParallel!(dRec, b, channels=channels)
##     end
##     rec[:,:] = convert(Array, dRec)
##     close(dRec)
##     #darray_closeall()

##     #return rec
## end

    function doWorkFilterContinuousParallel!(rec::AbstractMatrix{T}, b::AbstractVector{S}; channels::Union{Q, AbstractVector{Q}}=collect(1:size(rec,1))) where {T<:Real, S<:Real, Q<:Int}
        #lIdx = collect(localindexes(rec)[1])
        lIdx = collect(DistributedArrays.localindices(rec)[1])
    nRows = length(lIdx)
    for i=1:nRows
        if in(lIdx[i], channels)
            if VERSION < v"0.5-"
                localpart(rec)[i,:] = fftconvolve(reshape(localpart(rec)[i,:], size(localpart(rec)[i,:], 2)), b, "same")
                localpart(rec)[i,:] = reverse(fftconvolve(reverse(reshape(localpart(rec)[i,:], size(localpart(rec)[i,:], 2)), dims=1), b, "same"), dims=1)
            else
                localpart(rec)[i,:] = fftconvolve(localpart(rec)[i,:], b, "same")
                localpart(rec)[i,:] = reverse(fftconvolve(reverse(localpart(rec)[i,:], dims=1), b, "same"), dims=1)
            end
        end
    end
end
    

"""
"""
function _centered(arr, newsize)
    # Return the center newsize portion of the array.
    currsize = size(arr)[1]
    startind = div((currsize - newsize), 2) + 1
    endind = startind + newsize -1 #check indexing is the same pyth julia?
    return arr[startind:endind]
end

"""
Convolve two 1-dimensional arrays using the FFT.

$(SIGNATURES)

##### Arguments

* `x::AbstractVector{T}`: First input.
* `y::AbstractVector{T}`: Second input. Should have the same number of dimensions as `x`; if sizes of `x` and `y` are not equal then `x` has to be the larger array.
* `mode::String`: A string indicating the size of the output:
    * "full": The output is the full discrete linear convolution of the inputs. (Default)
    * "valid": The output consists only of those elements that do not rely on the zero-padding.
    * "same": The output is the same size as `x`, centered with respect to the "full" output.

##### Returns

* `out::AbstractVector{T}`: An 1-dimensional array containing a subset of the discrete linear convolution of `x` with `y`.

##### Examples

```julia
    x = rand(1:10, 10)
    y = rand(1:10, 10)
    fftconvolve(x, y, "same")
```

"""
function fftconvolve(x::AbstractVector{T}, y::AbstractVector{R}, mode::String) where {T<:Real, R<:Real}
    s1 = size(x)[1]#check if array has two dim?
    s2 = size(y)[1]
    #println(typeof(x), typeof(y))
    convArray = conv(x,y)
    if mode == "full"
        return convArray
    elseif mode == "same"
        return _centered(convArray, s1)
    elseif mode == "valid"
        return _centered(convArray, abs(s1 - s2) + 1)
    end
end

"""
Find peaks in a waveform. This function computes the first derivative of the waveform to find
critical points. Then finds the peaks by selecting the critical points at which the second
derivative is negative.

$(SIGNATURES)

##### Arguments

* `y::Union{AbstractMatrix{Real},AbstractVector{Real}}`: Waveform of which peaks are sought. `y` is converted to
                                                         a vector before processing, so sensible inputs would be a vector,
                                                         a Nx1 matrix, or a 1xN matrix.
* `sampRate::Real`: Waveform sampling rate.
* `epochStart::Real`: the time, in seconds, at which the epoch starts relative to the starting time of the stimulus
                      (i.e. if the epoch starts before the start of the stimulus `epochStart` is negative).

##### Returns

* `peakPoints::AbstractVector{Real}`: The indexes of the peak points.
* `peakTimes::AbstractVector{Real}`: The times of the peak points.

##### Examples

```julia
    sampRate = 256; nTaps=64
    rec, evtTab = simulateRecording(nChans=1, dur=4, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "lowpass", nTaps, [2], channels=[1], transitionWidth=0.2)
    pkPnts, pkTimes = findPeaks(rec[1,:], sampRate)
    ## not run
    ## using PlotlyJS
    ## s1 = scatter(;x=collect(1:length(rec[1,:])), y=rec[1,:])
    ## s2 = scatter(;x=pkPnts, y=rec[1, pkPnts], mode="markers")
    ## plot([s1, s2])
```
"""
function findPeaks(y::Union{AbstractMatrix{T}, AbstractVector{T}}, sampRate::Real; epochStart::Real=0) where {T<:Real}
    peakPnts = (Int)[] #peak points
    crPnts = (Int)[] #critical points
    
    y = vec(y)
    dy = diff(y)
    ddy = diff(dy)
    #    crPnts = find(dy[1:end-1].*dy[2:end] .<= 0) +1
    sgnch = dy[1:end-1].*dy[2:end]
    for i=1:length(sgnch)
        if signbit(sgnch[i]) == true
            push!(crPnts, i+1)
        end
    end

    for n=1:length(crPnts)
        if ddy[crPnts[n]-1] .< 0
            push!(peakPnts, crPnts[n])
        end     
    end

    peakTimes = zeros(length(peakPnts))
    
    for i=1:length(peakPnts)
        peakTimes[i] = (peakPnts[i]-1)/sampRate
    end
    peakTimes = peakTimes .+ epochStart
    
    return peakPnts, peakTimes
end

"""
Find troughs in a waveform. This function computes the first derivative of the waveform to find
critical points. Then finds the troughs by selecting the critical points at which the second
derivative is positive.

$(SIGNATURES)

##### Arguments

* `y::Union{AbstractMatrix{Real},AbstractVector{Real}}`: Waveform of which troughs are sought. `y` is converted to
                                                         a vector before processing, so sensible inputs would be a vector,
                                                         a Nx1 matrix, or a 1xN matrix.
* `sampRate::Real`: Waveform sampling rate.
* `epochStart::Real`: the time, in seconds, at which the epoch starts relative to the starting time of the stimulus
                      (i.e. if the epoch starts before the start of the stimulus `epochStart` is negative).

##### Returns

* `troughsPoints::AbstractVector{Real}`: The indexes of the troughs points.
* `troughsTimes::AbstractVector{Real}`: The times of the troughs points.

##### Examples

```julia
    sampRate = 256; nTaps=64
    rec, evtTab = simulateRecording(nChans=1, dur=4, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "lowpass", nTaps, [2], channels=[1], transitionWidth=0.2)
    trPnts, trTimes = findTroughs(rec[1,:], sampRate)
    ## not run
    ## using PlotlyJS
    ## s1 = scatter(;x=collect(1:length(rec[1,:])), y=rec[1,:])
    ## s2 = scatter(;x=trPnts, y=rec[1, trPnts], mode="markers")
    ## plot([s1, s2])
```
"""
function findTroughs(y::Union{AbstractMatrix{T}, AbstractVector{T}}, sampRate::Real; epochStart::Real=0) where {T<:Real}
    troughPnts = (Int)[]
    crPnts = (Int)[]
    
    y = vec(y)
    dy = diff(y)
    ddy = diff(dy)
    #crPnts = find(dy[1:end-1].*dy[2:end] .<= 0) +1
    sgnch = dy[1:end-1].*dy[2:end]
    for i=1:length(sgnch)
        if signbit(sgnch[i]) == true
            push!(crPnts, i+1)
        end
    end
    for n=1:length(crPnts)
        if ddy[crPnts[n]-1] .> 0
            push!(troughPnts, crPnts[n])
        end     
    end

    troughTimes = zeros(length(troughPnts))
    
    for i=1:length(troughPnts)
        troughTimes[i] = (troughPnts[i]-1)/sampRate
    end
    troughTimes = troughTimes .+ epochStart

    return troughPnts, troughTimes
end

"""
Find inflection points in a waveform. This function finds inflection points by looking for sign changes in the second derivative. 
derivative is positive.

$(SIGNATURES)

##### Arguments

* `y::Union{AbstractMatrix{Real},AbstractVector{Real}}`: Waveform of which inflection points are sought. `y` is converted to
                                                         a vector before processing, so sensible inputs would be a vector,
                                                         a Nx1 matrix, or a 1xN matrix.
* `sampRate::Real`: Waveform sampling rate.
* `epochStart::Real`: the time, in seconds, at which the epoch starts relative to the starting time of the stimulus
                      (i.e. if the epoch starts before the start of the stimulus `epochStart` is negative).

##### Returns

* `inflPnts::AbstractVector{Real}`: The indexes of the inflection points.
* `inflTimes::AbstractVector{Real}`: The times of the inflection points.

##### Examples

```julia
    sampRate = 256; nTaps=64
    rec, evtTab = simulateRecording(nChans=1, dur=4, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "lowpass", nTaps, [2], channels=[1], transitionWidth=0.2)
    inflPnts, inflTimes = findInflections(rec[1,:], sampRate)
    ## not run
    ## using PlotlyJS
    ## s1 = scatter(;x=collect(1:length(rec[1,:])), y=rec[1,:])
    ## s2 = scatter(;x=inflPnts, y=rec[1,inflPnts], mode="markers")
    ## plot([s1, s2])
```
"""
function findInflections(y::Union{AbstractMatrix{T}, AbstractVector{T}}, sampRate::Real; epochStart::Real=0) where {T<:Real}
    inflPnts = (Int)[]
    y = vec(y)
    dy = diff(y)
    ddy = diff(dy)
    sgnch = ddy[1:end-1].*ddy[2:end]
    for i=1:length(sgnch)
        if signbit(sgnch[i]) == true
            push!(inflPnts, i+2)
        end
    end

    inflTimes = zeros(length(inflPnts))
    
    for i=1:length(inflPnts)
        inflTimes[i] = (inflPnts[i]-1)/sampRate
    end
    inflTimes = inflTimes .+ epochStart
    
    return inflPnts, inflTimes

end


"""
Find epochs with voltage values exceeding a given threshold.

$(SIGNATURES)

##### Arguments

* `rec::Dict{String, Array{Real, 3}}`: The segmented recording.
* `thresh::Union{Real, AbstractVector{Real}}`: The threshold value(s).
* `chans::AbstractVector{Int}`: The indexes of the channels on which to find artefacts.
* `chanlabels::AbstractVector{String}`: The labels of the channels on which to find artefacts.
* `chanList::AbstractVector{String}`: The names of all the channels.

##### Returns

* `segsToReject::Dict{String,Array{Int64,1}}`: dictionary containing the list of segments to reject for each condition.

##### Notes

If neither channel indexes (`chans`) nor channel labels (`chanLabels`)
for the channels on which to check artefacts are provided, then artefacts
will be checked for on all channels.
If channel labels (`chanLabels`) are given for the channels on which to
check for artefacts, then an ordered list containing the names of all available
channels (`chanList`) must be provided as well.

`thresh` should be a list of threshold values, one for each channel to check.
If `thresh` contains only one value, it is assumed that this is the desired
threshold for all of the channels to check. If `thresh` contains more than one
value, its length must match the number of channels to check.

##### Examples

```julia
    epochDur=0.5; preDur=0.15; events=[1,2]; sampRate=256;
    rec, evtTab = simulateRecording(dur=120, epochDur=epochDur, preDur=preDur, events=events)
    segs, nRaw = segment(rec, evtTab, -preDur, epochDur, sampRate)
    # on all channels
    badSegs = findArtefactThresh(segs, 65)
    # on channels 1 and 2
    badSegs = findArtefactThresh(segs, 65, [1,2])
    # on channels FP1 and F4 #not run
    #findArtefactThresh(segs, 20, ["Fp1", "F4"], chanLabels)
```

"""
function findArtefactThresh(rec::Dict{String, Array{T, 3}}, thresh::Union{P, AbstractVector{P}},
                            channels::Union{Q, AbstractVector{Q}}=collect(1:size(rec[collect(keys(rec))[1]], 1))) where {T<:Real, P<:Real, Q<:Integer}

    eventList = collect(keys(rec))
    ## if chans != nothing
    ##     channels = chans
    ## elseif chanLabels != nothing
    ##     channels = (Int)[]
    ##     for i=1:length(chanLabels)
    ##         push!(channels, find(chanList .== chanLabels[i])[1])
    ##     end
    ## else
    ##     channels = [1:size(rec[eventList[1]])[1]] #assume n channels the same for all dict entries
    ## end

    if length(channels) != length(thresh)
        if length(thresh) == 1
            thresh = [thresh[1] for i=1:length(channels)]
        else
            error("The number of thresholds must be equal to 1 or to the number of channels \n")
            return
        end
    end

    segsToReject = Dict{String,Array{Int,1}}()#(String => Array{Int,1})[]
    for i=1:length(eventList)
        segsToReject[eventList[i]] = []
        for j=1:size(rec[eventList[i]])[3]
            for k=1:length(channels) #size(rec[eventList[i]])[1]
                thisChan = channels[k]
                if (maximum(rec[eventList[i]][thisChan,:,j]) > thresh[k] || minimum(rec[eventList[i]][thisChan,:,j]) < -thresh[k]) == true
                    #println(thresh[k])
                    #println(min(rec[eventList[i]][k,:,j]))
                    segsToReject[eventList[i]] = vcat(segsToReject[eventList[i]], j)
                end
            end
        end
    end



    for i=1:length(eventList)
        #segment may be flagged by detection in more than one channel
        segsToReject[eventList[i]] = unique(segsToReject[eventList[i]])
    end
    return segsToReject

end

function findArtefactThresh(rec::Dict{String, Array{T, 3}}, thresh::Union{P, AbstractVector{P}}, chanLabels::AbstractVector{S}, chanList::AbstractVector{R}) where {T<:Real, P<:Real, R<:String, S<:String}
    channels = (Int)[]
    for i=1:length(chanLabels)
        push!(channels, find(chanList .== chanLabels[i])[1])
    end
    segsToReject = findArtefactThresh(rec, thresh, channels)
    return segsToReject
end

"""
Find the time point at which a waveform reaches a maximum or a minimum.

$(SIGNATURES)

##### Arguments:

* `wave::Union{AbstractVector{Real}, AbstractMatrix{Real}}`: the waveform for which the extremum should be found.
* `searchStart::Real`: the starting time point, in seconds, of the window in which to search for the extremum.
* `searchStop::Real`: the stopping time point, in seconds, of the window in which to search for the extremum.
* `extremumSign::String`: whether the sought extremum is of `positive` or `negative` sign.
* `epochStart::Real`: the time, in seconds, at which the epoch starts.
* `samprate::Real`: the sampling rate of the signal.

##### Returns

* `extremumPnt::Real`: the sample number at which the extremum occurs.
* `extremumTime::Real`: the time, in seconds, at which the extremum occurs.

##### Examples

```julia
    ## not run
    ## P2SearchStart = 0.150
    ## P2SearchStop = 0.250
    ## epochStart = -0.150
    ## sampRate = 8192
    ## sampPnt, timePnt = findExtremum(wave1, P2SearchStart, P2SearchStop, "positive", epochStart, sampRate)

    # contrived example
    #using Winston
    sampRate = 256
    dur = 0.6
    epochStart = -0.15
    P2SearchStart = 0.150
    P2SearchStop = 0.250
    nSamp = round(Int, dur*sampRate)
    freq = 2
    phase = 1/4*(-pi)
    tArr = collect(0:nSamp-1)/sampRate + epochStart
    wave1 = sin(2*pi*freq*tArr+phase)
    wave1[1:round(Int, abs(epochStart)*sampRate)] = 0
    sampPnt, timePnt = findExtremum(wave1, P2SearchStart, P2SearchStop, "positive", epochStart, sampRate)
    #p = plot(tArr, wave1)
    #l1 = LineX(timePnt, color="red")
    #add(p, l1)
    #display(p)
```

"""
function findExtremum(wave::Union{AbstractVector{T}, AbstractMatrix{T}}, searchStart::Real, searchStop::Real, extremumSign::String, epochStart::Real, sampRate::Real) where {T<:Real}

    if ndims(wave) > 1
        if in(1, size(wave)) == false
            error("Only 1-dimensional arrays, or 1xN dimensional arrays are allowed")
        end
        wave = vec(wave)
    end

    searchStartPnt = round(Int, (searchStart - epochStart)*sampRate)
    searchStopPnt = round(Int, (searchStop - epochStart)*sampRate)

    searchWin = wave[searchStartPnt:searchStopPnt]
    if extremumSign == "positive"
        extremumPntRel = find(searchWin .== maximum(searchWin))
    elseif extremumSign == "negative"
        extremumPntRel = find(searchWin .== minimum(searchWin))
    end

    if length(extremumPntRel) > 1
        println("Warning: more than one extrema detected with the same amplitude. Selecting the first one...")
    end

    tArr = collect(0:length(wave)-1)/sampRate + epochStart
    extremumPnt = extremumPntRel[1] + searchStartPnt-1
    #extremumTime = (extremumPnt / sampRate) + epochStart
    extremumTime = tArr[extremumPnt]

    return extremumPnt, extremumTime
end

"""

$(SIGNATURES)

##### References

* Elberling, C., & Don, M. (1984). Quality estimation of averaged auditory brainstem responses. Scandinavian Audiology, 13(3), 187–197.
* Ozdamar, O., & Delgado, R. E. (1996). Measurement of signal and noise characteristics in ongoing auditory brainstem response averaging. Annals of Biomedical Engineering, 24(6), 702–715.
"""
function FMP(segs::Array{T,3}; at=-1) where {T<:Real}
    if at == -1
        at = size(segs)[3]
    end
    nChans = size(segs)[1]
    nSamples = size(segs)[2]
    #nSweeps = size(segs)[3]
    n = length(at)
    FMPNum=zeros(nChans, n)
    FMPDen=zeros(nChans, n)

    for i=1:n
        thisIdx = at[i]
        FMPNum[:,i] = var(squeeze(mean(segs[:,:,1:thisIdx],3),3),2) #Variance across time of the averaged sweeps. The variance is calculated for each channel and.
        tmpFMPDen = zeros(nChans, nSamples)
        for j=1:nSamples
            tmpFMPDen[:,j] = squeeze(var(segs[:,j,1:thisIdx], 3), (2,3))/thisIdx #Variance across sweeps for each sample. The variance is calculated for each channel.
        end
        FMPDen[:,i] = mean(tmpFMPDen,2) #Average the variance across samples.
    end
    FMPVals = FMPNum ./ FMPDen
    return FMPVals, FMPNum, FMPDen

end

function FMP(segs::Dict{String,Array{T,3}}; at=-1) where {T<:Real}

    eventList = collect(keys(segs))
    FMPVals = Dict{String,Array{Float64,2}}()
    FMPNum = Dict{String,Array{Float64,2}}()
    FMPDen = Dict{String,Array{Float64,2}}()

    for i=1:length(eventList)
        FMPVals[eventList[i]], FMPNum[eventList[i]], FMPDen[eventList[i]] = FMP(segs[eventList[i]], at=at)
    end

    return FMPVals, FMPNum, FMPDen
end


"""
$(SIGNATURES)
"""
function FMPIterativeWeighted(sweeps::Array{T,3}; at=-1, noiseEstimate::String="global") where {T<:Real}
    if at == -1
        at = size(sweeps)[3]
    end
    n=length(at)
    nChans = size(sweeps)[1]
    nSamp = size(sweeps)[2]
    nSweeps = size(sweeps)[3]

    weighted_ave = zeros(nChans, nSamp)
    weighted_ave_num = zeros(nChans, nSamp)
    weighted_ave_den = zeros(nChans, nSamp)

    FMPNum=zeros(nChans, n)
    FMPDen=zeros(nChans, n)

    if noiseEstimate == "global"
        sweepWeights = zeros(nSweeps)
        for i=1:nSweeps
            if i==1
                est_sig = zeros(nChans, nSamp)
            else
                est_sig = weighted_ave
            end
            noisePow = mean((sweeps[:,:,i] .- est_sig).^2)
            weighted_ave_num = weighted_ave_num .+ sweeps[:,:,i].* (1/noisePow)

            sweepWeights[i] = 1 ./ noisePow
            weighted_ave_den = weighted_ave_den .+ (1/noisePow)
            weighted_ave = weighted_ave_num ./ weighted_ave_den

            if in(i, at)
                thisIdx = find(at.==i)[1]
                FMPNum[:,thisIdx] = var(weighted_ave,2) #variance of the averaged sweeps
                resNoise = zeros(nChans, nSamp) #residual noise
                for ch=1:nChans
                    for s=1:nSamp
                        resNoise[ch,s] = sum((squeeze(sweeps[ch,s,1:i].-weighted_ave[ch,s], (1,2)).^2).*sweepWeights[1:i]) ./ ((i-1)*sum(sweepWeights[1:i]))
                    end
                end
                FMPDen[:,thisIdx] = mean(resNoise, 2)
            end

        end

    elseif noiseEstimate == "byChannel"
        sweepWeights = zeros(nChans, nSweeps)
        for i=1:nSweeps
            if i==1
                est_sig = zeros(nChans, nSamp)
            else
                est_sig = weighted_ave
            end
            for c=1:nChans
                noisePow = mean((sweeps[c,:,i] .- est_sig[c,:]).^2)
                sweepWeights[c, i] = 1 ./noisePow
                weighted_ave_num[c,:] = weighted_ave_num[c,:] .+ sweeps[c,:,i].* (1/noisePow)
                weighted_ave_den[c,:] = weighted_ave_den[c,:] .+ (1/noisePow)
            end
            weighted_ave = weighted_ave_num ./ weighted_ave_den

            if in(i, at)
                thisIdx = find(at.==i)[1]
                FMPNum[:,thisIdx] = var(weighted_ave,2) #variance of the averaged sweeps
                resNoise = zeros(nChans, nSamp) #residual noise
                for c=1:nChans
                    for s=1:nSamp
                        resNoise[c,s] = sum((squeeze(sweeps[c,s,1:i].-weighted_ave[c,s], (1,2)).^2).*squeeze(sweepWeights[c,1:i],1)) ./ ((i-1)*sum(sweepWeights[c,1:i]))
                    end
                end
                FMPDen[:,thisIdx] = mean(resNoise, 2)
            end

        end
    end
    FMPVals = FMPNum ./ FMPDen
    return FMPVals, FMPNum, FMPDen
end


function FMPIterativeWeighted(segs::Dict{String,Array{T,3}}; at=-1, noiseEstimate::String="global") where {T<:Real}

    eventList = collect(keys(segs))
    FMPVals = Dict{String,Array{Float64,2}}()
    FMPNum = Dict{String,Array{Float64,2}}()
    FMPDen = Dict{String,Array{Float64,2}}()

    for i=1:length(eventList)
        FMPVals[eventList[i]], FMPNum[eventList[i]], FMPDen[eventList[i]] = FMPIterativeWeighted(segs[eventList[i]], at=at, noiseEstimate=noiseEstimate)
    end

    return FMPVals, FMPNum, FMPDen
end



"""
Compute the autocorrelation function of a 1-dimensional signal.

$(SIGNATURES)

##### Arguments:

* `sig::Union{AbstractVector{Real}, AbstractMatrix{Real}}`: the signal for which the autocorrelation should be computed.
* `samprate::Real`: the sampling rate of the signal.
* `maxLag::Real`: the maximum lag (1/f) for which the autocorrelation function should be computed.
* `normalize::Bool`: whether the autocorrelation should be scaled between [-1, 1].
* `window::Function`: The type of window to apply to the signal before computing its ACF (see DSP.jl).
                      Choose `rect` if you don't want to apply any window.

##### Returns

* `acf::Array{Real,1}`: the autocorrelation function.
* `lags::Array{Real,1}`: the time lags for which the autocorrelation function was computed.

##### Examples

```julia
    using DSP
    sampRate = 48000
    dur = 1
    nSamp = round(Int, sampRate*dur)
    tArr = collect(0:nSamp-1)/sampRate
    freq = 440
    sig = sin(2*pi*freq*tArr)
    maxLag = 1/200
    acf, lags = getACF(sig, sampRate, maxLag, normalize=true, window=hamming)
```

"""
function getACF(sig::Union{AbstractVector{T}, AbstractMatrix{T}}, sampRate::Real, maxLag::Real; normalize::Bool=true, window::Function=rect) where {T<:Real}
##     """
## n = length(sig)
## acfArray = zeros(n*2)
## acfArray[1:n] = sig
## out = zeros(n)

## maxLagPnt = round(Int, maxLag*sampRate)
## if maxLagPnt > n
## maxLagPnt = n
## end

## for i = 1:maxLagPnt
## out[i] = sum(acfArray[1:n] .* acfArray[i:(n+i-1)])
## end

## lags = [1:maxLagPnt] ./ sampRate

## if normalize == true
## out = out ./ maximum(out)
## end
## return out, lags
## """

    if ndims(sig) > 1
        if in(1, size(sig)) == false
            error("Only 1-dimensional arrays, or 1xN dimensional arrays are allowed")
        end
        sig = vec(sig)
    end

    n = length(sig)
    w = window(n)
    sig = sig.*w

    maxLagPnt = round(Int, maxLag*sampRate)
    if maxLagPnt > n
        maxLagPnt = n
    end
    out = xcorr(sig, sig)[n:n+maxLagPnt-1]

    lags = collect(1:maxLagPnt) ./ sampRate

    if normalize == true
        out = out ./ maximum(out)
    end
    return out, lags

end

function getACF2(sig::Union{AbstractVector{T}, AbstractMatrix{T}}, sampRate::Real, maxLag::Real; normalize::Bool=true, window::Function=rect) where {T<:Real}

    if ndims(sig) > 1
        if in(1, size(sig)) == false
            error("Only 1-dimensional arrays, or 1xN dimensional arrays are allowed")
        end
        sig = vec(sig)
    end

    n = length(sig)
    w = window(n)
    sig = sig.*w

    acfArray = zeros(n*2)
    acfArray[1:n] = sig
    out = zeros(n)

    maxLagPnt = round(Int, maxLag*sampRate)
    if maxLagPnt > n
        maxLagPnt = n
    end

    for i = 1:maxLagPnt
        out[i] = sum(acfArray[1:n] .* acfArray[i:(n+i-1)])
    end

    out = out[1:maxLagPnt]
    lags = [1:maxLagPnt] ./ sampRate

    if normalize == true
        out = out ./ maximum(out)
    end
    return out, lags
end

"""
Compute the autocorrelogram of a 1-dimensional array.

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractVector{Real}, AbstractMatrix{Real}}` The signal of which the autocorrelogram should be computed.
* `sampRate::Real`: The sampling rate of the signal.
* `winLength::Real`: The length of the window over which to take the ACF, in seconds.
* `overlap::Real`: The percent of overlap between successive windows (useful for smoothing the autocorrelogram).
* `maxLag::Real`: the maximum lag (1/f) for which the autocorrelation function should be computed.
* `normalize::Bool`: whether the autocorrelation should be scaled between [-1, 1].
* `window::Function`: The type of window to apply to the signal before computing its ACF (see DSP.jl).
                      Choose `rect` if you don't want to apply any window.

##### Returns

* `acfMatrix::Array{Real,2}`: the autocorrelogram.
* `lags::Array{Real,1}`: the ACF lags.
* `timeArray::Array{Real,1}`: The time axis.

##### Examples

```julia
    sig = rand(512)
    acg, lags, t = getAutocorrelogram(sig, 256, 0.02, 30, 0.01)
```

"""
function getAutocorrelogram(sig::Union{AbstractVector{T}, AbstractMatrix{T}}, sampRate::Real, winLength::Real, overlap::Real, maxLag::Real; normalize::Bool=true, window::Function=rect) where {T<:Real}
    winLengthPnt = floor(Int, winLength * sampRate)
    stepSize = winLengthPnt - round(Int, winLengthPnt * overlap / 100)
    ind = collect(1:stepSize:length(sig) - winLengthPnt)
    n = length(ind)
    acf, lags = getACF(sig[ind[1]:ind[1]+winLengthPnt], sampRate, maxLag, normalize=normalize, window=window)

    acfMatrix = zeros(length(acf), n)
    acfMatrix[:,1] = acf
    for i=2:n
        acf, lags = getACF(sig[ind[i]:ind[i]+winLengthPnt], sampRate, maxLag, normalize=normalize, window=window)
        acfMatrix[:,i] = acf
    end

    #timeInd = arange(0, len(sig), stepSize)
    #timeArray = 1 ./ sampRate * (timeInd)
    timeArray3 = linspace(0, (ind[end]+winLengthPnt-1)/sampRate, n+1)
    return acfMatrix, lags, timeArray3
end

function getAutocorrelogram2(sig::Union{AbstractVector{T}, AbstractMatrix{T}}, sampRate::Real, winLength::Real, overlap::Real, maxLag::Real; normalize::Bool=true, window::Function=rect) where {T<:Real}
    winLengthPnt = floor(Int, winLength * sampRate)
    stepSize = winLengthPnt - round(Int, winLengthPnt * overlap / 100)
    ind = collect(1:stepSize:length(sig) - winLengthPnt)
    indStop = ind + winLengthPnt
    #indStop[end] = length(sig)
    n = length(ind)
    acf, lags = getACF2(sig[ind[1]:indStop[1]], sampRate, maxLag, normalize=normalize, window=window)
    nACF = length(acf)
    acfMatrix = zeros(length(acf), n)
    acfMatrix[:,1] = acf
    for i=2:n
        acf, lags = getACF2(sig[ind[i]:indStop[i]], sampRate, maxLag, normalize=normalize, window=window)
        acfMatrix[:,i] = acf[1:nACF]
    end

    #timeInd = arange(0, len(sig), stepSize)
    #timeArray = 1 ./ sampRate * (timeInd)
    timeArray3 = linspace(0, (ind[end]+winLengthPnt-1)/sampRate, n+1)
    return acfMatrix, lags, timeArray3
end

"""
Compute the signal-to-noise ratio at a given frequency in the power spectrum of a recording.

$(SIGNATURES)

##### Arguments

* `spec::AbstractVector{Real}` The power spectrum of the recording.
* `freqArr::AbstractVector{Real}`: the FFT frequency array.
* `sigFreq::Real`: the signal frequency of interest.
* `nSideComp::Integer`: the number of components adjacent to the signal used to estimate the noise.
* `nExclude::Integer`: the number of components closest to the signal to exclude from the noise estimate.

##### Returns

* `snr::Real`: The signal-to-noise ratio at the target frequency.

##### Examples

```julia
    sig = rand(512)
    p, f = getSpectrum(sig, 512)
    snr = getSNR(p, f, 140, 10, 1)
```

"""
function getSNR(spec::AbstractVector{T}, freqArr::AbstractVector{R}, sigFreq::Real, nSideComp::Integer, nExclude::Integer) where {T<:Real, R<:Real}

    sigIdx = find(abs.(freqArr .- sigFreq) .== minimum(abs.(freqArr .- sigFreq)))[1]
    sigMag = spec[sigIdx]
    loNoiseMag = spec[sigIdx-nExclude-1-nSideComp+1:sigIdx-nExclude-1]
    hiNoiseMag = spec[sigIdx+nExclude+1:sigIdx+nExclude+1+nSideComp-1]
    noiseMag = mean([loNoiseMag; hiNoiseMag])
    snr = 10*log10(sigMag ./ noiseMag)
    return snr
end

"""
Compute the signal-to-noise ratio at a given frequency in the power spectrum of a recording.
This function is the same as `getSNR`, but it additionaly returns the signal and noise magnitudes separately.

$(SIGNATURES)

##### Arguments

* `spec::AbstractVector{Real}` The power spectrum of the recording.
* `freqArr::AbstractVector{Real}`: the FFT frequency array.
* `sigFreq::Real`: the signal frequency of interest.
* `nSideComp::Integer`: the number of components adjacent to the signal used to estimate the noise.
* `nExclude::Integer`: the number of components closest to the signal to exclude from the noise estimate.

##### Returns

* `snr::Real`: The signal-to-noise ratio at the target frequency.
* `sigMag::Real`: The signal magnitude.
* `noiseMag::Real`: The noise magnitude.

##### Examples

```julia
    sig = rand(512)
    p, f = getSpectrum(sig, 512)
    snr, sigMag, noiseMag = getSNR2(p, f, 140, 10, 1)
```

"""
function getSNR2(spec::AbstractVector{T}, freqArr::AbstractVector{R}, sigFreq::Real, nSideComp::Integer, nExclude::Integer) where {T<:Real, R<:Real}

    sigIdx = findall(abs.(freqArr .- sigFreq) .== minimum(abs.(freqArr .- sigFreq)))[1]
    sigMag = spec[sigIdx]
    loNoiseMag = spec[sigIdx-nExclude-1-nSideComp+1:sigIdx-nExclude-1]
    hiNoiseMag = spec[sigIdx+nExclude+1:sigIdx+nExclude+1+nSideComp-1]
    noiseMag = mean([loNoiseMag; hiNoiseMag])
    snr = 10*log10(sigMag ./ noiseMag)
    return snr, sigMag, noiseMag
end

"""
Compute the spectrogram of a 1-dimensional array.

$(SIGNATURES)

##### Arguments
* `sig::Union{AbstractVector{Real}, AbstractMatrix{Real}}` The signal of which the spectrum should be computed.
* `sampRate::Real`: The sampling rate of the signal.
* `winLength::Real`: The length of the window over which to take the FFTs, in seconds.
* `overlap::Real`: The percent of overlap between successive windows (useful for smoothing the spectrogram).
* `window::Function`: The type of window to apply to the signal before computing its FFT (see DSP.jl).
                      Choose `rect` if you don't want to apply any window.
* `powerOfTwo::Bool`: If `true` `sig` will be padded with zeros (if necessary) so that its length is a power of two.

##### Returns

* `powerMatrix::Array{Real, 2}`: the power spectrum for each time window.
* `freqArray::Array{Real, 1}`: The frequency axis.
* `timeArray::Array{Real, 1}`: The time axis.

##### Notes

If the signal length is not a multiple of the window length it is trucated.

##### Examples

```julia
    sig = rand(512)
    spec, f, t = getSpectrogram(sig, 256, 0.02, 30)
```

"""
function getSpectrogram(sig::Union{AbstractVector{T}, AbstractMatrix{T}}, sampRate::Real, winLength::Real, overlap::Real; window::Function=rect, powerOfTwo::Bool=false) where {T<:Real}
    winLengthPnt = floor(Int, winLength * sampRate)

    stepSize = winLengthPnt - round(Int, winLengthPnt * overlap / 100)
    ind = collect(1:stepSize:length(sig) - winLengthPnt)
    n = length(ind)
    p, freqArray = getSpectrum(sig[ind[1]:ind[1]+winLengthPnt], sampRate, window=window, powerOfTwo=powerOfTwo)

    powerMatrix = zeros(length(freqArray), n)
    powerMatrix[:,1] = p
    for i=2:n
        p, freqArray = getSpectrum(sig[ind[i]:ind[i]+winLengthPnt], sampRate, window=window, powerOfTwo=powerOfTwo)
        powerMatrix[:,i] = p
    end
    timeArray = linspace(0, (ind[end]+winLengthPnt-1)/sampRate, n+1)

    return powerMatrix, freqArray, timeArray
end

"""
Compute the power spectrum of a 1-dimensional array.

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractVector{Real}, AbstractMatrix{Real}}`: The signal of which the spectrum should be computed.
* `sampRate::Real`: The sampling rate of the signal.
* `window::Function`: The type of window to apply to the signal before computing its FFT (see DSP.jl).
                      Choose `rect` if you don't want to apply any window.
* `powerOfTwo::Bool`: If `true` `sig` will be padded with zeros (if necessary) so that its length is a power of two.

##### Returns

* `p::Array{Real,1}`: the power spectrum of the signal.
* `freqArray::Array{Real,1}`: The FFT frequencies.

##### Examples

```julia
    sig = rand(512)
    p, f = getSpectrum(sig, 256)
```

"""
function getSpectrum(sig::Union{AbstractVector{T}, AbstractMatrix{T}}, sampRate::Integer; window::Function=rect, powerOfTwo::Bool=false) where {T<:Real}
    if ndims(sig) > 1
        if in(1, size(sig)) == false
            error("Only 1-dimensional arrays, or 1xN dimensional arrays are allowed")
        end
        sig = vec(sig)
    end
    n = length(sig)
    w = window(n)
    sig = sig.*w
    if powerOfTwo == true
        nfft = nextpow2(n) #2^nextPowTwo(n)
        if nfft != n #zero-pad
            sig = vcat(sig, zeros(eltype(sig), nfft-n))
        end
    else
        nfft = n
    end

    p = fft(sig)#, nfft) # take the fourier transform

    nUniquePts = ceil(Int, (nfft+1)/2)
    p = p[1:nUniquePts]
    p = abs.(p)
    p = p ./ n  # scale by the number of points so that
    # the magnitude does not depend on the length
    # of the signal or on its sampling frequency
    p = p.^2  # square it to get the power

    # multiply by two (see technical document for details)
    # odd nfft excludes Nyquist point
    if nfft % 2 > 0 # we"ve got odd number of points fft
         p[2:end] = p[2:end] * 2
    else
        p[2:(end-1)] = p[2:(end-1)] * 2 # we"ve got even number of points fft
    end

    freqArray = collect(0:(nUniquePts-1)) * (sampRate / nfft)

    return p, freqArray
end

"""
Compute the phase spectrum of a 1-dimensional array.

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractVector{Real}, AbstractMatrix{Real}}`: The signal of which the phase spectrum should be computed.
* `sampRate::Real`: The sampling rate of the signal.
* `window::Function`: The type of window to apply to the signal before computing its FFT (see DSP.jl).
                      Choose `rect` if you don't want to apply any window.
* `powerOfTwo::Bool`: If `true` `sig` will be padded with zeros (if necessary) so that its length is a power of two.

##### Returns

* `p::Array{Real,1}`: the phase spectrum of the signal.
* `freqArray::Array{Real,1}`: The FFT frequencies.

##### Examples

```julia
    sig = rand(512)
    p, f = getPhaseSpectrum(sig, 256)
```

"""
function getPhaseSpectrum(sig::Union{AbstractVector{T}, AbstractMatrix{T}}, sampRate::Real; window::Function=rect, powerOfTwo::Bool=false) where {T<:Real}
    if ndims(sig) > 1
        if in(1, size(sig)) == false
            error("Only 1-dimensional arrays, or 1xN dimensional arrays are allowed")
        end
        sig = vec(sig)
    end
    n = length(sig)
    w = window(n)
    sig = sig.*w
    if powerOfTwo == true
        nfft = nextpow2(n) #2^nextPowTwo(n)
        if nfft != n #zero-pad
            sig = vcat(sig, zeros(eltype(sig), nfft-n))
        end
    else
        nfft = n
    end

    p = fft(sig)#, nfft) # take the fourier transform

    nUniquePts = ceil(Int, (nfft+1)/2)
    p = p[1:nUniquePts]

    p = angle(p)
    freqArray = collect(0:(nUniquePts-1)) * (sampRate / nfft)

    return p, freqArray
end

"""
Compute an iterative weighted average from a matrix of segmented EEG recordings.

$(SIGNATURES)

##### Arguments

* `sweeps`::Array{T,3}`: The matrix containing the segmented recording from which the
average should be computed.
* `noiseEstimate::String`: if `global` the estimate of the noise used to weight individual
segments is derived from all the channels. If `byChannel` an noise estimate is derived for each
channel and segments are weighted differently for each channel depending on the channel noise
estimate.
* `noiseWinStart`::Real: Time in seconds at which the noise estimate should start relative to the start of the epoch.
* `noiseWinStop`::Real: Time in seconds at which the noise estimate should stop relative to the start of the epoch.
* `preDur::Real`: Duration of recording before the experimental event, in seconds.
* `sampRate::Real`: The samplig rate of the EEG recording.

##### Returns

* `weighted_ave`: The weighted average of the segments.

##### Examples

nChans = 3; nSamp=1000; nSweeps=500
sweeps = rand(nChans, nSamp, nSweeps)
iterativeWeightedAverage(sweeps, noiseEstimate="global")
iterativeWeightedAverage(sweeps, noiseEstimate="byChannel")
iterativeWeightedAverage(sweeps, 0, 0.2, 0.1, 1000, noiseEstimate="global")

##### References

* Riedel, H., Granzow, M., & Kollmeier, B. (2001). Single-sweep-based methods to improve the quality of auditory brain stem responses Part II: Averaging methods. Z Audiol, 40(2), 62–85.

"""
function iterativeWeightedAverage(sweeps::Array{T,3}; noiseEstimate::String="global") where {T<:Real}

    if in(noiseEstimate, ["global", "byChannel"]) == false
        error("`noiseEstimate` must be either `global`, or `byChannel`")
    end

    nSweeps = size(sweeps)[3]
    nSamp = size(sweeps)[2]
    nChans = size(sweeps)[1]

    weighted_ave = zeros(nChans, nSamp)
    weighted_ave_num = zeros(nChans, nSamp)
    weighted_ave_den = zeros(nChans, nSamp)

    if noiseEstimate == "global"
        for i=1:nSweeps
            if i==1
                est_sig = zeros(nChans, nSamp)
            else
                est_sig = weighted_ave
            end
            noisePow = mean((sweeps[:,:,i] .- est_sig).^2)
            weighted_ave_num = weighted_ave_num .+ sweeps[:,:,i].* (1/noisePow)
            weighted_ave_den = weighted_ave_den .+ (1/noisePow)
            weighted_ave = weighted_ave_num ./ weighted_ave_den
        end
    elseif noiseEstimate == "byChannel"
        for i=1:nSweeps
            if i==1
                est_sig = zeros(nChans, nSamp)
            else
                est_sig = weighted_ave
            end
            for c=1:nChans
                noisePow = mean((sweeps[c,:,i] .- est_sig[c,:]).^2)
                weighted_ave_num[c,:] = weighted_ave_num[c,:] .+ sweeps[c,:,i].* (1/noisePow)
                weighted_ave_den[c,:] = weighted_ave_den[c,:] .+ (1/noisePow)
            end
            weighted_ave = weighted_ave_num ./ weighted_ave_den
        end
    end

    return weighted_ave
end

function iterativeWeightedAverage(sweeps::Array{T,3}, noiseWinStart::Real, noiseWinStop::Real, preDur::Real, sampRate::Real; noiseEstimate::String="global") where {T<:Real}
    #noiseWinStart - time in seconds at which the noise estimate should start relative to the start of the epoch
    #noiseWinStop - time in seconds at which the noise estimate should stop relative to the start of the epoch

    if in(noiseEstimate, ["global", "byChannel"]) == false
        error("`noiseEstimate` must be either `global`, or `byChannel`")
    end

    nSweeps = size(sweeps)[3]
    nSamp = size(sweeps)[2]
    nChans = size(sweeps)[1]

    noiseWinStartSample = round(Int, (noiseWinStart+preDur)*sampRate)+1
    noiseWinStopSample = round(Int, (noiseWinStop+preDur)*sampRate)
    if noiseWinStartSample < 1
        noiseWinStartSample = 1
    end
    if noiseWinStopSample > nSamp
        noiseWinStopSample = nSamp
    end


    weighted_ave = zeros(nChans, nSamp)
    weighted_ave_num = zeros(nChans, nSamp)
    weighted_ave_den = zeros(nChans, nSamp)

    if noiseEstimate == "global"
        for i=1:nSweeps
            if i==1
                est_sig = zeros(nChans, nSamp)
            else
                est_sig = weighted_ave
            end
            noisePow = mean((sweeps[:,noiseWinStartSample:noiseWinStopSample,i] .- est_sig[:,noiseWinStartSample:noiseWinStopSample]).^2)
            weighted_ave_num = weighted_ave_num .+ sweeps[:,:,i].* (1/noisePow)
            weighted_ave_den = weighted_ave_den .+ (1/noisePow)
            weighted_ave = weighted_ave_num ./ weighted_ave_den
        end
    elseif noiseEstimate == "byChannel"
        for i=1:nSweeps
            if i==1
                est_sig = zeros(nChans, nSamp)
            else
                est_sig = weighted_ave
            end
            for c=1:nChans
                noisePow = mean((sweeps[c,noiseWinStartSample:noiseWinStopSample,i] .- est_sig[c,noiseWinStartSample:noiseWinStopSample]).^2)
                weighted_ave_num[c,:] = weighted_ave_num[c,:] .+ sweeps[c,:,i].* (1/noisePow)
                weighted_ave_den[c,:] = weighted_ave_den[c,:] .+ (1/noisePow)
            end
            weighted_ave = weighted_ave_num ./ weighted_ave_den
        end
    end

    return weighted_ave
end

"""
Compute the mean amplitude of an ERP waveform in a time window centered on a given point.

$(SIGNATURES)

#### Arguments:

* `wave::Union{AbstractVector{Real}, AbstractMatrix{Real}}`: the waveform for which the mean amplitude should be computed.
* `center::Real`: the center of the window in which the mean amplitude should be computed. `center` can be specified either
                  in terms of sample point number, or in terms of time in seconds. If center is specified in terms of time
                  in seconds, the time at which the epoch starts, should be passed as an additional argument to the function.
* `centerType::String`: whether the `center` is specified in terms of sample point number `point`, or interms of time in seconds `time`.
* `winLength::String`: the length of the window, in seconds, over which to compute the meam amplitude.
* `samprate::Real`: the sampling rate of the signal.
* `epochStart::Real`: the time, in seconds, at which the epoch starts.

#### Returns

* `meanAmp::Real`: the mean amplitude of the waveform in the time window centered at the specified point.

### Examples

    ```julia
    ## not run
    ## P2WinLength = 0.050
    ## centerPoint = 2512
    ## sampRate = 8192
    ## meanAmp =  meanERPAmplitude(wave, centerPoint, "point", P2WinLength, sampRate)

    ## epochStart = -0.150
    ## centerPointTm = 0.156
    ## meanAmp =  meanERPAmplitude(wave, centerPointTm, "time", P2WinLength, sampRate)

    # contrived example
    #using Winston
    sampRate = 256
    dur = 0.6
    epochStart = -0.15
    P2SearchStart = 0.150
    P2SearchStop = 0.250
    nSamp = round(Int, dur*sampRate)
    freq = 2
    phase = 1/4*(-pi)
    tArr = collect(0:nSamp-1)/sampRate + epochStart
    wave1 = sin(2*pi*freq*tArr+phase)
    wave1[1:round(Int, abs(epochStart)*sampRate)] = 0
    sampPnt, timePnt = findExtremum(wave1, P2SearchStart, P2SearchStop, "positive", epochStart, sampRate)
    #p = plot(tArr, wave1)
    #l1 = LineX(timePnt, color="red")
    #add(p, l1)
    #display(p)
    meanAmp =  meanERPAmplitude(wave1, sampPnt, "point", 0.05, sampRate)
    ```

"""
function meanERPAmplitude(wave::Union{AbstractVector{T}, AbstractMatrix{T}}, center::Real, centerType::String, winLength::Real, sampRate::Real) where {T<:Real}

    if ndims(wave) > 1
        if in(1, size(wave)) == false
            error("Only 1-dimensional arrays, or 1xN dimensional arrays are allowed")
        end
        wave = vec(wave)
    end

    if centerType == "point"
        centerPnt = center
    elseif centerType == "time"
        error("If `centerType` is `time` you need to specify `epochStart`") #centerPnt = (center-epochStart)*sampRate
    end

    startPnt = centerPnt - round(Int, winLength/2*sampRate)
    stopPnt = centerPnt + round(Int, winLength/2*sampRate)

    meanAmp = mean(wave[startPnt:stopPnt])
    return meanAmp
end


function meanERPAmplitude(wave::Union{AbstractVector{T}, AbstractMatrix{T}}, center::Real, centerType::String, winLength::Real, sampRate::Real, epochStart::Real) where {T<:Real}
    #centerType: point or time

    if ndims(wave) > 1
        if in(1, size(wave)) == false
            error("Only 1-dimensional arrays, or 1xN dimensional arrays are allowed")
        end
        wave = vec(wave)
    end

    if centerType == "point"
        centerPnt = center
    elseif centerType == "time"
        centerPnt = (center-epochStart)*sampRate
    end

    startPnt = centerPnt - round(Int, winLength/2*sampRate)
    stopPnt = centerPnt + round(Int, winLength/2*sampRate)

    meanAmp = mean(wave[startPnt:stopPnt])
    return meanAmp
end

"""
Substitute the event table triggers listed in `trigList`
with newTrig

$(SIGNATURES)

#### Arguments

* `eventTable::Dict{String,Any}`: The event table.
* `trigList::AbstractVector{Integer}`: The list of triggers to substitute.
* `newTrig::Integer`: The new trigger used to substitute the triggers in `trigList`.

#### Examples

    ```julia
    rec, evtTab= simulateRecording(events=[1,2,3,4])
    mergeEventTableCodes!(evtTab, [1,2], 1)
    mergeEventTableCodes!(evtTab, [3,4], 3)
    ```

"""
function mergeEventTableCodes!(eventTable::Dict{String,Any}, trigList::AbstractVector{T}, newTrig::Integer) where {T<:Integer}
    eventTable["code"][findin(eventTable["code"], trigList)] = newTrig
    return
end

"""
Find the exponent of the next power of 2 closest to `x`.

$(SIGNATURES)

#### Arguments

* `x::Real`

#### Examples

    ```julia
    nextPowTwo(6)
    nextPowTwo(511)
    isequal(2^(nextPowTwo(6)), 2^3)
    ```

"""
function nextPowTwo(x::Real)
    out = round(Int, ceil(log2(x)))
    return out
end

"""
Remove epochs from a segmented recording.

$(SIGNATURES)

##### Arguments

* `rec::Dict{String,Array{Real,3}}`: The segmented recording
* `toRemove::Dict{String,Array{P,1}}`: List of epochs to remove for each condition

##### Examples

```julia
    epochDur=0.5; preDur=0.15; events=[1,2]; sampRate=256;
    rec, evtTab = simulateRecording(dur=120, epochDur=epochDur, preDur=preDur, events=events)
    segs, nRaw = segment(rec, evtTab, -preDur, epochDur, sampRate)
    segsToReject = Dict{String,Array{Int,1}}()
    segsToReject["1"] = [3,5]
    segsToReject["2"] = [1,2]
    #toRemove = @compat Dict("1" => [3,5], "2" => [2])
    removeEpochs!(segs, segsToReject)
```

"""
function removeEpochs!(rec::Dict{String,Array{T,3}}, toRemove::Dict{String,Array{P,1}}) where {T<:Real, P<:Integer}
    eventList = collect(keys(rec))
    for i=1:length(eventList)
        code = eventList[i]
        rec[code] = deleteSlice3D(rec[code], toRemove[code], 3)
    end
    #return rec
end

"""
Attempt to remove spurious triggers from an event table.

$(SIGNATURES)

##### Arguments

* `eventTable`: The event table. A dictionary with three fields
    * code: trigger codes
    * idx: trigger indexes
    * dur: trigger durations
* `sentTrigs`: The triggers that were actually sent.
* `minTrigDur`: The minimum duration of legitimate triggers.
* `minTrig`: The minimum trigger code.
* `maxTrig`: The maximum trigger code.

##### Returns

* `res_info`: A dictionary with the following fields:
    * `lenFound`: The number of event table triggers after removal of spurious triggers.
    * `lenSent`: The number of triggers actually sent.
    * `match`: Whether the event table triggers match the sent triggers after removal of spurious triggers.

##### Notes

Spurious triggers may be cause by hardware mulfunction. They may also occurr
if the triggers are not sent digitally (e.g. through the parallel port), but are
sent through an analogous device (e.g. through the soundcard to synchronize
them directly with sound onset). Spurious triggers may have completely different
codes (e.g. 177, when the triggers that were actually sent could only take the
values of 100 and 120), in which case they are easy to find and remove. However,
they may also have the same code as legitimate triggers. In this case they are
more difficult to find. This function can find them if they have a shorter
duration than legitimate triggers.

##### Examples

```julia
    using BrainWave, Compat

    sentTrigs = [1,1,1,2,2,1,2,2,1,1] #triggers that were actually sent
    evtTab = Dict{String,Any}() #fictitious event table
    evtTab["code"] = [1,1,5,1,7,2,5,2,1,2,8,2,1,1,7] #with spurious triggers
    evtTab["idx"] = collect(1:500:500*15)
    evtTab["dur"] = [0.001 for i=1:15]
    res_info = removeSpuriousTriggers!(evtTab, sentTrigs, 0.0004, 192, 254)
    println(res_info)
    assert(res_info["match"] == true)
    #sometimes spurious triggers have the same code as legitimate ones
    #but they can still be found if they differ from legitimate
    #triggers in duration
    sentTrigs = [1,2,3,4,5,6,7,8,9,10] #triggers that were actually sent
    evtTab = Dict{String,Any}() #fictitious event table
    evtTab["code"] = [1,1,1,2,3,4,4,5,6,7,7,7,8,9,10] #with spurious triggers
    evtTab["idx"] = collect(1:500:500*15)
    evtTab["dur"] = [0.001 for i=1:15]
    evtTab["dur"][[2,3,7,11,12]] = 0.0001 #spurious triggers have shorter duration

    res_info = removeSpuriousTriggers!(evtTab, sentTrigs, 0.0004, 192, 254)
    println(res_info)
    assert(res_info["match"] == true)

```

"""
function removeSpuriousTriggers!(eventTable::Dict{String, Any}, sentTrigs::Array{Int}, minTrigDur::Real, minTrig::Real, maxTrig::Real)
    recTrigs = eventTable["code"]
    recTrigsStart = eventTable["idx"]
    recTrigsDur = eventTable["dur"]

    orig_len = length(recTrigs[(recTrigs .<maxTrig) .& (recTrigs .>minTrig)])

    allowedTrigs = round.(Int16, unique(sentTrigs))
    #allowedIdx = findin(recTrigs, allowedTrigs)
    allowedIdx = findall((in)(allowedTrigs), recTrigs)

    recTrigsDur = recTrigsDur[allowedIdx]
    recTrigsStart = recTrigsStart[allowedIdx]
    recTrigs = recTrigs[allowedIdx]

    len_after_notsent_removed = length(recTrigs)

    durCondition = recTrigsDur .>= minTrigDur
    recTrigs = recTrigs[durCondition]
    recTrigsStart = recTrigsStart[durCondition]
    recTrigsDur = recTrigsDur[durCondition]

    foo = zeros(Int16, size(sentTrigs)[1])
    for i=1:size(sentTrigs)[1]
        foo[i] = round(Int16, sentTrigs[i])
    end

    if recTrigs == foo
        match_found = true
    else
        match_found = false
    end


    eventTable["code"] = recTrigs
    eventTable["dur"] = recTrigsDur
    eventTable["idx"] = recTrigsStart

    resInfo = Dict{String,Any}() #(String => Any)[]
    resInfo["match"] = match_found
    resInfo["lenSent"] = length(sentTrigs)
    resInfo["lenFound"] = length(recTrigs)
    #resInfo["origLength"] = orig_len
    #resInfo["len_after_notsent_removed"] = len_after_notsent_removed

    return resInfo
end


"""
$(SIGNATURES)
"""
function RMS(ave::Array{T,2}, winStart::Real, winStop::Real, preDur::Real, sampRate::Real) where {T<:Real}
    nSamp = size(ave)[2]
    nChans = size(ave)[1]

    winStartSample = round(Int, (winStart+preDur)*sampRate)+1
    winStopSample = round(Int, (winStop+preDur)*sampRate)
    if winStartSample < 1
        winStartSample = 1
    end
    if winStopSample > nSamp
        winStopSample = nSamp
    end
    ## println(winStartSample)
    ## println(winStopSample)
    RMSVal = zeros(nChans)
    for n=1:nChans
        RMSVal[n] = sqrt(mean(ave[n, winStartSample:winStopSample].^2))
    end
    return RMSVal
end


"""
Rereference channels in a continuous recording.

$(SIGNATURES)

##### Arguments

* `rec::AbstractMatrix{Real}`: EEG recording.
* `refChan::Integer`: The reference channel number.
* `channels::Union{P, AbstractVector{P}}`: The channel(s) to be rereferenced.

##### Examples

```julia
    rec, evtTab = simulateRecording(nChans=4)
    rerefCnt!(rec, 4, channels=[1, 2, 3])
```

"""
function rerefCnt!(rec::AbstractMatrix{T}, refChan::Integer; channels::Union{P, AbstractVector{P}}=collect(1:size(rec, 1))) where {T<:Real, P<:Integer}

    for i=1:length(channels)
        if channels[i] != refChan
            rec[channels[i],:] = rec[channels[i],:] - rec[refChan,:]
        end
    end

    if in(refChan, channels)
     rec[refChan,:] = rec[refChan,:] - rec[refChan,:]
    end

    return
end


"""
Segment a continuous EEG recording into discrete event-related epochs.

$(SIGNATURES)

##### Arguments

* `rec::AbstractMatrix{T<:Real}`: The nChannelsXnSamples array with the EEG data.
* `eventTable`: dictionary with the following keys
    * `code::AbstractVector{T<:Integer}`: The list of triggers in the EEG recording.
    * `idx::AbstractVector{T<:Integer}`: The indexes of `trigs` in the EEG recording.
* epochStart::FloatingPoint`: The time at which the epoch starts relative to the trigger code, in seconds.
* epochEnd::FloatingPoint`: The time at which the epoch ends relative to the trigger code, in seconds.
* sampRate::Integer`: The sampling rate of the EEG recording.
* eventList::AbstractVector{T<:Integer}`: The list of events for which epochs should be extracted.
        If no list is given epochs will be extracted for all the trigger
        codes present in the event table.

##### Returns

* `segs::Dict{String,Array{T,3}}`: The segmented recording.
        The dictionary has a key for each condition.
        The corresponding key value is a 3D array with dimensions
        nChannels x nSamples x nSegments
* `nSegs::Dict{String,Int64}`: The number of segments for each condition.

##### Examples

```julia
    epochDur=0.5; preDur=0.2; events=[1,2]; sampRate=256;
    rec, evtTab = simulateRecording(dur=120, epochDur=epochDur, preDur=preDur, events=events, sampRate=sampRate)
    segs, nSegs = segment(rec, evtTab, -preDur, epochDur, sampRate, eventsList=[1, 2], eventsLabelsList=["cnd1", "cnd2"])
```

"""
function segment(rec::AbstractMatrix{T},
                 eventTable::Dict{String, Any},
                 epochStart::Real, epochEnd::Real,
                 sampRate::Integer;
                 eventsList::AbstractVector{P}=unique(eventTable["code"]),
                 eventsLabelsList::AbstractVector{S}=String[string(eventsList[i]) for i=1:length(eventsList)]) where {T<:Real, P<:Integer, S<:String}

    trigs = eventTable["code"]
    trigs_pos = eventTable["idx"]

    epochStartSample = round(Int, epochStart*sampRate)
    epochEndSample = round(Int, epochEnd*sampRate) - 1

    nSamples = epochEndSample - epochStartSample + 1
    segs = Dict{String,Array{eltype(rec),3}}() 
    for i=1:length(eventsList)
        idx = trigs_pos[trigs .== eventsList[i]]

        segs[eventsLabelsList[i]] = zeros(eltype(rec), size(rec)[1], nSamples, length(idx))
        for j=1:length(idx)
            thisStartPnt = (idx[j]+epochStartSample)
            #print(thisStartPnt)
            thisStopPnt = (idx[j]+epochEndSample)
            if thisStartPnt < 0 || thisStopPnt > size(rec)[2]
                if thisStartPnt < 0
                    print(idx[j], " Epoch starts before start of recording. \n")
                end
                if thisStopPnt > size(rec)[2]
                    print(idx[j], " Epoch ends after end of recording. \n")
                end

            else
                segs[eventsLabelsList[i]][:,:,j] = rec[:, thisStartPnt:thisStopPnt]

            end
        end
    end
    nSegs = Dict{String,Int}() 
    for i=1:length(eventsList) #count
        nSegs[eventsLabelsList[i]] = size(segs[eventsLabelsList[i]])[3]
    end

    return segs, nSegs

end

"""
Generate a simulated EEG recordings. Not physiologically plausible.
Mainly useful for testing purposes.

$(SIGNATURES)

##### Arguments

* `nChans`: Number of channels.
* `dur`: duration of the recordings, in seconds.
* `sampRate`: Sampling rate.
* `events`: A list of event codes.
* `epochDur`: Duration of an ERP epoch.
* `preDur`: Duration of the pre-stimulus baseline.
* `minVolt`: Minimum possible voltage value.
* `maxVolt`: Maximum possible voltage value.

##### Returns

* `rec`: A nChannelsXnSamples data matrix with values randomly drawn from
     a uniform distribution between minVolt and maxVolts.
* `evtTab`: The event table with event codes and indexes of the samples at which
    events start.

##### Examples

```julia
    rec, evtTab = simulateRecording()
    rec, evtTab = simulateRecording(dur=180, events=[1,2,3])
```

"""
function simulateRecording(;nChans::Int=16, dur::Real=120, sampRate::Real=256,
                           events=[1,2], epochDur::Real=0.5, preDur::Real=0.2,
                           minVolt::Real=-150, maxVolt::Real=150)

    if dur<(preDur+epochDur+1)*2
        error("`dur` is too short")
    end

    rec = rand(nChans, sampRate*dur)*(maxVolt-minVolt) .+ minVolt
    startPoints = collect(ceil(Int, preDur*sampRate):(ceil(Int, preDur*sampRate)+ceil(Int, epochDur*sampRate)):(size(rec)[2]-round(Int, sampRate*1)))
    #nEvents = round(Int, dur/(epochDur+preDur)) - 3
    #nCodes = length(events)
    #nEventsPerCode = floor(Int, nEvents/nCodes)
    #evt = repeat(events, outer=[nEventsPerCode])
    #shuffle!(evt)
    evt = zeros(Int, length(startPoints))
    for i=1:length(evt)
        evt[i] = events[rand(1:length(events))]
    end

    evtTab =  Dict{String,Any}("code" => evt,
                               "idx" => startPoints)
    return rec, evtTab
end

"""
`input(prompt::AbstractString="")`

$(SIGNATURES)

Read a string from STDIN. The trailing newline is stripped.

The prompt string, if given, is printed to standard output without a
trailing newline before reading input.
"""
function input(prompt::AbstractString="")
    print(prompt)
    return chomp(readline())
end

function plotRawEEG(EEG::RawEEG; winDur::Real=8, startTime::Real=0, uVRange::Real=200, startChan::Int=1, nChansToPlot::Int=-1,
                    lineWidth::Real=0.5)
    badSegs = plotRawEEG(EEG.data, EEG.sampRate, winDur=winDur, startTime=startTime, uVRange=uVRange, chanLabels=EEG.chanLabels)
    return badSegs
end

function(dats::AbstractMatrix{T}, sampRate::Real;
                             winDur::Real=8, startTime::Real=0, uVRange::Real=200,
                             startChan::Int=1, nChansToPlot::Int=-1,
                             chanLabels=["" for i=1:size(dats[1])],
                             lineWidth::Real=0.5) where T <: Real
    badSegmentStart = (Real)[]
    badSegmentEnd = (Real)[]
    function doPlot()
        nChansToPlot=min(nChansToPlot, nChans)
        if startChan+nChansToPlot-1 > nChans
            startChan=nChans-(nChansToPlot-1)
        end
        startChan = min(max(1, startChan), nChans)
        chansToPlot = reverse(collect(startChan:startChan+nChansToPlot-1), dims=1)
        #chansToPlot = chansToPlot[chansToPlot .<= nChans]
        nChansToPlot = length(chansToPlot)
        println(chansToPlot)
        plt.cla()
        for i=1:nChansToPlot
            yRaw = dats[chansToPlot[i], idxStart:idxEnd]'
            y = yRaw-mean(yRaw)+uVRange*(i-1)
            plt.plot(tArr[idxStart:idxEnd], y, linewidth=lineWidth)
        end
        plt.xlim([startTime, idxEnd/fs])
        plt.gcf()[:gca]()[:get_xaxis]()[:get_major_formatter]()[:set_useOffset](false)
        plt.xlabel("Time (s)")
        plt.yticks(collect(0:uVRange:nChansToPlot*uVRange), chanLabels[chansToPlot])
        plt.grid()
        plt.draw()
    end
    nChans = size(dats)[1]
    if nChansToPlot == -1
        nChansToPlot = nChans
    end
    nChansToPlot=min(nChansToPlot, nChans)
    nSamp = size(dats)[2]
    tArr = collect(0:nSamp-1)/fs
    nSampToShow = round(Int, winDur*fs)
    idxStart = round(Int, startTime*fs)+1
    idxEnd = idxStart+nSampToShow
 
    doPlot()
    
    theEnd = false
    println("Press:\n 'n' for next window \n 'p' for previous window \n 'x' to exit \n 'j' to jump to time \n 'w' to change window length \n 'r' to change aplitude range \n 'm' to mark bad segments")
    while theEnd == false
        inp = input()
        if inp == "x"
            theEnd = true
        else
            if in(inp, ["n", "\e[C"])
                startTime = startTime+winDur
            elseif inp == "nn"
                startTime = startTime+winDur*2
            elseif inp == "nn"
                startTime = startTime+winDur*3
            elseif in(inp, ["p", "\e[D"])
                startTime = startTime-winDur
            elseif inp == "pp"
                startTime = startTime-winDur*2
            elseif inp == "ppp"
                startTime = startTime-winDur*3
            elseif inp == "s"
                ui = input("Input desired start channel:")
                startChan = parse(Int, ui)
            elseif inp == "c"
                ui = input("Input desired number of channels to plot:")
                nChansToPlot = parse(Int, ui)
            elseif in(inp, ["o", "\e[A"])
                startChan = startChan-1
            elseif in(inp, ["k", "\e[B"])
                startChan = startChan+1
            elseif inp == "j"
                ui = input("Input desired start time (s):")
                startTime = parse(Float32, ui)
            elseif inp == "w"
                ui = input("Input desired window duration (s):")
                winDur = parse(Float32, ui)
                nSampToShow = round(Int, winDur*fs)
            elseif inp == "r"
                ui = input("Input desired range (uV):")
                uVRange = parse(Float32, ui)
            elseif inp == "l"
                ui = input("Input desired line width:")
                lineWidth = parse(Float32, ui)
            elseif inp == "m"
                ui = input("Bad segment start time (s):")
                badStart = parse(Float32, ui)
                ui = input("Bad segment end time (s):")
                badEnd = parse(Float32, ui)
                if badEnd <= badStart
                    printl("Bad segment end time cannot be <= than bad segment start time")
                else
                    push!(badSegmentStart, badStart)
                    push!(badSegmentEnd, badEnd)
                end
           end
       end
       idxStart = max(round(Int, startTime*fs)+1, 1)
       idxEnd = idxStart+nSampToShow
       if idxEnd > nSamp
           idxEnd = nSamp
           idxStart = idxEnd-nSampToShow
           startTime = idxStart/fs
       end
       doPlot()
          
   end
   badSegs = Dict{String, AbstractVector{Real}}()
   badSegs["startTime"] = badSegmentStart
   badSegs["endTime"] = badSegmentEnd
   plt.close()

   return badSegs
end
  

end #Module
