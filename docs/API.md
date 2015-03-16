# ElectroJulia

## Exported
---

#### averageAverages{T<:Real, P<:Integer}(aveList::Array{Dict{String, Array{T<:Real, 2}}, N}, nSegments::Array{Dict{String, P<:Integer}, N})
Perform a weighted average of a list of averages. The weight of
each average in the list is determined by the number of segments
from which it was obtained.
    
#### Arguments

* `aveListArray::Array{Dict{String, Array{Real, 2}}}`: The list of averages for each experimental condition.
* `nSegments::Array{Dict{String, Integer}}`: The number of epochs on which each average is based.

#### Returns

* `weightedAve::Dict{String,Array{Real,2}}`: The weighted average of the averages in the list.

#### Examples

```julia
aveAll, nSegsAll = averageAverages(aveList, nCleanByBlock)
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:38](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### averageEpochs{T<:Real}(rec::Dict{String, Array{T<:Real, 3}})
Average the epochs of a segmented recording.

#### Arguments

* `rec::Dict{String,Array{T,3}}`: Dictionary containing the segmented recordings for each condition.
        The segmented recordings consist of 3-dimensional arrays (n_channels x n_samples x n_epochs).

#### Returns

* `ave::Dict{String,Array{Real,2}}`: The averaged epochs for each condition.
* `n_segs::Dict{String,Integer}`: The number of epochs averaged for each condition.
        
#### Examples

```julia
ave, nSegs = averageEpochs(rec)
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:80](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### baselineCorrect!{T<:Real}(rec::Dict{String, Array{T<:Real, 3}}, baselineStart::Real, preDur::Real, sampRate::Integer)
Perform baseline correction by subtracting the average pre-event
voltage from each channel of a segmented recording.

#### Arguments

* `rec::Dict{String,Array{T,3}}`: The segmented recording.
* `baselineStart::Real`: Start time of the baseline window relative to the event onset, in seconds.
                          The absolute value of `baselineStart` cannot be greater than `preDur`.
                          In practice `baselineStart` allows you to define a baseline window shorter
                          than the time window before the experimental event (`preDur`).
* `preDur::Real`: Duration of recording before the experimental event, in seconds.
* `sampRate::Integer`: The samplig rate of the EEG recording.
    
#### Examples

```julia
#baseline window has the same duration of pre_dur
baselineCorrect(rec, -0.2, 0.2, 512)
#now with a baseline shorter than pre_dur
baselineCorrect(rec, -0.15, 0.2, 512)
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:116](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### chainSegments(rec, nChunks::Integer, sampRate::Integer, startTime::Real, endTime::Real, baselineDur::Real, window)
    ## Take a dictionary containing in each key a list of segments, and chain these segments
    ## into chunks of length nChunks
    ## baselineDur is for determining what is the zero point
    ## startTime and endTime are given with reference to the zero point


**source:**
[ElectroJulia/src/ElectroJulia.jl:180](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### deleteSlice2D{T, P<:Integer}(x::AbstractArray{T, 2}, toRemove::Union(P<:Integer, AbstractArray{P<:Integer, 1}), axis::Integer)
Delete a row or a column from a 2-dimensional array.

#### Args

* `x::AbstractMatrix{T}`: the 2-dimensional array from which rows of columns should be deleted.
* `toRemove::Union(Integer, AbstractVector{Integer})`: an integer or a list of integers indicating the rows to remove.
* `axis::Integer`: an integer indicating whether rows or columns should be removed. 1 corresponds to rows, 2 corresponds to columns.

#### Returns

* `x::AbstractMatrix{T}`: a 2-dimansional array with the selected row or columns removed.

#### Examples

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


**source:**
[ElectroJulia/src/ElectroJulia.jl:289](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### deleteSlice3D(x, toRemove, axis)
Delete a slice from a 3-dimensional array.


**source:**
[ElectroJulia/src/ElectroJulia.jl:311](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### detrendEEG!{T<:Real}(rec::AbstractArray{T<:Real, 2})
Remove the mean value from each channel of an EEG recording.

#### Arguments
* `rec::AbstractMatrix{T}`: The EEG recording.

#### Examples

```julia
x = [1 2 3; 4 5 6]
detrendEEG!(x)
```
    


**source:**
[ElectroJulia/src/ElectroJulia.jl:371](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### fftconvolve(x, y, mode)



**source:**
[ElectroJulia/src/ElectroJulia.jl:464](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### filterContinuous!{T<:Real, P<:Real}(rec::AbstractArray{T<:Real, 2}, sampRate::Integer, filterType::String, nTaps::Integer, cutoffs::Union(P<:Real, AbstractArray{P<:Real, 1}))
Filter a continuous EEG recording.
    
#### Arguments

* rec::AbstractMatrix{Real}`: The nChannelsXnSamples array with the EEG recording.
* `sampRate::Integer`: The EEG recording sampling rate.
* `filterType::String`:  The filter type., one of "lowpass", "highpass", or "bandpass".
* `nTaps::Integer`: The number of filter taps.
* `cutoffs::Union(Real, AbstractVector{Real}`:: The filter cutoffs. If "filterType" is "lowpass" or "highpass"
        the "cutoffs" array should contain a single value. If "filterType"
        is bandpass the "cutoffs" array should contain the lower and
        the upper cutoffs in increasing order.
* `channels::Union(Q, AbstractVector{Q})`: The channel number, or the list of channels numbers that should be filtered.
* `transitionWidth::Real`: The width of the filter transition region, normalized between 0-1.
        For a lower cutoff the nominal transition region will go from
        `(1-transitionWidth)*cutoff` to `cutoff`. For a higher cutoff
        the nominal transition region will go from cutoff to
        `(1+transitionWidth)*cutoff`.
        
#### Examples

```julia
filterContinuous(rec, 2048, "highpass", 512, [30], channels=[0,1,2,3], transitionWidth=0.2)
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:406](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### findArtefactThresh{T<:Real, P<:Real}(rec::Dict{String, Array{T<:Real, 3}}, thresh::Union(P<:Real, AbstractArray{P<:Real, 1}))

Find epochs with voltage values exceeding a given threshold.
    
#### Args

* `rec::Dict{String, Array{Real, 3}}`: The segmented recording.
* `thresh::Union(Real, AbstractVector{Real})`: The threshold value(s).
* `chans::array of ints`: The indexes of the channels on which to find artefacts.
* `chanlabels::array of strings`: The labels of the channels on which to find artefacts.
*  `chanList::array of strings`: The names of all the channels.
    
#### Returns

    
#### Notes

If neither channel indexes (`chans`) nor channel labels (`chanLabels`)
for the channels on which to check artefacts are provided, then artefacts
will be checked for on all channels.
If channel indexes (`chans`) are provided, then channel labels
(`chanLabels`) will be ignored.
If channel labels (`chanLabels`) are given for the channels on which to
check for artefacts, then a list containing the names of all available
channels (`chanList`) must be provided as well.

`thresh` should be a list of threshold values, one for each channel to check.
If `thresh` contains only one value, it is assumed that this is the desired
threshold for all of the channels to check. If `thresh` contains more than one
value, its length must match the number of channels to check.
    
#### Examples




**source:**
[ElectroJulia/src/ElectroJulia.jl:513](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### getACF(sig, sampRate::Real, maxLag::Real)
Compute the autocorrelation function
Arguments:
sig: the signal for which the autocorrelation should be computed
samprate: the sampling rate of the signal
maxLag: the maximum lag (1/f) for which the autocorrelation function should be computed
normalize: whether the autocorrelation should be scaled between [-1, 1]

Returns
acf: the autocorrelation function
lags: the time lags for which the autocorrelation function was computed

n = length(sig)
acfArray = zeros(n*2)
acfArray[1:n] = sig
out = zeros(n)

maxLagPnt = int(round(maxLag*sampRate))
if maxLagPnt > n
maxLagPnt = n
end

for i = 1:maxLagPnt
out[i] = sum(acfArray[1:n] .* acfArray[i:(n+i-1)])
end

lags = [1:maxLagPnt]./sampRate

if normalize == true
out = out ./ maximum(out)
end
return out, lags


**source:**
[ElectroJulia/src/ElectroJulia.jl:597](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### getAutocorrelogram(sig, sampRate::Integer, winLength, overlap, maxLag)
sig: the signal for which the autocorrelogram should be computed
sampRate: the sampling rate of the signal
winLength = the length of the sliding window over which to take the autocorrelations
overlap: overlap between successive windows, in percent
maxLag: the maximum lag for which to compute the autocorrelations
normalize: if `true` divide the output by the maximum ACF value so that ACF values range between 0 and 1
window: the window to be applied to each segment before computing the ACF, defauls to `rect` which does nothing


**source:**
[ElectroJulia/src/ElectroJulia.jl:629](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### getFRatios(ffts, freqs, nSideComp, nExcludedComp, otherExclude)


**source:**
[ElectroJulia/src/ElectroJulia.jl:654](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### getNoiseSidebands(freqs, nCompSide, nExcludedComp, fftDict, otherExclude)
the 2 has the possibility to exclude extra components, useful for distortion products
components: a list containing the indexes of the target components
nCompSide: number of components used for each side band
n_exclude_side: number of components adjacent to to the target components to exclude
fft_array: array containing the fft values


**source:**
[ElectroJulia/src/ElectroJulia.jl:709](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### getPhaseSpectrum(sig, sampRate::Integer)


**source:**
[ElectroJulia/src/ElectroJulia.jl:856](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### getSNR(spec, freqArr, sigFreq, nSideComp, nExclude)


**source:**
[ElectroJulia/src/ElectroJulia.jl:768](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### getSNR2(spec, freqArr, sigFreq, nSideComp, nExclude)
like getSNR, but return signal and noise magnitude separately



**source:**
[ElectroJulia/src/ElectroJulia.jl:783](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### getSpectrogram(sig, sampRate::Integer, winLength::Real, overlap::Real)
winLength in seconds
overlap in percent
if the signal length is not a multiple of the window length it is trucated


**source:**
[ElectroJulia/src/ElectroJulia.jl:799](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### getSpectrum(sig, sampRate::Integer)


**source:**
[ElectroJulia/src/ElectroJulia.jl:820](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### mergeEventTableCodes!{T<:Integer}(eventTable::Dict{String, Any}, trigList::AbstractArray{T<:Integer, 1}, newTrig::Integer)
Substitute the event table triggers listed in `trigList`
with newTrig

#### Parameters

* `eventTable::Dict{String,Any}`: The event table.
* `trigList::AbstractVector{Integer}`: The list of triggers to substitute.
* `newTrig::Integer`: The new trigger used to substitute the triggers in `trigList`.

#### Examples

```julia
mergeEventTableCodes!(evtTab, [200, 220], 999)
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:893](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### nextPowTwo(x::Real)
Find the exponent to which 2 should be raise to find the number corresponding to the next power of 2 closest to `x`.

#### Arguments

* `x::Real`

#### Examples

nextPowTwo(6)
nextPowTwo(511)
isequal(2^(nextPowTwo(6)), 2^3)


**source:**
[ElectroJulia/src/ElectroJulia.jl:911](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### removeEpochs!{T<:Real, P<:Integer}(rec::Dict{String, Array{T<:Real, 3}}, toRemove::Dict{String, Array{P<:Integer, 1}})
Remove epochs from a segmented recording.
    
#### Parameters

* `rec::Dict{String,Array{Real,3}}`: The segmented recording
* `toRemove::Dict{String,Array{P,1}}`: List of epochs to remove for each condition

#### Examples

```julia
removeEpochs!(segs, toRemoveDict)
```



**source:**
[ElectroJulia/src/ElectroJulia.jl:931](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### removeSpuriousTriggers!(eventTable::Dict{String, Any}, sentTrigs::Array{Int64, N}, minTrigDur::Real)

#### Examples

```julia
res_info = removeSpuriousTriggers!(evtTab, behav_trigs, 0.0004)
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:948](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### rerefCnt!{T<:Real}(rec::AbstractArray{T<:Real, 2}, refChan::Integer)
Rereference channels in a continuous recording.

#### Parameters

* `rec::AbstractMatrix{Real}`: EEG recording.
* `refChan::Integer`: The reference channel number.
* `channels::Union(P, AbstractVector{P})`: The channel(s) to be rereferenced.
        
#### Examples

```julia
rerefCnt!(dats, refChan=4, channels=[1, 2, 3])
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:1010](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### segment{T<:Real}(rec::AbstractArray{T<:Real, 2}, eventTable::Dict{String, Any}, epochStart::Real, epochEnd::Real, sampRate::Integer)
Segment a continuous EEG recording into discrete event-related epochs.
    
#### Parameters

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

#### Returns

* `segs::Dict{String,Array{T,3}}`: The segmented recording.
        The dictionary has a key for each condition.
        The corresponding key value is a 3D array with dimensions
        nChannels x nSamples x nSegments
* `nSegs::Dict{String,Int64}`: The number of segments for each condition.
        
#### Examples

```julia
segs, nSegs = segment(dats, evtTab, -0.2, 0.8, 512, eventsList=[200, 201], eventsLabelsList=["cnd1", "cnd2"])
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:1056](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

## Internal
---

#### _centered(arr, newsize)


**source:**
[ElectroJulia/src/ElectroJulia.jl:453](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

---

#### baselineCorrectloop!{T<:Real}(rec::Dict{String, Array{T<:Real, 3}}, baselineStart::Real, preDur::Real, sampRate::Integer)
Perform baseline correction by subtracting the average pre-event
voltage from each channel of a segmented recording.

#### Parameters
* `rec::Dict{String,Array{Real,3}}`: The segmented recording.
* `baselineStart::Real`: Start time of the baseline window relative to the event onset, in seconds.
                          The absolute value of `baselineStart` cannot be greater than `preDur`.
                          In practice `baselineStart` allows you to define a baseline window shorter
                          than the time window before the experimental event (`preDur`).
* `preDur::Real`: Duration of recording before the experimental event, in seconds.
* `sampRate::Integer`: The samplig rate of the EEG recording.
    
#### Examples

```julia
#baseline window has the same duration of pre_dur
baselineCorrect(rec, -0.2, 0.2, 512)
#now with a baseline shorter than pre_dur
baselineCorrect(rec, -0.15, 0.2, 512)
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:155](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)


