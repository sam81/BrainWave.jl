# BrainWave

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

    aveAll, nSegsAll = averageAverages(aveList, nCleanByBlock)



**source:**
[BrainWave/src/BrainWave.jl:37](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

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


    ave, nSegs = averageEpochs(rec)



**source:**
[BrainWave/src/BrainWave.jl:79](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

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


    #baseline window has the same duration of preDur
    baselineCorrect(rec, -0.2, 0.2, 512)
    #now with a baseline shorter than preDur
    baselineCorrect(rec, -0.15, 0.2, 512)



**source:**
[BrainWave/src/BrainWave.jl:115](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### deleteSlice2D{T, P<:Integer}(x::AbstractArray{T, 2}, toRemove::Union(P<:Integer, AbstractArray{P<:Integer, 1}), axis::Integer)
Delete a row or a column from a 2-dimensional array.

#### Args

* `x::AbstractMatrix{T}`: the 2-dimensional array from which rows of columns should be deleted.
* `toRemove::Union(Integer, AbstractVector{Integer})`: an integer or a list of integers indicating the rows/columns to remove.
* `axis::Integer`: an integer indicating whether rows or columns should be removed. 1 corresponds to rows, 2 corresponds to columns.

#### Returns

* `x::AbstractMatrix{T}`: a 2-dimensional array with the selected row or columns removed.

#### Examples

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



**source:**
[BrainWave/src/BrainWave.jl:306](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### deleteSlice3D{T, P<:Integer}(x::Array{T, 3}, toRemove::Union(P<:Integer, AbstractArray{P<:Integer, 1}), axis::Integer)
Delete a slice from a 3-dimensional array.

#### Args

* `x::Array{T, 3}`: the 3-dimensional array from which rows of columns should be deleted.
* `toRemove::Union(Integer, AbstractVector{Integer})`: an integer or a list of integers indicating the slices to remove.
* `axis::Integer`: an integer indicating the axis along which slices should be removed.

#### Returns

* `x::Array{T, 3}`: a 3-dimensional array with the selected row or columns removed.

#### Examples

    x = reshape([1:27], 3,3,3)
    deleteSlice3D(x, 2, 1)
    deleteSlice3D(x, [2,3], 3)
    isequal(deleteSlice3D(x, [2,3], 3), x[:,:, [1]])



**source:**
[BrainWave/src/BrainWave.jl:346](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### detrendEEG!{T<:Real}(rec::AbstractArray{T<:Real, 2})
Remove the mean value from each channel of an EEG recording.

#### Arguments

* `rec::AbstractMatrix{T}`: The EEG recording.

#### Examples

    x = [1 2 3; 4 5 6]
    detrendEEG!(x)

    


**source:**
[BrainWave/src/BrainWave.jl:385](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### fftconvolve{T<:Real, R<:Real}(x::AbstractArray{T<:Real, 1}, y::AbstractArray{R<:Real, 1}, mode::String)
Convolve two 1-dimensional arrays using the FFT.

#### Arguments

* `x::AbstractVector{T}`: First input.
* `y::AbstractVector{T}`: Second input. Should have the same number of dimensions as `x`; if sizes of `x` and `y` are not equal then `x` has to be the larger array.
* `mode::String`: A string indicating the size of the output:
    * "full": The output is the full discrete linear convolution of the inputs. (Default)
    * "valid": The output consists only of those elements that do not rely on the zero-padding.
    * "same": The output is the same size as `x`, centered with respect to the "full" output.

#### Returns

* `out::AbstractVector{T}`: An 1-dimensional array containing a subset of the discrete linear convolution of `x` with `y`.

#### Examples

    x = rand(1:10, 10)
    y = rand(1:10, 10)
    fftconvolve(x, y, "same")



**source:**
[BrainWave/src/BrainWave.jl:498](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### filterContinuous!{T<:Real, P<:Real}(rec::AbstractArray{T<:Real, 2}, sampRate::Integer, filterType::String, nTaps::Integer, cutoffs::Union(P<:Real, AbstractArray{P<:Real, 1}))
Filter a continuous EEG recording.
    
#### Arguments

* `rec::AbstractMatrix{Real}`: The nChannelsXnSamples array with the EEG recording.
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


    filterContinuous!(rec, 2048, "highpass", 512, [30], channels=[0,1,2,3], transitionWidth=0.2)



**source:**
[BrainWave/src/BrainWave.jl:420](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### findArtefactThresh{T<:Real, P<:Real, Q<:Integer}(rec::Dict{String, Array{T<:Real, 3}}, thresh::Union(P<:Real, AbstractArray{P<:Real, 1}), channels::Union(Q<:Integer, AbstractArray{Q<:Integer, 1}))

Find epochs with voltage values exceeding a given threshold.
    
#### Args

* `rec::Dict{String, Array{Real, 3}}`: The segmented recording.
* `thresh::Union(Real, AbstractVector{Real})`: The threshold value(s).
* `chans::AbstractVector{Int}`: The indexes of the channels on which to find artefacts.
* `chanlabels::AbstractVector{String}`: The labels of the channels on which to find artefacts.
* `chanList::AbstractVector{String}`: The names of all the channels.
    
#### Returns

* `segsToReject::Dict{String,Array{Int64,1}}`: dictionary containing the list of segments to reject for each condition.
    
#### Notes

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
    
#### Examples

    # on all channels
    findArtefactThresh(segs, 65)
    # on channels 1 and 2
    findArtefactThresh(segs, 65, [1,2])
    # on channels FP1 and F4
    findArtefactThresh(segs, 20, ["Fp1", "F4"], chanLabels)




**source:**
[BrainWave/src/BrainWave.jl:553](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### findArtefactThresh{T<:Real, P<:Real}(rec::Dict{String, Array{T<:Real, 3}}, thresh::Union(P<:Real, AbstractArray{P<:Real, 1}))

Find epochs with voltage values exceeding a given threshold.
    
#### Args

* `rec::Dict{String, Array{Real, 3}}`: The segmented recording.
* `thresh::Union(Real, AbstractVector{Real})`: The threshold value(s).
* `chans::AbstractVector{Int}`: The indexes of the channels on which to find artefacts.
* `chanlabels::AbstractVector{String}`: The labels of the channels on which to find artefacts.
* `chanList::AbstractVector{String}`: The names of all the channels.
    
#### Returns

* `segsToReject::Dict{String,Array{Int64,1}}`: dictionary containing the list of segments to reject for each condition.
    
#### Notes

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
    
#### Examples

    # on all channels
    findArtefactThresh(segs, 65)
    # on channels 1 and 2
    findArtefactThresh(segs, 65, [1,2])
    # on channels FP1 and F4
    findArtefactThresh(segs, 20, ["Fp1", "F4"], chanLabels)




**source:**
[BrainWave/src/BrainWave.jl:553](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### getACF{T<:Real}(sig::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), sampRate::Real, maxLag::Real)
Compute the autocorrelation function of a 1-dimensional signal.

#### Arguments:

* `sig::Union(AbstractVector{Real}, AbstractMatrix{Real})`: the signal for which the autocorrelation should be computed.
* `samprate::Real`: the sampling rate of the signal.
* `maxLag::Real`: the maximum lag (1/f) for which the autocorrelation function should be computed.
* `normalize::Bool`: whether the autocorrelation should be scaled between [-1, 1].
* `window::Function`: The type of window to apply to the signal before computing its ACF (see DSP.jl).
                      Choose `rect` if you don't want to apply any window.

#### Returns

* `acf::Array{Real,1}`: the autocorrelation function.
* `lags::Array{Real,1}`: the time lags for which the autocorrelation function was computed.

### Examples

    acf, lags = getACF(sig, sampRate, maxLag, normalize=true, window=hamming)



**source:**
[BrainWave/src/BrainWave.jl:634](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### getAutocorrelogram{T<:Real}(sig::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), sampRate::Real, winLength::Real, overlap::Real, maxLag::Real)
Compute the autocorrelogram of a 1-dimensional array.
    
#### Parameters
* `sig::Union(AbstractVector{Real}, AbstractMatrix{Real})` The signal of which the autocorrelogram should be computed.
* `sampRate::Real`: The sampling rate of the signal.
* `winLength::Real`: The length of the window over which to take the ACF, in seconds.
* `overlap::Real`: The percent of overlap between successive windows (useful for smoothing the autocorrelogram).
* `maxLag::Real`: the maximum lag (1/f) for which the autocorrelation function should be computed.
* `normalize::Bool`: whether the autocorrelation should be scaled between [-1, 1].
* `window::Function`: The type of window to apply to the signal before computing its ACF (see DSP.jl).
                      Choose `rect` if you don't want to apply any window.
        
#### Returns

* `acfMatrix::Array{Real,2}`: the autocorrelogram.
* `lags::Array{Real,1}`: the ACF lags.
* `timeArray::Array{Real,1}`: The time axis.

### Examples

    sig = rand(512)
    acg, lags, t = getAutocorrelogram(sig, 256, 0.02, 30, 0.01)



**source:**
[BrainWave/src/BrainWave.jl:708](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### getPhaseSpectrum{T<:Real}(sig::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), sampRate::Real)
Compute the phase spectrum of a 1-dimensional array.
    
#### Arguments

* `sig::Union(AbstractVector{Real}, AbstractMatrix{Real})`: The signal of which the phase spectrum should be computed.
* `sampRate::Real`: The sampling rate of the signal.
* `window::Function`: The type of window to apply to the signal before computing its FFT (see DSP.jl).
                      Choose `rect` if you don't want to apply any window.
* `powerOfTwo::Bool`: If `true` `sig` will be padded with zeros (if necessary) so that its length is a power of two.
        
#### Returns

* `p::Array{Real,1}`: the phase spectrum of the signal.
* `freqArray::Array{Real,1}`: The FFT frequencies.

#### Examples

    sig = rand(512)
    p, f = getPhaseSpectrum(sig, 256)



**source:**
[BrainWave/src/BrainWave.jl:1142](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### getSNR2{T<:Real, R<:Real}(spec::AbstractArray{T<:Real, 1}, freqArr::AbstractArray{R<:Real, 1}, sigFreq::Real, nSideComp::Integer, nExclude::Integer)
Compute the signal-to-noise ratio at a given frequency in the power spectrum of a recording.
This function is the same as `getSNR`, but it additionaly returns the signal and noise magnitudes separately.

#### Arguments

* `spec::AbstractVector{Real}` The power spectrum of the recording.
* `freqArr::AbstractVector{Real}`: the FFT frequency array.
* `sigFreq::Real`: the signal frequency of interest.
* `nSideComp::Integer`: the number of components adjacent to the signal used to estimate the noise.
* `nExclude::Integer`: the number of components closest to the signal to exclude from the noise estimate.

#### Returns

* `snr::Real`: The signal-to-noise ratio at the target frequency.
* `sigMag::Real`: The signal magnitude.
* `noiseMag::Real`: The noise magnitude.

#### Examples

    snr, sigMag, noiseMag = getSNR2(pspec, freq, 140, 10, 1)



**source:**
[BrainWave/src/BrainWave.jl:1000](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### getSNR{T<:Real, R<:Real}(spec::AbstractArray{T<:Real, 1}, freqArr::AbstractArray{R<:Real, 1}, sigFreq::Real, nSideComp::Integer, nExclude::Integer)
Compute the signal-to-noise ratio at a given frequency in the power spectrum of a recording.

#### Arguments

* `spec::AbstractVector{Real}` The power spectrum of the recording.
* `freqArr::AbstractVector{Real}`: the FFT frequency array.
* `sigFreq::Real`: the signal frequency of interest.
* `nSideComp::Integer`: the number of components adjacent to the signal used to estimate the noise.
* `nExclude::Integer`: the number of components closest to the signal to exclude from the noise estimate.

#### Returns

* `snr::Real`: The signal-to-noise ratio at the target frequency.

#### Examples

    getSNR(pspec, freq, 140, 10, 1)



**source:**
[BrainWave/src/BrainWave.jl:966](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### getSpectrogram{T<:Real}(sig::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), sampRate::Real, winLength::Real, overlap::Real)
Compute the spectrogram of a 1-dimensional array.
    
#### Parameters
* `sig::Union(AbstractVector{Real}, AbstractMatrix{Real})` The signal of which the spectrum should be computed.
* `sampRate::Real`: The sampling rate of the signal.
* `winLength::Real`: The length of the window over which to take the FFTs, in seconds.
* `overlap::Real`: The percent of overlap between successive windows (useful for smoothing the spectrogram).
* `window::Function`: The type of window to apply to the signal before computing its FFT (see DSP.jl).
                      Choose `rect` if you don't want to apply any window.
* `powerOfTwo::Bool`: If `true` `sig` will be padded with zeros (if necessary) so that its length is a power of two.
        
#### Returns

* `powerMatrix::Array{Real, 2}`: the power spectrum for each time window.
* `freqArray::Array{Real, 1}`: The frequency axis.
* `timeArray::Array{Real, 1}`: The time axis.

#### Notes

If the signal length is not a multiple of the window length it is trucated.

#### Examples

    sig = rand(512)
    spec, f, t = getSpectrogram(sig, 256, 0.02, 30)
    


**source:**
[BrainWave/src/BrainWave.jl:1039](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### getSpectrum{T<:Real}(sig::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), sampRate::Integer)
Compute the power spectrum of a 1-dimensional array.
    
#### Arguments

* `sig::Union(AbstractVector{Real}, AbstractMatrix{Real})`: The signal of which the spectrum should be computed.
* `sampRate::Real`: The sampling rate of the signal.
* `window::Function`: The type of window to apply to the signal before computing its FFT (see DSP.jl).
                      Choose `rect` if you don't want to apply any window.
* `powerOfTwo::Bool`: If `true` `sig` will be padded with zeros (if necessary) so that its length is a power of two.
        
#### Returns

* `p::Array{Real,1}`: the power spectrum of the signal.
* `freqArray::Array{Real,1}`: The FFT frequencies.

#### Examples

    sig = rand(512)
    p, f = getSpectrum(sig, 256)



**source:**
[BrainWave/src/BrainWave.jl:1080](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

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
[BrainWave/src/BrainWave.jl:1186](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### nextPowTwo(x::Real)
Find the exponent of the next power of 2 closest to `x`.

#### Arguments

* `x::Real`

#### Examples

```julia
nextPowTwo(6)
nextPowTwo(511)
isequal(2^(nextPowTwo(6)), 2^3)
```


**source:**
[BrainWave/src/BrainWave.jl:1206](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

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
[BrainWave/src/BrainWave.jl:1226](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

#### removeSpuriousTriggers!(eventTable::Dict{String, Any}, sentTrigs::Array{Int64, N}, minTrigDur::Real)

#### Examples

```julia
res_info = removeSpuriousTriggers!(evtTab, behav_trigs, 0.0004)
```


**source:**
[BrainWave/src/BrainWave.jl:1243](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

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
[BrainWave/src/BrainWave.jl:1305](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

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
[BrainWave/src/BrainWave.jl:1351](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

## Internal
---

#### _centered(arr, newsize)


**source:**
[BrainWave/src/BrainWave.jl:467](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

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


    #baseline window has the same duration of pre_dur
    baselineCorrect(rec, -0.2, 0.2, 512)
    #now with a baseline shorter than pre_dur
    baselineCorrect(rec, -0.15, 0.2, 512)



**source:**
[BrainWave/src/BrainWave.jl:154](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)


