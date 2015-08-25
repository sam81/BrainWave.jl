# BrainWave

## Exported

---

<a id="method__averageaverages.1" class="lexicon_definition"></a>
#### averageAverages{T<:Real, P<:Integer}(aveList::Array{Dict{String, Array{T<:Real, 2}}, N}, nSegments::Array{Dict{String, P<:Integer}, N}) [¶](#method__averageaverages.1)
Perform a weighted average of a list of averages. The weight of
each average in the list is determined by the number of segments
from which it was obtained.
    
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



*source:*
[BrainWave/src/BrainWave.jl:48](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__averageepochs.1" class="lexicon_definition"></a>
#### averageEpochs{T<:Real}(rec::Dict{String, Array{T<:Real, 3}}) [¶](#method__averageepochs.1)
Average the epochs of a segmented recording.

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



*source:*
[BrainWave/src/BrainWave.jl:94](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__baselinecorrect.1" class="lexicon_definition"></a>
#### baselineCorrect!{T<:Real}(rec::Dict{String, Array{T<:Real, 3}}, baselineStart::Real, preDur::Real, sampRate::Integer) [¶](#method__baselinecorrect.1)
Perform baseline correction by subtracting the average pre-event
voltage from each channel of a segmented recording.

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



*source:*
[BrainWave/src/BrainWave.jl:135](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__deleteslice2d.1" class="lexicon_definition"></a>
#### deleteSlice2D{T, P<:Integer}(x::AbstractArray{T, 2}, toRemove::Union(AbstractArray{P<:Integer, 1}, P<:Integer), axis::Integer) [¶](#method__deleteslice2d.1)
Delete a row or a column from a 2-dimensional array.

##### Arguments

* `x::AbstractMatrix{T}`: the 2-dimensional array from which rows of columns should be deleted.
* `toRemove::Union(Integer, AbstractVector{Integer})`: an integer or a list of integers indicating the rows/columns to remove.
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



*source:*
[BrainWave/src/BrainWave.jl:184](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__deleteslice3d.1" class="lexicon_definition"></a>
#### deleteSlice3D{T, P<:Integer}(x::Array{T, 3}, toRemove::Union(AbstractArray{P<:Integer, 1}, P<:Integer), axis::Integer) [¶](#method__deleteslice3d.1)
Delete a slice from a 3-dimensional array.

##### Arguments

* `x::Array{T, 3}`: the 3-dimensional array from which rows of columns should be deleted.
* `toRemove::Union(Integer, AbstractVector{Integer})`: an integer or a list of integers indicating the slices to remove.
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



*source:*
[BrainWave/src/BrainWave.jl:228](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__detrendeeg.1" class="lexicon_definition"></a>
#### detrendEEG!{T<:Real}(rec::AbstractArray{T<:Real, 2}) [¶](#method__detrendeeg.1)
Remove the mean value from each channel of an EEG recording.

##### Arguments

* `rec::AbstractMatrix{T}`: The EEG recording.

##### Examples

```julia
    x = [1 2 3; 4 5 6]
    detrendEEG!(x)
```



*source:*
[BrainWave/src/BrainWave.jl:268](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__fftconvolve.1" class="lexicon_definition"></a>
#### fftconvolve{T<:Real, R<:Real}(x::AbstractArray{T<:Real, 1}, y::AbstractArray{R<:Real, 1}, mode::String) [¶](#method__fftconvolve.1)
Convolve two 1-dimensional arrays using the FFT.

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



*source:*
[BrainWave/src/BrainWave.jl:386](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__filtercontinuous.1" class="lexicon_definition"></a>
#### filterContinuous!{T<:Real, P<:Real}(rec::AbstractArray{T<:Real, 2}, sampRate::Integer, filterType::String, nTaps::Integer, cutoffs::Union(AbstractArray{P<:Real, 1}, P<:Real)) [¶](#method__filtercontinuous.1)
Filter a continuous EEG recording.
    
##### Arguments

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
        
##### Examples

```julia
    sampRate = 2048; nTaps=512
    rec, evtTab = simulateRecording(nChans=4, dur=120, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "highpass", nTaps, [30], channels=[1,2,3,4], transitionWidth=0.2)
```



*source:*
[BrainWave/src/BrainWave.jl:306](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__findartefactthresh.1" class="lexicon_definition"></a>
#### findArtefactThresh{T<:Real, P<:Real, Q<:Integer}(rec::Dict{String, Array{T<:Real, 3}}, thresh::Union(AbstractArray{P<:Real, 1}, P<:Real), channels::Union(Q<:Integer, AbstractArray{Q<:Integer, 1})) [¶](#method__findartefactthresh.1)

Find epochs with voltage values exceeding a given threshold.
    
##### Arguments

* `rec::Dict{String, Array{Real, 3}}`: The segmented recording.
* `thresh::Union(Real, AbstractVector{Real})`: The threshold value(s).
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



*source:*
[BrainWave/src/BrainWave.jl:445](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__findartefactthresh.2" class="lexicon_definition"></a>
#### findArtefactThresh{T<:Real, P<:Real}(rec::Dict{String, Array{T<:Real, 3}}, thresh::Union(AbstractArray{P<:Real, 1}, P<:Real)) [¶](#method__findartefactthresh.2)

Find epochs with voltage values exceeding a given threshold.
    
##### Arguments

* `rec::Dict{String, Array{Real, 3}}`: The segmented recording.
* `thresh::Union(Real, AbstractVector{Real})`: The threshold value(s).
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



*source:*
[BrainWave/src/BrainWave.jl:445](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__findextremum.1" class="lexicon_definition"></a>
#### findExtremum{T<:Real}(wave::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), searchStart::Real, searchStop::Real, extremumSign::String, epochStart::Real, sampRate::Real) [¶](#method__findextremum.1)
Find the time point at which a waveform reaches a maximum or a minimum.

##### Arguments:

* `wave::Union(AbstractVector{Real}, AbstractMatrix{Real})`: the waveform for which the extremum should be found.
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
    using Winston
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
    p = plot(tArr, wave1)
    l1 = LineX(timePnt, color="red")
    add(p, l1)
    display(p)
```



*source:*
[BrainWave/src/BrainWave.jl:551](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__getacf.1" class="lexicon_definition"></a>
#### getACF{T<:Real}(sig::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), sampRate::Real, maxLag::Real) [¶](#method__getacf.1)
Compute the autocorrelation function of a 1-dimensional signal.

##### Arguments:

* `sig::Union(AbstractVector{Real}, AbstractMatrix{Real})`: the signal for which the autocorrelation should be computed.
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



*source:*
[BrainWave/src/BrainWave.jl:615](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__getautocorrelogram.1" class="lexicon_definition"></a>
#### getAutocorrelogram{T<:Real}(sig::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), sampRate::Real, winLength::Real, overlap::Real, maxLag::Real) [¶](#method__getautocorrelogram.1)
Compute the autocorrelogram of a 1-dimensional array.
    
##### Arguments

* `sig::Union(AbstractVector{Real}, AbstractMatrix{Real})` The signal of which the autocorrelogram should be computed.
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



*source:*
[BrainWave/src/BrainWave.jl:693](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__getphasespectrum.1" class="lexicon_definition"></a>
#### getPhaseSpectrum{T<:Real}(sig::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), sampRate::Real) [¶](#method__getphasespectrum.1)
Compute the phase spectrum of a 1-dimensional array.
    
##### Arguments

* `sig::Union(AbstractVector{Real}, AbstractMatrix{Real})`: The signal of which the phase spectrum should be computed.
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



*source:*
[BrainWave/src/BrainWave.jl:924](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__getsnr2.1" class="lexicon_definition"></a>
#### getSNR2{T<:Real, R<:Real}(spec::AbstractArray{T<:Real, 1}, freqArr::AbstractArray{R<:Real, 1}, sigFreq::Real, nSideComp::Integer, nExclude::Integer) [¶](#method__getsnr2.1)
Compute the signal-to-noise ratio at a given frequency in the power spectrum of a recording.
This function is the same as `getSNR`, but it additionaly returns the signal and noise magnitudes separately.

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



*source:*
[BrainWave/src/BrainWave.jl:776](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__getsnr.1" class="lexicon_definition"></a>
#### getSNR{T<:Real, R<:Real}(spec::AbstractArray{T<:Real, 1}, freqArr::AbstractArray{R<:Real, 1}, sigFreq::Real, nSideComp::Integer, nExclude::Integer) [¶](#method__getsnr.1)
Compute the signal-to-noise ratio at a given frequency in the power spectrum of a recording.

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



*source:*
[BrainWave/src/BrainWave.jl:738](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__getspectrogram.1" class="lexicon_definition"></a>
#### getSpectrogram{T<:Real}(sig::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), sampRate::Real, winLength::Real, overlap::Real) [¶](#method__getspectrogram.1)
Compute the spectrogram of a 1-dimensional array.
    
##### Arguments
* `sig::Union(AbstractVector{Real}, AbstractMatrix{Real})` The signal of which the spectrum should be computed.
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
    


*source:*
[BrainWave/src/BrainWave.jl:817](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__getspectrum.1" class="lexicon_definition"></a>
#### getSpectrum{T<:Real}(sig::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), sampRate::Integer) [¶](#method__getspectrum.1)
Compute the power spectrum of a 1-dimensional array.
    
##### Arguments

* `sig::Union(AbstractVector{Real}, AbstractMatrix{Real})`: The signal of which the spectrum should be computed.
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



*source:*
[BrainWave/src/BrainWave.jl:860](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__meanerpamplitude.1" class="lexicon_definition"></a>
#### meanERPAmplitude{T<:Real}(wave::Union(AbstractArray{T<:Real, 1}, AbstractArray{T<:Real, 2}), center::Real, centerType::String, winLength::Real, sampRate::Real) [¶](#method__meanerpamplitude.1)
Compute the mean amplitude of an ERP waveform in a time window centered on a given point.

#### Arguments:

* `wave::Union(AbstractVector{Real}, AbstractMatrix{Real})`: the waveform for which the mean amplitude should be computed.
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
    using Winston
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
    p = plot(tArr, wave1)
    l1 = LineX(timePnt, color="red")
    add(p, l1)
    display(p)
    meanAmp =  meanERPAmplitude(wave1, sampPnt, "point", 0.05, sampRate)
    ```



*source:*
[BrainWave/src/BrainWave.jl:1005](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__mergeeventtablecodes.1" class="lexicon_definition"></a>
#### mergeEventTableCodes!{T<:Integer}(eventTable::Dict{String, Any}, trigList::AbstractArray{T<:Integer, 1}, newTrig::Integer) [¶](#method__mergeeventtablecodes.1)
Substitute the event table triggers listed in `trigList`
with newTrig

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



*source:*
[BrainWave/src/BrainWave.jl:1070](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__nextpowtwo.1" class="lexicon_definition"></a>
#### nextPowTwo(x::Real) [¶](#method__nextpowtwo.1)
Find the exponent of the next power of 2 closest to `x`.

#### Arguments

* `x::Real`

#### Examples

    ```julia
    nextPowTwo(6)
    nextPowTwo(511)
    isequal(2^(nextPowTwo(6)), 2^3)
    ```



*source:*
[BrainWave/src/BrainWave.jl:1091](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__removeepochs.1" class="lexicon_definition"></a>
#### removeEpochs!{T<:Real, P<:Integer}(rec::Dict{String, Array{T<:Real, 3}}, toRemove::Dict{String, Array{P<:Integer, 1}}) [¶](#method__removeepochs.1)
Remove epochs from a segmented recording.
    
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



*source:*
[BrainWave/src/BrainWave.jl:1118](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__removespurioustriggers.1" class="lexicon_definition"></a>
#### removeSpuriousTriggers!(eventTable::Dict{String, Any}, sentTrigs::Array{Int64, N}, minTrigDur::Real, minTrig::Real, maxTrig::Real) [¶](#method__removespurioustriggers.1)

Attempt to remove spurious triggers from an event table.

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
    evtTab = Dict{AbstractString,Any}() #fictitious event table
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
    evtTab = Dict{AbstractString,Any}() #fictitious event table
    evtTab["code"] = [1,1,1,2,3,4,4,5,6,7,7,7,8,9,10] #with spurious triggers
    evtTab["idx"] = collect(1:500:500*15)
    evtTab["dur"] = [0.001 for i=1:15]
    evtTab["dur"][[2,3,7,11,12]] = 0.0001 #spurious triggers have shorter duration

    res_info = removeSpuriousTriggers!(evtTab, sentTrigs, 0.0004, 192, 254)
    println(res_info)
    assert(res_info["match"] == true)

```



*source:*
[BrainWave/src/BrainWave.jl:1191](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__rerefcnt.1" class="lexicon_definition"></a>
#### rerefCnt!{T<:Real}(rec::AbstractArray{T<:Real, 2}, refChan::Integer) [¶](#method__rerefcnt.1)
Rereference channels in a continuous recording.

##### Arguments

* `rec::AbstractMatrix{Real}`: EEG recording.
* `refChan::Integer`: The reference channel number.
* `channels::Union(P, AbstractVector{P})`: The channel(s) to be rereferenced.
        
##### Examples

```julia
    rec, evtTab = simulateRecording(nChans=4)
    rerefCnt!(rec, 4, channels=[1, 2, 3])
```



*source:*
[BrainWave/src/BrainWave.jl:1255](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__segment.1" class="lexicon_definition"></a>
#### segment{T<:Real}(rec::AbstractArray{T<:Real, 2}, eventTable::Dict{String, Any}, epochStart::Real, epochEnd::Real, sampRate::Integer) [¶](#method__segment.1)
Segment a continuous EEG recording into discrete event-related epochs.
    
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



*source:*
[BrainWave/src/BrainWave.jl:1304](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

---

<a id="method__simulaterecording.1" class="lexicon_definition"></a>
#### simulateRecording() [¶](#method__simulaterecording.1)
Generate a simulated EEG recordings. Not physiologically plausible.
Mainly useful for testing purposes.

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



*source:*
[BrainWave/src/BrainWave.jl:1378](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

## Internal

---

<a id="method___centered.1" class="lexicon_definition"></a>
#### _centered(arr, newsize) [¶](#method___centered.1)


*source:*
[BrainWave/src/BrainWave.jl:353](file:///home/sam/.julia/v0.3/BrainWave/src/BrainWave.jl)

