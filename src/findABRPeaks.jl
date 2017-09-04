
"""
Select the largest peak among a set of peak candidates
(which must be passed as an argument to the function along with their time of occurrence),
that are within a minimum and maximum latency window. If no suitable
peak is found the function returns `NaN`.

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the waveform for which the peak is sought. 
* `peakPnts::AbstractVector{Real}`: vector with candidate peak points. Peak candidates can be found with the `findPeaks` function.
* `peakTimes::AbstractVector{Real}`: vector with candidate peak times. Peak candidates can be found with the `findPeaks` function.
* `minLat::Real`: minimum latency of the peak.
* `maxLat::Real`: maximum latency of the  peak.

##### Returns

* `peakPoint::Real`: index of the point at which the highest peak is found. If no peak is found `NaN` is returned.
* `peakTime`: time at which the highest peak is found, in seconds. If no peak is found `NaN` is returned.

##### Examples

```julia
    sampRate = 256; nTaps=64; minLat=0.5; maxLat=2.5
    rec, evtTab = simulateRecording(nChans=1, dur=4, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "lowpass", nTaps, [2], channels=[1], transitionWidth=0.2)
    pkPnts, pkTimes = findPeaks(rec[1,:], sampRate)
    peakPnt, peakTime = selectLargestPeakInWindow(rec[1,:], pkPnts, pkTimes, minLat, maxLat)
    ## not run
    ## using PlotlyJS
    ## tArr = collect(1:length(rec[1,:]))/sampRate
    ## s1 = scatter(;x=tArr, y=rec[1,:], name="Waveform")
    ## s2 = scatter(;x=tArr[pkPnts], y=rec[1, pkPnts], mode="markers", name="Peak Candidates")
    ## s3 = scatter(;x=[tArr[peakPnt]], y=[rec[1, peakPnt]], mode="markers", marker_size=10, name="Largest Peak")
    ## shapes = rect([minLat], [maxLat], [0], [1], fillcolor="gray", opacity=0.2, line_width=0, xref="x", yref="paper")
    ## plot([s1, s2, s3], Layout(shapes=shapes, xaxis_title="Time (s)", yaxis_title="Amplitude (a.u.)"))


```

"""
function selectLargestPeakInWindow{T<:Real, P<:Real, S<:Real}(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, peakPnts::AbstractVector{P}, peakTimes::AbstractVector{S}, minLat::Real, maxLat::Real)

    peakCandidates = find((peakTimes .<= maxLat ) .& (peakTimes .>= minLat))
    nCandidates = length(peakCandidates)
    if nCandidates == 0
        peakPnt = NaN
        peakTime = NaN
    elseif nCandidates == 1
        peakPnt = peakPnts[peakCandidates[1]]
        peakTime = peakTimes[peakCandidates[1]]
    else
        idx = find(sig[peakPnts[peakCandidates]] .== maximum(sig[peakPnts[peakCandidates]]))[1]
        peakPnt = peakPnts[peakCandidates[idx]]
        peakTime = peakTimes[peakCandidates[idx]]
    end

    return peakPnt, peakTime

end

"""
Select the largest trough among a set of trough candidates
(which must be passed as an argument to the function along with their time of occurrence),
that are within a minimum and maximum latency window. If no suitable
trough is found the function returns `NaN`.

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the waveform for which the trough is sought. 
* `troughPnts::AbstractVector{Real}`: vector with candidate trough points. Trough candidates can be found with the `findTroughs` function.
* `troughTimes::AbstractVector{Real}`: vector with candidate trough times. Trough candidates can be found with the `findTroughs` function.
* `minLat::Real`: minimum latency of the trough.
* `maxLat::Real`: maximum latency of the  trough.

##### Returns

* `troughPoint::Real`: index of the point at which the highest trough is found. If no trough is found `NaN` is returned.
* `troughTime`: time at which the highest trough is found, in seconds. If no trough is found `NaN` is returned.

##### Examples

```julia
    sampRate = 256; nTaps=64; minLat=0.5; maxLat=2.5
    rec, evtTab = simulateRecording(nChans=1, dur=4, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "lowpass", nTaps, [2], channels=[1], transitionWidth=0.2)
    trPnts, trTimes = findTroughs(rec[1,:], sampRate)
    troughPnt, troughTime = selectLargestTroughInWindow(rec[1,:], trPnts, trTimes, minLat, maxLat)
    ## not run
    ## using PlotlyJS
    ## tArr = collect(1:length(rec[1,:]))/sampRate
    ## s1 = scatter(;x=tArr, y=rec[1,:], name="Waveform")
    ## s2 = scatter(;x=tArr[pkPnts], y=rec[1, pkPnts], mode="markers", name="Trough Candidates")
    ## s3 = scatter(;x=[tArr[troughPnt]], y=[rec[1, troughPnt]], mode="markers", marker_size=10, name="Largest Trough")
    ## shapes = rect([minLat], [maxLat], [0], [1], fillcolor="gray", opacity=0.2, line_width=0, xref="x", yref="paper")
    ## plot([s1, s2, s3], Layout(shapes=shapes, xaxis_title="Time (s)", yaxis_title="Amplitude (a.u.)"))


```

"""

function selectLargestTroughInWindow{T<:Real, P<:Real, S<:Real}(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, troughPnts::AbstractVector{P}, troughTimes::AbstractVector{S}, minLat::Real, maxLat::Real, maxAmp::Real=Inf)

    troughCandidates = find((troughTimes .<= maxLat ) .& (troughTimes .>= minLat))
    #make sure trough amplitude is smaller than peak amplitude
    troughCandidates = troughCandidates[find(sig[troughPnts[troughCandidates]] .< maxAmp)]
    nCandidates = length(troughCandidates)
    if nCandidates == 0
        troughPnt = NaN
        troughTime = NaN
    elseif nCandidates == 1
        troughPnt = troughPnts[troughCandidates[1]]
        troughTime = troughTimes[troughCandidates[1]]
    else
        idx = find(sig[troughPnts[troughCandidates]] .== minimum(sig[troughPnts[troughCandidates]]))[1]
        troughPnt = troughPnts[troughCandidates[idx]]
        troughTime = troughTimes[troughCandidates[idx]]
    end

    return troughPnt, troughTime

end

"""
Select the strongest inflection point among a set of inflection point candidates
(which must be passed as an argument to the function along with their time of occurrence),
that are within a minimum and maximum latency window. If no suitable
inflection point is found the function returns `NaN`.

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the waveform for which the inflection point is sought.
* `dy::AbstractVector{Real}`: first derivative of `sig`.
* `inflectionPnts::AbstractVector{Real}`: vector with candidate inflection points. Inflection point candidates can be found with the `findInflections` function.
* `inflectionTimes::AbstractVector{Real}`: vector with candidate inflection times. Inflection point candidates can be found with the `findInflections` function.
* `minLat::Real`: minimum latency of the inflection.
* `maxLat::Real`: maximum latency of the  inflection.

##### Returns

* `inflectionPoint::Real`: index of the point at which the highest inflection point is found. If no inflection point is found `NaN` is returned.
* `inflectionTime`: time at which the highest inflection point is found, in seconds. If no inflection point is found `NaN` is returned.

##### Examples

```julia
    sampRate = 256; nTaps=64; minLat=1; maxLat=2
    rec, evtTab = simulateRecording(nChans=1, dur=4, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "lowpass", nTaps, [2], channels=[1], transitionWidth=0.2)
    inflPnts, inflTimes = findInflections(rec[1,:], sampRate)
    dy = zeros(size(rec[1,:])); dy[2:end] = diff(rec[1,:])
    inflectionPnt, inflectionTime = selectStrongestInflectionInWindow(rec[1,:], dy, inflPnts, inflTimes, minLat, maxLat)
    ## not run
    ## using PlotlyJS
    ## tArr = collect(1:length(rec[1,:]))/sampRate
    ## s1 = scatter(;x=tArr, y=rec[1,:], name="Waveform")
    ## s2 = scatter(;x=tArr[pkPnts], y=rec[1, pkPnts], mode="markers", name="Inflection Candidates")
    ## s3 = scatter(;x=[tArr[inflectionPnt]], y=[rec[1, inflectionPnt]], mode="markers", marker_size=10, name="Strongest Inflection")
    ## shapes = rect([minLat], [maxLat], [0], [1], fillcolor="gray", opacity=0.2, line_width=0, xref="x", yref="paper")
    ## plot([s1, s2, s3], Layout(shapes=shapes, xaxis_title="Time (s)", yaxis_title="Amplitude (a.u.)"))


```

"""

function selectStrongestInflectionInWindow{T<:Real, Q<:Real, P<:Real, S<:Real}(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, dy::Union{AbstractMatrix{Q}, AbstractVector{Q}}, inflPnts::AbstractVector{P}, inflTimes::AbstractVector{S}, minLat::Real, maxLat::Real; maxAmp::Real=Inf)


    inflectionCandidates = find((inflTimes .<= maxLat ) .& (inflTimes .>= minLat))
    #make sure inflection amplitude is smaller than peak amplitude
    inflectionCandidates = inflectionCandidates[find(sig[inflPnts[inflectionCandidates]] .< maxAmp)]
    nCandidates = length(inflectionCandidates)
    if nCandidates == 0
        inflectionPnt = NaN
        inflectionTime = NaN
    elseif nCandidates == 1
        inflectionPnt = inflPnts[inflectionCandidates[1]]
        inflectionTime = inflTimes[inflectionCandidates[1]]
    else
        #strategy 1: find the minimum of signal (useful as surrogate of finding troughs)
        #idx = find(sig[inflPnts[inflectionCandidates]] .== minimum(sig[inflPnts[inflectionCandidates]]))[1]
        # strategy 2: find point at which first derivative is closer to zero (closer to a trough)
        idx = find(abs.(dy[inflPnts[inflectionCandidates]]) .== minimum(abs.(dy[inflPnts[inflectionCandidates]])))[1]
        inflectionPnt = inflPnts[inflectionCandidates[idx]]
        inflectionTime = inflTimes[inflectionCandidates[idx]]
    end

    return inflectionPnt, inflectionTime

end

"""
Attempt to find peaks and troughs of the ABR response for a click of a given level. The algorithm is
largely based on Bradley and Wilson (2004).

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the ABR waveform for which the peaks and troughs are sought. 
* `stimLevel::Real`: the level of the click used to evoke the ABR response.
* `sampRate::Real`: the sampling rate of the ABR recording.
* `epochStart::Real`: the time, in seconds, at which the epoch starts.
* `waveVLatencyOffset::Real`: additional wave V latency delay (e.g. if stimulus is different than a standard click you may want to specify an additional delay on top of the delay computed by the formula)
* `dBRef::String`: whether the stimulus level is specified in dB SPL `SPL`, or in dB normal hearing level `nHL`.

##### Returns

A dataframe with the following columns:

* `wave::String`: wave label.
* `peakPoint::Real`: index of the point at which the peak is detected in the waveform.
* `troughPoint::Real`: index of the point at which the peak is detected in the waveform.
* `peakLatency::Real`: latency of the peak, in seconds.
* `troughLatency::Real`: latency of the trough, in seconds.
* `peakAmp::Real`: amplitude of the peak, in microvolts.
* `troughAmp::Real`: amplitude of the trough, in microvolts.
* `peakTroughAmp::Real`: peak-to-trough amplitude, in microvolts.
* `minPeakLat::Real`: minimum peak latency, in seconds, used by the algorithm to find the peak.
* `maxPeakLat::Real`: maximum peak latency, in seconds, used by the algorithm to find the peak.
* `minTroughLat::Real`: minimum trough latency, in seconds, used by the algorithm to find the trough.
* `maxTroughLat::Real`: maximum trough latency, in seconds, used by the algorithm to find the trough.



##### References

* Bradley, A. P., & Wilson, W. J. (2004). Automated Analysis of the Auditory Brainstem Response. Proceedings of the 2004 Intelligent Sensors, Sensor Networks and Information Processing Conference, 2004., 541–546. http://doi.org/10.1109/ISSNIP.2004.1417519

* Prosser, S., & Arslan, E. (1987). Prediction of auditory brainstem wave V latency as a diagnostic tool of sensorineural hearing loss. Audiology : Official Organ of the International Society of Audiology, 26(3), 179–87.

##### Examples

```julia

```

"""
function findABRWaves{T<:Real}(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, stimLevel::Real, sampRate::Real; epochStart::Real=0, waveVLatencyOffset::Real=0, dBRef::String="SPL")
    sig = vec(sig)
    if dBRef == "SPL"
        stimLevel = stimLevel - 31
    elseif dBRef == "nHL"
        stimLevel = stimLevel
    else
        error("`stimLevel` needs to be specified in `SPL` or `nHL`")
    end

    #equation derived from Prosser, S., & Arslan, E. (1987). Prediction of auditory brainstem wave V latency as a diagnostic tool of sensorineural hearing loss. Audiology : Official Organ of the International Society of Audiology, 26(3), 179–87. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/3662941
    avPeakVLat = (-5.859e-08*stimLevel^4 + 1.274e-05*stimLevel^3 - 0.0005424*stimLevel^2 - 0.05297*stimLevel + 9.3)/1000
    avPeakVLat = avPeakVLat + waveVLatencyOffset
    
    minPeakVLat = avPeakVLat-(0.00024*3)
    maxPeakVLat = avPeakVLat+(0.00024*3)
    #todo 
    #shorten peak II trough II latency

    ##data from Picton book
    avPeakVpeakI_IPL = 0.00395
    avPeakVpeakIII_IPL = 0.00207

    minPeakVpeakI_IPL = avPeakVpeakI_IPL - (0.00023*2)
    maxPeakVpeakI_IPL = avPeakVpeakI_IPL + (0.00023*2)
    minPeakVpeakIII_IPL = avPeakVpeakIII_IPL - (0.00019*2)
    maxPeakVpeakIII_IPL = avPeakVpeakIII_IPL + (0.00019*2)

    ## guesstimates
    minPeakITroughILat = 0.25/1000
    maxPeakITroughILat = minPeakITroughILat+0.75/1000

    minPeakIIITroughIIILat = 0.25/1000
    maxPeakIIITroughIIILat = minPeakIIITroughIIILat+0.75/1000

    minPeakVTroughVLat = 0.25/1000
    maxPeakVTroughVLat = 1.75/1000

    minPeakIITroughIILat = 0.125/1000
    minPeakIVTroughIVLat = 0.125/1000

    # find all peaks, troughs and inflection points in the waveform
    peakPnts, peakTimes = findPeaks(sig, sampRate, epochStart=epochStart)
    troughPnts, troughTimes = findTroughs(sig, sampRate, epochStart=epochStart)
    inflPnts, inflTimes = findInflections(sig, sampRate, epochStart=epochStart)
    #dy = zeros(size(sig)); dy[:,2:size(sig)[2]] = diff(sig, 2)
    dy = zeros(size(sig)); dy[2:end] = diff(sig)


    peakVPnt, peakVTime = selectLargestPeakInWindow(sig, peakPnts, peakTimes, minPeakVLat, maxPeakVLat)
    if isnan(peakVPnt) == true
        peakIPnt = NaN
        peakITime = NaN
        peakIIPnt = NaN
        peakIITime = NaN
        peakIIIPnt = NaN
        peakIIITime = NaN
        peakIVPnt = NaN
        peakIVTime = NaN

        troughIPnt = NaN
        troughITime = NaN
        troughIIPnt = NaN
        troughIITime = NaN
        troughIIIPnt = NaN
        troughIIITime = NaN
        troughIVPnt = NaN
        troughIVTime = NaN
        troughVPnt = NaN
        troughVTime = NaN
        minPeakILat = NaN
        maxPeakILat = NaN
        minPeakIILat = NaN
        maxPeakIILat = NaN
        minPeakIIILat = NaN
        maxPeakIIILat = NaN
        minPeakIVLat = NaN
        maxPeakIVLat = NaN

        minTroughILat = NaN
        maxTroughILat = NaN
        minTroughIILat = NaN
        maxTroughIILat = NaN
        minTroughIIILat = NaN
        maxTroughIIILat = NaN
        minTroughIVLat = NaN
        maxTroughIVLat = NaN
        minTroughVLat = NaN
        maxTroughVLat = NaN
    else
        minPeakILat = peakVTime - maxPeakVpeakI_IPL
        maxPeakILat = peakVTime - minPeakVpeakI_IPL
        minPeakIIILat = peakVTime - maxPeakVpeakIII_IPL
        maxPeakIIILat = peakVTime - minPeakVpeakIII_IPL

        minTroughVLat = peakVTime + minPeakVTroughVLat
        maxTroughVLat = peakVTime + maxPeakVTroughVLat

        #Peak I
        peakIPnt, peakITime = selectLargestPeakInWindow(sig, peakPnts, peakTimes, minPeakILat, maxPeakILat)
        #Peak III
        peakIIIPnt, peakIIITime = selectLargestPeakInWindow(sig, peakPnts, peakTimes, minPeakIIILat, maxPeakIIILat)
        
        #Trough V
        troughVPnt, troughVTime = selectLargestTroughInWindow(sig, troughPnts, troughTimes, minTroughVLat, maxTroughVLat, sig[peakVPnt])
        if isnan(troughVPnt)
            troughVPnt, troughVTime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minTroughVLat, maxTroughVLat, maxAmp=sig[peakVPnt])
        end

        #Trough I
        if isnan(peakIPnt) == true
            troughIPnt = NaN
            troughITime = NaN
            minTroughILat = NaN
            maxTroughILat = NaN
        else
            minTroughILat = peakITime + minPeakITroughILat
            maxTroughILat = peakITime + maxPeakITroughILat
            troughIPnt, troughITime = selectLargestTroughInWindow(sig, troughPnts, troughTimes, minTroughILat, maxTroughILat, sig[peakIPnt])
            if isnan(troughIPnt)
                troughIPnt, troughITime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minTroughILat, maxTroughILat, maxAmp=sig[peakIPnt])
            end
        end

        #Trough III
        if isnan(peakIIIPnt) == true
            troughIIIPnt = NaN
            troughIIITime = NaN
            minTroughIIILat = NaN
            maxTroughIIILat = NaN
        else
            minTroughIIILat = peakIIITime + minPeakIIITroughIIILat
            maxTroughIIILat = peakIIITime + maxPeakIIITroughIIILat
            troughIIIPnt, troughIIITime = selectLargestTroughInWindow(sig, troughPnts, troughTimes, minTroughIIILat, maxTroughIIILat, sig[peakIIIPnt])
            if isnan(troughIIIPnt)
                troughIIIPnt, troughIIITime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minTroughIIILat, maxTroughIIILat, maxAmp=sig[peakIIIPnt])
            end
        end



        #Peak II
        minPeakIILat = NaN
        maxPeakIILat = NaN
        if isnan(peakIIITime) == false
            maxPeakIILat = peakIIITime - 0.25/1000
        end
        if isnan(troughITime) == false
            minPeakIILat = troughITime + 0.25/1000
        end
        if isnan(maxPeakIILat) == true
            if isnan(peakVTime) == false
                if isnan(minPeakIILat) == false
                    maxPeakIILat = (peakVTime + minPeakIILat)/2.5
                end
            end
        end

        peakIIPnt=NaN; peakIITime=NaN
        if ((isnan(minPeakIILat) == false) .& (isnan(maxPeakIILat) == false))
            peakIIPnt, peakIITime = selectLargestPeakInWindow(sig, peakPnts, peakTimes, minPeakIILat, maxPeakIILat)
            if isnan(peakIIPnt)
                peakIIPnt, peakIITime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minPeakIILat, maxPeakIILat)
            end
        end
            
         
        #Trough II
        minTroughIILat = NaN
        maxTroughIILat = NaN
        if isnan(peakIIITime) == false
            maxTroughIILat = peakIIITime - 0.25/1000
        end
        if isnan(peakIITime) == false
            minTroughIILat = peakIITime + minPeakIITroughIILat
        end

        troughIIPnt=NaN; troughIITime=NaN
        if ((isnan(minTroughIILat) == false) .& (isnan(maxTroughIILat) == false))
            troughIIPnt, troughIITime = selectLargestTroughInWindow(sig, troughPnts, troughTimes, minTroughIILat, maxTroughIILat, sig[peakIIPnt])
            if isnan(troughIIPnt)
                troughIIPnt, troughIITime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minTroughIILat, maxTroughIILat, maxAmp=sig[peakIIPnt])
            end
        end
            
        #peak IV
        minPeakIVLat = NaN
        maxPeakIVLat = NaN
        if isnan(peakVTime) == false
            maxPeakIVLat = peakVTime - 0.25/1000
        end
        if isnan(troughIIITime) == false
            minPeakIVLat = troughIIITime + 0.25/1000
        end

        peakIVPnt=NaN; peakIVTime=NaN
        if ((isnan(minPeakIVLat) == false) .& (isnan(maxPeakIVLat) == false))
            peakIVPnt, peakIVTime = selectLargestPeakInWindow(sig, peakPnts, peakTimes, minPeakIVLat, maxPeakIVLat)
            if isnan(peakIVPnt)
                peakIVPnt, peakIVTime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minPeakIVLat, maxPeakIVLat)
            end
        end
            
         

        #trough IV
        minTroughIVLat = NaN
        maxTroughIVLat = NaN
        if isnan(peakVTime) == false
            maxTroughIVLat = peakVTime - 0.25/1000
        end
        if isnan(peakIVTime) == false
            minTroughIVLat = peakIVTime + minPeakIVTroughIVLat
        end

        troughIVPnt=NaN; troughIVTime=NaN
        if ((isnan(minTroughIVLat) == false) .& (isnan(maxTroughIVLat) == false))
            troughIVPnt, troughIVTime = selectLargestTroughInWindow(sig, troughPnts, troughTimes, minTroughIVLat, maxTroughIVLat, sig[peakIVPnt])
            if isnan(troughIVPnt)
                troughIVPnt, troughIVTime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minTroughIVLat, maxTroughIVLat, maxAmp=sig[peakIVPnt])
            end
        end
    end

    peakPoints = [peakIPnt; peakIIPnt; peakIIIPnt; peakIVPnt; peakVPnt]
    peakLatencies = [peakITime; peakIITime; peakIIITime; peakIVTime; peakVTime]
    troughPoints = [troughIPnt; troughIIPnt; troughIIIPnt; troughIVPnt; troughVPnt]
    troughLatencies = [troughITime; troughIITime; troughIIITime; troughIVTime; troughVTime]

    peakAmps = zeros(length(peakPoints))
    troughAmps = zeros(length(peakPoints))
    #peakAmpsEstMax = zeros(length(peakPoints))
    #troughAmpsEstMax = zeros(length(peakPoints))
    for i=1:length(peakAmps)
        if isnan(peakPoints[i]) == false
            peakAmps[i] = sig[round(Int, peakPoints[i])]
            #peakAmpsEstMax[i] = sig[round(Int, peakPoints[i])]
        else
            peakAmps[i] = NaN
            #peakAmpsEstMax[i] = maximum(sig[,)
        end

        if isnan(troughPoints[i]) == false
            troughAmps[i] = sig[round(Int, troughPoints[i])]
        else
            troughAmps[i] = NaN
        end
    end


    peakTroughAmps = peakAmps-troughAmps

    waveLabels = ["I", "II", "III", "IV", "V"]

    df = DataFrame(wave=waveLabels, peakPoint=peakPoints, troughPoint=troughPoints, peakLatency=peakLatencies, troughLatency=troughLatencies, peakAmp=peakAmps, troughAmp=troughAmps, peakTroughAmp=peakTroughAmps,
                   minPeakLat=[minPeakILat, minPeakIILat, minPeakIIILat, minPeakIVLat, minPeakVLat],
                   maxPeakLat=[maxPeakILat, maxPeakIILat, maxPeakIIILat, maxPeakIVLat, maxPeakVLat],
                   minTroughLat=[minTroughILat, minTroughIILat, minTroughIIILat, minTroughIVLat, minTroughVLat],
                   maxTroughLat=[maxTroughILat, maxTroughIILat, maxTroughIIILat, maxTroughIVLat, maxTroughVLat])

    return df
end




## """
## Attempt to find peaks and troughs of the ABR response for a click of a given level. The algorithm is
## largely based on Bradley and Wilson (2004).

## $(SIGNATURES)

## ##### Arguments

## * `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the ABR waveform for which the peaks and troughs are sought. 
## * `stimLevel::Real`: the level of the click used to evoke the ABR response.
## * `sampRate::Real`: the sampling rate of the ABR recording.
## * `epochStart::Real`: the time, in seconds, at which the epoch starts.
## * `minInterPeakLatency::Real``: the minimum allowed latency between peaks.
## * `dBRef::String`: whether the stimulus level is specified in dB SPL `SPL`, or in dB normal hearing level `nHL`.

## ##### Returns

## * `peakPoints::`: array containing the points at which peaks I-V are detected
## * `troughPoints`
## * `peakLatencies`
## * `troughLatencies`
## * `peakAmps`
## * `troughAmps`
## * `peakTroughAmps`
## * `prms`


## ##### References

## * Bradley, A. P., & Wilson, W. J. (2004). Automated Analysis of the Auditory Brainstem Response. Proceedings of the 2004 Intelligent Sensors, Sensor Networks and Information Processing Conference, 2004., 541–546. http://doi.org/10.1109/ISSNIP.2004.1417519

## ##### Examples

## ```julia

## ```

## """
## #min interpeak latency = 0.5 ms
## #0 dB nHL = 31 dB ppe SPL
## function findABRPeaks{T<:Real}(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, stimLevel::Real, sampRate::Real; epochStart::Real=0, minInterPeakLatency::Real=0.0005, dBRef::String="SPL")

##     if dBRef == "SPL"
##         stimLevel = stimLevel - 31
##     elseif dBRef == "nHL"
##         stimLevel = stimLevel
##     else
##         error("`stimLevel` needs to be specified in `SPL` or `nHL`")
##     end

##     #equation derived from Prosser, S., & Arslan, E. (1987). Prediction of auditory brainstem wave V latency as a diagnostic tool of sensorineural hearing loss. Audiology : Official Organ of the International Society of Audiology, 26(3), 179–87. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/3662941
##     avPeakVLat = (-5.859e-08*stimLevel^4 + 1.274e-05*stimLevel^3 - 0.0005424*stimLevel^2 - 0.05297*stimLevel + 9.3)/1000
##     minPeakVLat = avPeakVLat-0.001
##     maxPeakVLat = avPeakVLat+0.001

##     #todo enforce troughs to be lower in amplitude than peaks
##     #shorten peak II trough II latency

##     avPeakVpeakI_IPL = 0.0038
##     avPeakVpeakIII_IPL = 0.0017

##     minPeakVpeakI_IPL = avPeakVpeakI_IPL - 0.0002*2
##     maxPeakVpeakI_IPL = avPeakVpeakI_IPL + 0.0002*2
##     minPeakVpeakIII_IPL = avPeakVpeakIII_IPL - 0.00015*2
##     maxPeakVpeakIII_IPL = avPeakVpeakIII_IPL + 0.00015*2

##     minPeakITroughILat = 0.25/1000
##     maxPeakITroughILat = minPeakITroughILat+0.75/1000

##     minPeakIIITroughIIILat = 0.25/1000
##     maxPeakIIITroughIIILat = minPeakIIITroughIIILat+0.75/1000

##     minPeakVTroughVLat = 0.25/1000
##     maxPeakVTroughVLat = minPeakVTroughVLat+1.2/1000

##     minPeakIITroughIILat = 0.12/1000
##     minPeakIVTroughIVLat = 0.12/1000


##     peakPnts, peakTimes = findPeaks(sig, sampRate, epochStart=epochStart)
##     troughPnts, troughTimes = findTroughs(sig, sampRate, epochStart=epochStart)
##     inflPnts, inflTimes = findInflections(sig, sampRate, epochStart=epochStart)
##     dy = zeros(size(sig)); dy[:,2:size(sig)[2]] = diff(sig, 2)
##     #ddy = zeros(size(sig)); ddy[:,2:size(sig)[2]] = diff(dy, 2)


##     peakVCandidates = find((peakTimes .< maxPeakVLat ) .& (peakTimes .> minPeakVLat))
##     nPeakVCandidates = length(peakVCandidates)
##     if nPeakVCandidates == 0
##         peakVPnt = NaN
##         peakVTime = NaN
##     elseif nPeakVCandidates == 1
##         peakVPnt = peakPnts[peakVCandidates[1]]
##         peakVTime = peakTimes[peakVCandidates[1]]
##     else
##         idx = find(sig[peakPnts[peakVCandidates]] .== maximum(sig[peakPnts[peakVCandidates]]))[1]
##         peakVPnt = peakPnts[peakVCandidates[idx]]
##         peakVTime = peakTimes[peakVCandidates[idx]]
##     end

##     if isnan(peakVPnt) == true
##         peakIPnt = NaN
##         peakITime = NaN
##         peakIIIPnt = NaN
##         peakIIITime = NaN
##         troughVPnt = NaN
##         troughVTime = NaN
##         minPeakILat = NaN
##         maxPeakILat = NaN
##         minPeakIIILat = NaN
##         maxPeakIIILat = NaN
##         minTroughVLat = NaN
##         maxTroughVLat = NaN
##     else
##         minPeakILat = peakVTime - maxPeakVpeakI_IPL
##         maxPeakILat = peakVTime - minPeakVpeakI_IPL
##         minPeakIIILat = peakVTime - maxPeakVpeakIII_IPL
##         maxPeakIIILat = peakVTime - minPeakVpeakIII_IPL

##         #Peak I
##         peakICandidates = find((peakTimes .< maxPeakILat ) .& (peakTimes .> minPeakILat))
##         nPeakICandidates = length(peakICandidates)
##         if nPeakICandidates == 0
##             peakIPnt = NaN
##             peakITime = NaN
##         elseif nPeakICandidates == 1
##             peakIPnt = peakPnts[peakICandidates[1]]
##             peakITime = peakTimes[peakICandidates[1]]
##         else
##             idx = find(sig[peakPnts[peakICandidates]] .== maximum(sig[peakPnts[peakICandidates]]))[1]
##             peakIPnt = peakPnts[peakICandidates[idx]]
##             peakITime = peakTimes[peakICandidates[idx]]
##         end

##         #Peak III
##         peakIIICandidates = find((peakTimes .< maxPeakIIILat ) .& (peakTimes .> minPeakIIILat))
##         nPeakIIICandidates = length(peakIIICandidates)
##         if nPeakIIICandidates == 0
##             peakIIIPnt = NaN
##             peakIIITime = NaN
##         elseif nPeakIIICandidates == 1
##             peakIIIPnt = peakPnts[peakIIICandidates[1]]
##             peakIIITime = peakTimes[peakIIICandidates[1]]
##         else
##             idx = find(sig[peakPnts[peakIIICandidates]] .== maximum(sig[peakPnts[peakIIICandidates]]))[1]
##             peakIIIPnt = peakPnts[peakIIICandidates[idx]]
##             peakIIITime = peakTimes[peakIIICandidates[idx]]
##         end

##         #Trough V
##         minTroughVLat = peakVTime + minPeakVTroughVLat
##         maxTroughVLat = peakVTime + maxPeakVTroughVLat

##         troughVCandidates = find((troughTimes .< maxTroughVLat ) .& (troughTimes .> minTroughVLat))
##         troughVCandidates = troughVCandidates[find(sig[troughPnts[troughVCandidates]] .< sig[peakVPnt])] #ensure trough amp is less than peak amp
##         nTroughVCandidates = length(troughVCandidates)
##         if nTroughVCandidates == 0
##             troughVCandidates2 = find((inflTimes .< maxTroughVLat ) .& (inflTimes .> minTroughVLat))
##             troughVCandidates2 = troughVCandidates2[find(sig[inflPnts[troughVCandidates2]] .< sig[peakVPnt])] #ensure trough amp is less than peak amp
##             nTroughVCandidates2 = length(troughVCandidates2)
##             if nTroughVCandidates2 == 0
##                 troughVPnt = NaN
##                 troughVTime = NaN
##             elseif nTroughVCandidates2 == 1
##                 troughVPnt = inflPnts[troughVCandidates2[1]]
##                 troughVTime = inflTimes[troughVCandidates2[1]]
##             else
##                 ##idx = find(sig[inflPnts[troughVCandidates2]] .== minimum(sig[inflPnts[troughVCandidates2]])) #need to check here which point to take
##                 ##idx = find(abs(diff(sig[inflPnts[troughVCandidates2]-1])) .== minimum(abs(diff(sig[inflPnts[troughVCandidates2]-1])))) #need to check here which point to take
##                 idx = find(abs.(dy[inflPnts[troughVCandidates2]]) .== minimum(abs.(dy[inflPnts[troughVCandidates2]])))[1]
##                 troughVPnt = inflPnts[troughVCandidates2[idx]]
##                 troughVTime = inflTimes[troughVCandidates2[idx]]
##             end
##         elseif nTroughVCandidates == 1
##             troughVPnt = troughPnts[troughVCandidates[1]]
##             troughVTime = troughTimes[troughVCandidates[1]]
##         else
##             idx = find(sig[troughPnts[troughVCandidates]] .== minimum(sig[troughPnts[troughVCandidates]]))[1]
##             troughVPnt = troughPnts[troughVCandidates[idx]]
##             troughVTime = troughTimes[troughVCandidates[idx]]
##         end

##     end

##     #Trough III
##     if isnan(peakIIIPnt) == true
##         troughIIIPnt = NaN
##         troughIIITime = NaN
##         minTroughIIILat = NaN
##         maxTroughIIILat = NaN
##     else
##         minTroughIIILat = peakIIITime + minPeakIIITroughIIILat
##         maxTroughIIILat = peakIIITime + maxPeakIIITroughIIILat

##         troughIIICandidates = find((troughTimes .< maxTroughIIILat ) .& (troughTimes .> minTroughIIILat))
##         troughIIICandidates = troughIIICandidates[find(sig[troughPnts[troughIIICandidates]] .< sig[peakIIIPnt])] #ensure trough amp is less than peak amp
##         nTroughIIICandidates = length(troughIIICandidates)
##         if nTroughIIICandidates == 0
##             troughIIICandidates2 = find((inflTimes .< maxTroughIIILat ) .& (inflTimes .> minTroughIIILat))
##             troughIIICandidates2 = troughIIICandidates2[find(sig[inflPnts[troughIIICandidates2]] .< sig[peakIIIPnt])] #ensure trough amp is less than peak amp
##             nTroughIIICandidates2 = length(troughIIICandidates2)
##             if nTroughIIICandidates2 == 0
##                 troughIIIPnt = NaN
##                 troughIIITime = NaN
##             elseif nTroughIIICandidates2 == 1
##                 troughIIIPnt = inflPnts[troughIIICandidates2[1]]
##                 troughIIITime = inflTimes[troughIIICandidates2[1]]
##             else
##                 ##idx = find(sig[inflPnts[troughIIICandidates2]] .== minimum(sig[inflPnts[troughIIICandidates2]])) #need to check here which point to take
##                 ##idx = find(sig[inflPnts[troughIIICandidates2]] .== minimum(diff(sig[inflPnts[troughIIICandidates2]-1]))) #need to check here which point to take
##                 idx = find(abs.(dy[inflPnts[troughIIICandidates2]]) .== minimum(abs.(dy[inflPnts[troughIIICandidates2]])))[1]
##                 troughIIIPnt = inflPnts[troughIIICandidates2[idx]]
##                 troughIIITime = inflTimes[troughIIICandidates2[idx]]
##             end
##         elseif nTroughIIICandidates == 1
##             troughIIIPnt = troughPnts[troughIIICandidates[1]]
##             troughIIITime = troughTimes[troughIIICandidates[1]]
##         else
##             idx = find(sig[troughPnts[troughIIICandidates]] .== minimum(sig[troughPnts[troughIIICandidates]]))[1]
##             troughIIIPnt = troughPnts[troughIIICandidates[idx]]
##             troughIIITime = troughTimes[troughIIICandidates[idx]]
##         end
##     end

##     #Trough I
##     if isnan(peakIPnt) == true
##         troughIPnt = NaN
##         troughITime = NaN
##         minTroughILat = NaN
##         maxTroughILat = NaN
##     else
##         minTroughILat = peakITime + minPeakITroughILat
##         maxTroughILat = peakITime + maxPeakITroughILat

##         troughICandidates = find((troughTimes .< maxTroughILat ) .& (troughTimes .> minTroughILat))
##         troughICandidates = troughICandidates[find(sig[troughPnts[troughICandidates]] .< sig[peakIPnt])] #ensure trough amp is less than peak amp
##         nTroughICandidates = length(troughICandidates)
##         if nTroughICandidates == 0
##             troughICandidates2 = find((inflTimes .< maxTroughILat ) .& (inflTimes .> minTroughILat))
##             troughICandidates2 = troughICandidates2[find(sig[inflPnts[troughICandidates2]] .< sig[peakIPnt])] #ensure trough amp is less than peak amp
##             nTroughICandidates2 = length(troughICandidates2)
##             if nTroughICandidates2 == 0
##                 troughIPnt = NaN
##                 troughITime = NaN
##             elseif nTroughICandidates2 == 1
##                 troughIPnt = inflPnts[troughICandidates2[1]]
##                 troughITime = inflTimes[troughICandidates2[1]]
##             else
##                 ##idx = find(sig[inflPnts[troughICandidates2]] .== minimum(sig[inflPnts[troughICandidates2]])) #need to check here which point to take
##                 ##idx = find(abs(diff(sig[inflPnts[troughICandidates2]-1])) .== minimum(abs(diff(sig[inflPnts[troughICandidates2]-1]))))[1] #need to check here which point to take
##                 idx = find(abs.(dy[inflPnts[troughICandidates2]]) .== minimum(abs.(dy[inflPnts[troughICandidates2]])))[1]
##                 troughIPnt = inflPnts[troughICandidates2[idx]]
##                 troughITime = inflTimes[troughICandidates2[idx]]
##             end

##         elseif nTroughICandidates == 1
##             troughIPnt = troughPnts[troughICandidates[1]]
##             troughITime = troughTimes[troughICandidates[1]]
##         else
##             idx = find(sig[troughPnts[troughICandidates]] .== minimum(sig[troughPnts[troughICandidates]]))[1]
##             troughIPnt = troughPnts[troughICandidates[idx]]
##             troughITime = troughTimes[troughICandidates[idx]]
##         end
##     end

##     #peak II
##     minPeakIILat = NaN
##     maxPeakIILat = NaN
##     if isnan(peakIIITime) == false
##         maxPeakIILat = peakIIITime - 0.25/1000
##     end
##     if isnan(troughITime) == false
##         minPeakIILat = troughITime + 0.25/1000
##     end

##     ## heuristic to try to find peak II when peak III is not found, this is experimental
##     #if ((isnan(maxPeakIILat) == true & isnan(peakVTime) == false) & (isnan(minPeakIILat) == false))
##     if isnan(maxPeakIILat) == true
##         if isnan(peakVTime) == false
##             if isnan(minPeakIILat) == false
##                 maxPeakIILat = (peakVTime + minPeakIILat)/2.5
##             end
##         end
##     end

##     if ((isnan(minPeakIILat) == false) .& (isnan(maxPeakIILat) == false))
##         peakIICandidates = find((peakTimes .< maxPeakIILat ) .& (peakTimes .> minPeakIILat))
##         nPeakIICandidates = length(peakIICandidates)
##         if nPeakIICandidates == 0
##             peakIICandidates2 = find((inflTimes .< maxPeakIILat ) .& (inflTimes .> minPeakIILat))
##             nTroughICandidates2 = length(peakIICandidates2)
##             if nTroughICandidates2 == 0
##                 peakIIPnt = NaN
##                 peakIITime = NaN
##             elseif nTroughICandidates2 == 1
##                 peakIIPnt = inflPnts[peakIICandidates2[1]]
##                 peakIITime = inflTimes[peakIICandidates2[1]]
##             else
##                 ##idx = find(sig[inflPnts[peakIICandidates2]] .== maximum(sig[inflPnts[peakIICandidates2]]))
##                 ##idx = find(sig[inflPnts[peakIICandidates2]] .== minimum(diff(sig[inflPnts[peakIICandidates2]-1]))) #need to check here which point to take
##                 idx = find(abs.(dy[inflPnts[peakIICandidates2]]) .== minimum(abs.(dy[inflPnts[peakIICandidates2]])))[1]
##                 peakIIPnt = inflPnts[peakIICandidates2[idx]]
##                 peakIITime = inflTimes[peakIICandidates2[idx]]
##             end
##         elseif nPeakIICandidates == 1
##             peakIIPnt = peakPnts[peakIICandidates[1]]
##             peakIITime = peakTimes[peakIICandidates[1]]
##         else
##             idx = find(sig[peakPnts[peakIICandidates]] .== maximum(sig[peakPnts[peakIICandidates]]))[1]
##             peakIIPnt = peakPnts[peakIICandidates[idx]]
##             peakIITime = peakTimes[peakIICandidates[idx]]
##         end
##     else
##         peakIIPnt = NaN
##         peakIITime = NaN
##     end

##     #trough II
##     minTroughIILat = NaN
##     maxTroughIILat = NaN
##     if isnan(peakIIITime) == false
##         maxTroughIILat = peakIIITime - 0.25/1000
##     end
##     if isnan(peakIITime) == false
##         minTroughIILat = peakIITime + minPeakIITroughIILat
##     end

##     if ((isnan(minTroughIILat) == false) .& (isnan(maxTroughIILat) == false))
##         troughIICandidates = find((troughTimes .< maxTroughIILat ) .& (troughTimes .> minTroughIILat))
##         troughIICandidates = troughIICandidates[find(sig[troughPnts[troughIICandidates]] .< sig[peakIIPnt])] #ensure trough amp is less than peak amp
##         nTroughIICandidates = length(troughIICandidates)
##         if nTroughIICandidates == 0
##             troughIICandidates2 = find((inflTimes .< maxTroughIILat ) .& (inflTimes .> minTroughIILat))
##             troughIICandidates2 = troughIICandidates2[find(sig[inflPnts[troughIICandidates2]] .< sig[peakIIPnt])] #ensure trough amp is less than peak amp
##             nTroughIICandidates2 = length(troughIICandidates2)
##             if nTroughIICandidates2 == 0
##                 troughIIPnt = NaN
##                 troughIITime = NaN
##             elseif nTroughIICandidates2 == 1
##                 troughIIPnt = inflPnts[troughIICandidates2[1]]
##                 troughIITime = inflTimes[troughIICandidates2[1]]
##             else
##                 ##idx = find(sig[inflPnts[troughIICandidates2]] .== minimum(sig[inflPnts[troughIICandidates2]])) #need to check here which point to take
##                 ##idx = find(sig[inflPnts[troughIICandidates2]] .== minimum(diff(sig[inflPnts[troughIICandidates2]-1]))) #need to check here which point to take
##                 idx = find(abs.(dy[inflPnts[troughIICandidates2]]) .== minimum(abs.(dy[inflPnts[troughIICandidates2]])))[1]
##                 troughIIPnt = inflPnts[troughIICandidates2[idx]]
##                 troughIITime = inflTimes[troughIICandidates2[idx]]
##             end

##         elseif nTroughIICandidates == 1
##             troughIIPnt = troughPnts[troughIICandidates[1]]
##             troughIITime = troughTimes[troughIICandidates[1]]
##         else
##             idx = find(sig[troughPnts[troughIICandidates]] .== minimum(sig[troughPnts[troughIICandidates]]))[1]
##             troughIIPnt = troughPnts[troughIICandidates[idx]]
##             troughIITime = troughTimes[troughIICandidates[idx]]
##         end
##     else
##         troughIIPnt = NaN
##         troughIITime = NaN
##     end


##     #peak IV
##     minPeakIVLat = NaN
##     maxPeakIVLat = NaN
##     if isnan(peakVTime) == false
##         maxPeakIVLat = peakVTime - 0.25/1000
##     end
##     if isnan(troughIIITime) == false
##         minPeakIVLat = troughIIITime + 0.25/1000
##     end

##     if ((isnan(minPeakIVLat) == false) .& (isnan(maxPeakIVLat) == false))
##         peakIVCandidates = find((peakTimes .< maxPeakIVLat ) .& (peakTimes .> minPeakIVLat))
##         nPeakIVCandidates = length(peakIVCandidates)
##         if nPeakIVCandidates == 0
##             peakIVCandidates2 = find((inflTimes .< maxPeakIVLat ) .& (inflTimes .> minPeakIVLat))
##             nPeakIVCandidates2 = length(peakIVCandidates2)
##             if nPeakIVCandidates2 == 0
##                 peakIVPnt = NaN
##                 peakIVTime = NaN
##             elseif nPeakIVCandidates2 == 1
##                 peakIVPnt = inflPnts[peakIVCandidates2[1]]
##                 peakIVTime = inflTimes[peakIVCandidates2[1]]
##             else
##                 ##idx = find(sig[inflPnts[peakIVCandidates2]] .== maximum(sig[inflPnts[peakIVCandidates2]]))
##                 ##idx = find(sig[inflPnts[peakIVCandidates2]] .== minimum(diff(sig[inflPnts[peakIVCandidates2]-1]))) #need to check here which point to take
##                 idx = find(abs.(dy[inflPnts[peakIVCandidates2]]) .== minimum(abs.(dy[inflPnts[peakIVCandidates2]])))[1]
##                 peakIVPnt = inflPnts[peakIVCandidates2[idx]]
##                 peakIVTime = inflTimes[peakIVCandidates2[idx]]
##             end
##         elseif nPeakIVCandidates == 1
##             peakIVPnt = peakPnts[peakIVCandidates[1]]
##             peakIVTime = peakTimes[peakIVCandidates[1]]
##         else
##             idx = find(sig[peakPnts[peakIVCandidates]] .== maximum(sig[peakPnts[peakIVCandidates]]))[1]
##             peakIVPnt = peakPnts[peakIVCandidates[idx]]
##             peakIVTime = peakTimes[peakIVCandidates[idx]]
##         end
##     else
##         peakIVPnt = NaN
##         peakIVTime = NaN
##     end

##     #trough IV
##     minTroughIVLat = NaN
##     maxTroughIVLat = NaN
##     if isnan(peakVTime) == false
##         maxTroughIVLat = peakVTime - 0.25/1000
##     end
##     if isnan(peakIVTime) == false
##         minTroughIVLat = peakIVTime + minPeakIVTroughIVLat
##     end

##     if ((isnan(minTroughIVLat) == false) .& (isnan(maxTroughIVLat) == false))
##         troughIVCandidates = find((troughTimes .< maxTroughIVLat ) .& (troughTimes .> minTroughIVLat))
##         troughIVCandidates = troughIVCandidates[find(sig[troughPnts[troughIVCandidates]] .< sig[peakIVPnt])] #ensure trough amp is less than peak amp
##         nTroughIVCandidates = length(troughIVCandidates)
##         if nTroughIVCandidates == 0
##             troughIVCandidates2 = find((inflTimes .< maxTroughIVLat ) .& (inflTimes .> minTroughIVLat))
##             troughIVCandidates2 = troughIVCandidates2[find(sig[inflPnts[troughIVCandidates2]] .< sig[peakIVPnt])] #ensure trough amp is less than peak amp
##             nTroughIVCandidates2 = length(troughIVCandidates2)
##             if nTroughIVCandidates2 == 0
##                 troughIVPnt = NaN
##                 troughIVTime = NaN
##             elseif nTroughIVCandidates2 == 1
##                 troughIVPnt = inflPnts[troughIVCandidates2[1]]
##                 troughIVTime = inflTimes[troughIVCandidates2[1]]
##             else
##                 ##idx = find(sig[inflPnts[troughIVCandidates2]] .== minimum(sig[inflPnts[troughIVCandidates2]])) #need to check here which point to take
##                 ##idx = find(sig[inflPnts[troughIVCandidates2]] .== minimum(diff(sig[inflPnts[troughIVCandidates2]-1]))) #need to check here which point to take
##                 idx = find(abs.(dy[inflPnts[troughIVCandidates2]]) .== minimum(abs.(dy[inflPnts[troughIVCandidates2]])))[1]
##                 troughIVPnt = inflPnts[troughIVCandidates2[idx]]
##                 troughIVTime = inflTimes[troughIVCandidates2[idx]]
##             end

##         elseif nTroughIVCandidates == 1
##             troughIVPnt = troughPnts[troughIVCandidates[1]]
##             troughIVTime = troughTimes[troughIVCandidates[1]]
##         else
##             idx = find(sig[troughPnts[troughIVCandidates]] .== minimum(sig[troughPnts[troughIVCandidates]]))[1]
##             troughIVPnt = troughPnts[troughIVCandidates[idx]]
##             troughIVTime = troughTimes[troughIVCandidates[idx]]
##         end
##     else
##         troughIVPnt = NaN
##         troughIVTime = NaN
##     end



##     peakPoints = [peakIPnt; peakIIPnt; peakIIIPnt; peakIVPnt; peakVPnt]
##     peakLatencies = [peakITime; peakIITime; peakIIITime; peakIVTime; peakVTime]
##     troughPoints = [troughIPnt; troughIIPnt; troughIIIPnt; troughIVPnt; troughVPnt]
##     troughLatencies = [troughITime; troughIITime; troughIIITime; troughIVTime; troughVTime]

##     peakAmps = zeros(length(peakPoints))
##     troughAmps = zeros(length(peakPoints))
##     for i=1:length(peakAmps)
##         if isnan(peakPoints[i]) == false
##             peakAmps[i] = sig[1, round(Int, peakPoints[i])]
##         else
##             peakAmps[i] = NaN
##         end

##         if isnan(troughPoints[i]) == false
##             troughAmps[i] = sig[1, round(Int, troughPoints[i])]
##         else
##             troughAmps[i] = NaN
##         end
##     end

##     peakTroughAmps = peakAmps-troughAmps

##     prms = Dict("minPeakLat" => [minPeakILat, minPeakIILat, minPeakIIILat, minPeakIVLat, minPeakVLat],
##                 "maxPeakLat" => [maxPeakILat, maxPeakIILat, maxPeakIIILat, maxPeakIVLat, maxPeakVLat],
##                 "minTroughLat" => [minTroughILat, minTroughIILat, minTroughIIILat, minTroughIVLat, minTroughVLat],
##                 "maxTroughLat" => [maxTroughILat, maxTroughIILat, maxTroughIIILat, maxTroughIVLat, maxTroughVLat])

##     return peakPoints, troughPoints, peakLatencies, troughLatencies, peakAmps, troughAmps, peakTroughAmps, prms
## end



