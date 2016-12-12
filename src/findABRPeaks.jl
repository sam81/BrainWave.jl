

"""
Attempt to find peaks and troughs of the ABR response for a click of a given level. The algorithm is
largely based on Bradley and Wilson (2004).

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the ABR waveform for which the peaks and troughs are sought. 
* `stimLevel::Real`: the level of the click used to evoke the ABR response.
* `sampRate::Real`: the sampling rate of the ABR recording.
* `epochStart::Real`: the time, in seconds, at which the epoch starts.
* `minInterPeakLatency::Real``: the minimum allowed latency between peaks.
* `dBRef::String`: whether the stimulus level is specified in dB SPL `SPL`, or in dB normal hearing level `nHL`.

##### Returns

* `peakPoints::`: array containing the points at which peaks I-V are detected
* `troughPoints`
* `peakLatencies`
* `troughLatencies`
* `peakAmps`
* `troughAmps`
* `peakTroughAmps`
* `prms`


##### References

* Bradley, A. P., & Wilson, W. J. (2004). Automated Analysis of the Auditory Brainstem Response. Proceedings of the 2004 Intelligent Sensors, Sensor Networks and Information Processing Conference, 2004., 541–546. http://doi.org/10.1109/ISSNIP.2004.1417519

##### Examples

```julia

```

"""
#min interpeak latency = 0.5 ms
#0 dB nHL = 31 dB ppe SPL
function findABRPeaks{T<:Real}(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, stimLevel::Real, sampRate::Real; epochStart::Real=0, minInterPeakLatency::Real=0.0005, dBRef::String="SPL")

    if dBRef == "SPL"
        stimLevel = stimLevel - 31
    elseif dBRef == "nHL"
        stimLevel = stimLevel
    else
        error("`stimLevel` needs to be specified in `SPL` or `nHL`")
    end

    #equation derived from Prosser, S., & Arslan, E. (1987). Prediction of auditory brainstem wave V latency as a diagnostic tool of sensorineural hearing loss. Audiology : Official Organ of the International Society of Audiology, 26(3), 179–87. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/3662941
    avPeakVLat = (-5.859e-08*stimLevel^4 + 1.274e-05*stimLevel^3 - 0.0005424*stimLevel^2 - 0.05297*stimLevel + 9.3)/1000
    minPeakVLat = avPeakVLat-0.001
    maxPeakVLat = avPeakVLat+0.001

    #todo enforce troughs to be lower in amplitude than peaks
    #shorten peak II trough II latency

    avPeakVpeakI_IPL = 0.0038
    avPeakVpeakIII_IPL = 0.0017

    minPeakVpeakI_IPL = avPeakVpeakI_IPL - 0.0002*2
    maxPeakVpeakI_IPL = avPeakVpeakI_IPL + 0.0002*2
    minPeakVpeakIII_IPL = avPeakVpeakIII_IPL - 0.00015*2
    maxPeakVpeakIII_IPL = avPeakVpeakIII_IPL + 0.00015*2

    minPeakITroughILat = 0.25/1000
    maxPeakITroughILat = minPeakITroughILat+0.75/1000

    minPeakIIITroughIIILat = 0.25/1000
    maxPeakIIITroughIIILat = minPeakIIITroughIIILat+0.75/1000

    minPeakVTroughVLat = 0.25/1000
    maxPeakVTroughVLat = minPeakVTroughVLat+1.2/1000

    minPeakIITroughIILat = 0.12/1000
    minPeakIVTroughIVLat = 0.12/1000


    peakPnts, peakTimes = findPeaks(sig, sampRate, epochStart=epochStart)
    troughPnts, troughTimes = findTroughs(sig, sampRate, epochStart=epochStart)
    inflPnts, inflTimes = findInflections(sig, sampRate, epochStart=epochStart)
    dy = zeros(size(sig)); dy[:,2:size(sig)[2]] = diff(sig, 2)
    ddy = zeros(size(sig)); ddy[:,2:size(sig)[2]] = diff(dy, 2)


    peakVCandidates = find((peakTimes .< maxPeakVLat ) & (peakTimes .> minPeakVLat))
    nPeakVCandidates = length(peakVCandidates)
    if nPeakVCandidates == 0
        peakVPnt = NaN
        peakVTime = NaN
    elseif nPeakVCandidates == 1
        peakVPnt = peakPnts[peakVCandidates[1]]
        peakVTime = peakTimes[peakVCandidates[1]]
    else
        idx = find(sig[peakPnts[peakVCandidates]] .== maximum(sig[peakPnts[peakVCandidates]]))[1]
        peakVPnt = peakPnts[peakVCandidates[idx]]
        peakVTime = peakTimes[peakVCandidates[idx]]
    end

    if isnan(peakVPnt) == true
        peakIPnt = NaN
        peakITime = NaN
        peakIIIPnt = NaN
        peakIIITime = NaN
        troughVPnt = NaN
        troughVTime = NaN
        minPeakILat = NaN
        maxPeakILat = NaN
        minPeakIIILat = NaN
        maxPeakIIILat = NaN
        minTroughVLat = NaN
        maxTroughVLat = NaN
    else
        minPeakILat = peakVTime - maxPeakVpeakI_IPL
        maxPeakILat = peakVTime - minPeakVpeakI_IPL
        minPeakIIILat = peakVTime - maxPeakVpeakIII_IPL
        maxPeakIIILat = peakVTime - minPeakVpeakIII_IPL

        #Peak I
        peakICandidates = find((peakTimes .< maxPeakILat ) & (peakTimes .> minPeakILat))
        nPeakICandidates = length(peakICandidates)
        if nPeakICandidates == 0
            peakIPnt = NaN
            peakITime = NaN
        elseif nPeakICandidates == 1
            peakIPnt = peakPnts[peakICandidates[1]]
            peakITime = peakTimes[peakICandidates[1]]
        else
            idx = find(sig[peakPnts[peakICandidates]] .== maximum(sig[peakPnts[peakICandidates]]))[1]
            peakIPnt = peakPnts[peakICandidates[idx]]
            peakITime = peakTimes[peakICandidates[idx]]
        end

        #Peak III
        peakIIICandidates = find((peakTimes .< maxPeakIIILat ) & (peakTimes .> minPeakIIILat))
        nPeakIIICandidates = length(peakIIICandidates)
        if nPeakIIICandidates == 0
            peakIIIPnt = NaN
            peakIIITime = NaN
        elseif nPeakIIICandidates == 1
            peakIIIPnt = peakPnts[peakIIICandidates[1]]
            peakIIITime = peakTimes[peakIIICandidates[1]]
        else
            idx = find(sig[peakPnts[peakIIICandidates]] .== maximum(sig[peakPnts[peakIIICandidates]]))[1]
            peakIIIPnt = peakPnts[peakIIICandidates[idx]]
            peakIIITime = peakTimes[peakIIICandidates[idx]]
        end

        #Trough V
        minTroughVLat = peakVTime + minPeakVTroughVLat
        maxTroughVLat = peakVTime + maxPeakVTroughVLat

        troughVCandidates = find((troughTimes .< maxTroughVLat ) & (troughTimes .> minTroughVLat))
        troughVCandidates = troughVCandidates[find(sig[troughPnts[troughVCandidates]] .< sig[peakVPnt])] #ensure trough amp is less than peak amp
        nTroughVCandidates = length(troughVCandidates)
        if nTroughVCandidates == 0
            troughVCandidates2 = find((inflTimes .< maxTroughVLat ) & (inflTimes .> minTroughVLat))
            troughVCandidates2 = troughVCandidates2[find(sig[inflPnts[troughVCandidates2]] .< sig[peakVPnt])] #ensure trough amp is less than peak amp
            nTroughVCandidates2 = length(troughVCandidates2)
            if nTroughVCandidates2 == 0
                troughVPnt = NaN
                troughVTime = NaN
            elseif nTroughVCandidates2 == 1
                troughVPnt = inflPnts[troughVCandidates2[1]]
                troughVTime = inflTimes[troughVCandidates2[1]]
            else
                ##idx = find(sig[inflPnts[troughVCandidates2]] .== minimum(sig[inflPnts[troughVCandidates2]])) #need to check here which point to take
                ##idx = find(abs(diff(sig[inflPnts[troughVCandidates2]-1])) .== minimum(abs(diff(sig[inflPnts[troughVCandidates2]-1])))) #need to check here which point to take
                idx = find(abs(dy[inflPnts[troughVCandidates2]]) .== minimum(abs(dy[inflPnts[troughVCandidates2]])))[1]
                troughVPnt = inflPnts[troughVCandidates2[idx]]
                troughVTime = inflTimes[troughVCandidates2[idx]]
            end
        elseif nTroughVCandidates == 1
            troughVPnt = troughPnts[troughVCandidates[1]]
            troughVTime = troughTimes[troughVCandidates[1]]
        else
            idx = find(sig[troughPnts[troughVCandidates]] .== minimum(sig[troughPnts[troughVCandidates]]))[1]
            troughVPnt = troughPnts[troughVCandidates[idx]]
            troughVTime = troughTimes[troughVCandidates[idx]]
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

        troughIIICandidates = find((troughTimes .< maxTroughIIILat ) & (troughTimes .> minTroughIIILat))
        troughIIICandidates = troughIIICandidates[find(sig[troughPnts[troughIIICandidates]] .< sig[peakIIIPnt])] #ensure trough amp is less than peak amp
        nTroughIIICandidates = length(troughIIICandidates)
        if nTroughIIICandidates == 0
            troughIIICandidates2 = find((inflTimes .< maxTroughIIILat ) & (inflTimes .> minTroughIIILat))
            troughIIICandidates2 = troughIIICandidates2[find(sig[inflPnts[troughIIICandidates2]] .< sig[peakIIIPnt])] #ensure trough amp is less than peak amp
            nTroughIIICandidates2 = length(troughIIICandidates2)
            if nTroughIIICandidates2 == 0
                troughIIIPnt = NaN
                troughIIITime = NaN
            elseif nTroughIIICandidates2 == 1
                troughIIIPnt = inflPnts[troughIIICandidates2[1]]
                troughIIITime = inflTimes[troughIIICandidates2[1]]
            else
                ##idx = find(sig[inflPnts[troughIIICandidates2]] .== minimum(sig[inflPnts[troughIIICandidates2]])) #need to check here which point to take
                ##idx = find(sig[inflPnts[troughIIICandidates2]] .== minimum(diff(sig[inflPnts[troughIIICandidates2]-1]))) #need to check here which point to take
                idx = find(abs(dy[inflPnts[troughIIICandidates2]]) .== minimum(abs(dy[inflPnts[troughIIICandidates2]])))[1]
                troughIIIPnt = inflPnts[troughIIICandidates2[idx]]
                troughIIITime = inflTimes[troughIIICandidates2[idx]]
            end
        elseif nTroughIIICandidates == 1
            troughIIIPnt = troughPnts[troughIIICandidates[1]]
            troughIIITime = troughTimes[troughIIICandidates[1]]
        else
            idx = find(sig[troughPnts[troughIIICandidates]] .== minimum(sig[troughPnts[troughIIICandidates]]))[1]
            troughIIIPnt = troughPnts[troughIIICandidates[idx]]
            troughIIITime = troughTimes[troughIIICandidates[idx]]
        end
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

        troughICandidates = find((troughTimes .< maxTroughILat ) & (troughTimes .> minTroughILat))
        troughICandidates = troughICandidates[find(sig[troughPnts[troughICandidates]] .< sig[peakIPnt])] #ensure trough amp is less than peak amp
        nTroughICandidates = length(troughICandidates)
        if nTroughICandidates == 0
            troughICandidates2 = find((inflTimes .< maxTroughILat ) & (inflTimes .> minTroughILat))
            troughICandidates2 = troughICandidates2[find(sig[inflPnts[troughICandidates2]] .< sig[peakIPnt])] #ensure trough amp is less than peak amp
            nTroughICandidates2 = length(troughICandidates2)
            if nTroughICandidates2 == 0
                troughIPnt = NaN
                troughITime = NaN
            elseif nTroughICandidates2 == 1
                troughIPnt = inflPnts[troughICandidates2[1]]
                troughITime = inflTimes[troughICandidates2[1]]
            else
                ##idx = find(sig[inflPnts[troughICandidates2]] .== minimum(sig[inflPnts[troughICandidates2]])) #need to check here which point to take
                ##idx = find(abs(diff(sig[inflPnts[troughICandidates2]-1])) .== minimum(abs(diff(sig[inflPnts[troughICandidates2]-1]))))[1] #need to check here which point to take
                idx = find(abs(dy[inflPnts[troughICandidates2]]) .== minimum(abs(dy[inflPnts[troughICandidates2]])))[1]
                troughIPnt = inflPnts[troughICandidates2[idx]]
                troughITime = inflTimes[troughICandidates2[idx]]
            end

        elseif nTroughICandidates == 1
            troughIPnt = troughPnts[troughICandidates[1]]
            troughITime = troughTimes[troughICandidates[1]]
        else
            idx = find(sig[troughPnts[troughICandidates]] .== minimum(sig[troughPnts[troughICandidates]]))[1]
            troughIPnt = troughPnts[troughICandidates[idx]]
            troughITime = troughTimes[troughICandidates[idx]]
        end
    end

    #peak II
    minPeakIILat = NaN
    maxPeakIILat = NaN
    if isnan(peakIIITime) == false
        maxPeakIILat = peakIIITime - 0.25/1000
    end
    if isnan(troughITime) == false
        minPeakIILat = troughITime + 0.25/1000
    end

    ## heuristic to try to find peak II when peak III is not found, this is experimental
    #if ((isnan(maxPeakIILat) == true & isnan(peakVTime) == false) & (isnan(minPeakIILat) == false))
    if isnan(maxPeakIILat) == true
        if isnan(peakVTime) == false
            if isnan(minPeakIILat) == false
                maxPeakIILat = (peakVTime + minPeakIILat)/2.5
            end
        end
    end

    if ((isnan(minPeakIILat) == false) & (isnan(maxPeakIILat) == false))
        peakIICandidates = find((peakTimes .< maxPeakIILat ) & (peakTimes .> minPeakIILat))
        nPeakIICandidates = length(peakIICandidates)
        if nPeakIICandidates == 0
            peakIICandidates2 = find((inflTimes .< maxPeakIILat ) & (inflTimes .> minPeakIILat))
            nTroughICandidates2 = length(peakIICandidates2)
            if nTroughICandidates2 == 0
                peakIIPnt = NaN
                peakIITime = NaN
            elseif nTroughICandidates2 == 1
                peakIIPnt = inflPnts[peakIICandidates2[1]]
                peakIITime = inflTimes[peakIICandidates2[1]]
            else
                ##idx = find(sig[inflPnts[peakIICandidates2]] .== maximum(sig[inflPnts[peakIICandidates2]]))
                ##idx = find(sig[inflPnts[peakIICandidates2]] .== minimum(diff(sig[inflPnts[peakIICandidates2]-1]))) #need to check here which point to take
                idx = find(abs(dy[inflPnts[peakIICandidates2]]) .== minimum(abs(dy[inflPnts[peakIICandidates2]])))[1]
                peakIIPnt = inflPnts[peakIICandidates2[idx]]
                peakIITime = inflTimes[peakIICandidates2[idx]]
            end
        elseif nPeakIICandidates == 1
            peakIIPnt = peakPnts[peakIICandidates[1]]
            peakIITime = peakTimes[peakIICandidates[1]]
        else
            idx = find(sig[peakPnts[peakIICandidates]] .== maximum(sig[peakPnts[peakIICandidates]]))[1]
            peakIIPnt = peakPnts[peakIICandidates[idx]]
            peakIITime = peakTimes[peakIICandidates[idx]]
        end
    else
        peakIIPnt = NaN
        peakIITime = NaN
    end

    #trough II
    minTroughIILat = NaN
    maxTroughIILat = NaN
    if isnan(peakIIITime) == false
        maxTroughIILat = peakIIITime - 0.25/1000
    end
    if isnan(peakIITime) == false
        minTroughIILat = peakIITime + minPeakIITroughIILat
    end

    if ((isnan(minTroughIILat) == false) & (isnan(maxTroughIILat) == false))
        troughIICandidates = find((troughTimes .< maxTroughIILat ) & (troughTimes .> minTroughIILat))
        troughIICandidates = troughIICandidates[find(sig[troughPnts[troughIICandidates]] .< sig[peakIIPnt])] #ensure trough amp is less than peak amp
        nTroughIICandidates = length(troughIICandidates)
        if nTroughIICandidates == 0
            troughIICandidates2 = find((inflTimes .< maxTroughIILat ) & (inflTimes .> minTroughIILat))
            troughIICandidates2 = troughIICandidates2[find(sig[inflPnts[troughIICandidates2]] .< sig[peakIIPnt])] #ensure trough amp is less than peak amp
            nTroughIICandidates2 = length(troughIICandidates2)
            if nTroughIICandidates2 == 0
                troughIIPnt = NaN
                troughIITime = NaN
            elseif nTroughIICandidates2 == 1
                troughIIPnt = inflPnts[troughIICandidates2[1]]
                troughIITime = inflTimes[troughIICandidates2[1]]
            else
                ##idx = find(sig[inflPnts[troughIICandidates2]] .== minimum(sig[inflPnts[troughIICandidates2]])) #need to check here which point to take
                ##idx = find(sig[inflPnts[troughIICandidates2]] .== minimum(diff(sig[inflPnts[troughIICandidates2]-1]))) #need to check here which point to take
                idx = find(abs(dy[inflPnts[troughIICandidates2]]) .== minimum(abs(dy[inflPnts[troughIICandidates2]])))[1]
                troughIIPnt = inflPnts[troughIICandidates2[idx]]
                troughIITime = inflTimes[troughIICandidates2[idx]]
            end

        elseif nTroughIICandidates == 1
            troughIIPnt = troughPnts[troughIICandidates[1]]
            troughIITime = troughTimes[troughIICandidates[1]]
        else
            idx = find(sig[troughPnts[troughIICandidates]] .== minimum(sig[troughPnts[troughIICandidates]]))[1]
            troughIIPnt = troughPnts[troughIICandidates[idx]]
            troughIITime = troughTimes[troughIICandidates[idx]]
        end
    else
        troughIIPnt = NaN
        troughIITime = NaN
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

    if ((isnan(minPeakIVLat) == false) & (isnan(maxPeakIVLat) == false))
        peakIVCandidates = find((peakTimes .< maxPeakIVLat ) & (peakTimes .> minPeakIVLat))
        nPeakIVCandidates = length(peakIVCandidates)
        if nPeakIVCandidates == 0
            peakIVCandidates2 = find((inflTimes .< maxPeakIVLat ) & (inflTimes .> minPeakIVLat))
            nPeakIVCandidates2 = length(peakIVCandidates2)
            if nPeakIVCandidates2 == 0
                peakIVPnt = NaN
                peakIVTime = NaN
            elseif nPeakIVCandidates2 == 1
                peakIVPnt = inflPnts[peakIVCandidates2[1]]
                peakIVTime = inflTimes[peakIVCandidates2[1]]
            else
                ##idx = find(sig[inflPnts[peakIVCandidates2]] .== maximum(sig[inflPnts[peakIVCandidates2]]))
                ##idx = find(sig[inflPnts[peakIVCandidates2]] .== minimum(diff(sig[inflPnts[peakIVCandidates2]-1]))) #need to check here which point to take
                idx = find(abs(dy[inflPnts[peakIVCandidates2]]) .== minimum(abs(dy[inflPnts[peakIVCandidates2]])))[1]
                peakIVPnt = inflPnts[peakIVCandidates2[idx]]
                peakIVTime = inflTimes[peakIVCandidates2[idx]]
            end
        elseif nPeakIVCandidates == 1
            peakIVPnt = peakPnts[peakIVCandidates[1]]
            peakIVTime = peakTimes[peakIVCandidates[1]]
        else
            idx = find(sig[peakPnts[peakIVCandidates]] .== maximum(sig[peakPnts[peakIVCandidates]]))[1]
            peakIVPnt = peakPnts[peakIVCandidates[idx]]
            peakIVTime = peakTimes[peakIVCandidates[idx]]
        end
    else
        peakIVPnt = NaN
        peakIVTime = NaN
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

    if ((isnan(minTroughIVLat) == false) & (isnan(maxTroughIVLat) == false))
        troughIVCandidates = find((troughTimes .< maxTroughIVLat ) & (troughTimes .> minTroughIVLat))
        troughIVCandidates = troughIVCandidates[find(sig[troughPnts[troughIVCandidates]] .< sig[peakIVPnt])] #ensure trough amp is less than peak amp
        nTroughIVCandidates = length(troughIVCandidates)
        if nTroughIVCandidates == 0
            troughIVCandidates2 = find((inflTimes .< maxTroughIVLat ) & (inflTimes .> minTroughIVLat))
            troughIVCandidates2 = troughIVCandidates2[find(sig[inflPnts[troughIVCandidates2]] .< sig[peakIVPnt])] #ensure trough amp is less than peak amp
            nTroughIVCandidates2 = length(troughIVCandidates2)
            if nTroughIVCandidates2 == 0
                troughIVPnt = NaN
                troughIVTime = NaN
            elseif nTroughIVCandidates2 == 1
                troughIVPnt = inflPnts[troughIVCandidates2[1]]
                troughIVTime = inflTimes[troughIVCandidates2[1]]
            else
                ##idx = find(sig[inflPnts[troughIVCandidates2]] .== minimum(sig[inflPnts[troughIVCandidates2]])) #need to check here which point to take
                ##idx = find(sig[inflPnts[troughIVCandidates2]] .== minimum(diff(sig[inflPnts[troughIVCandidates2]-1]))) #need to check here which point to take
                idx = find(abs(dy[inflPnts[troughIVCandidates2]]) .== minimum(abs(dy[inflPnts[troughIVCandidates2]])))[1]
                troughIVPnt = inflPnts[troughIVCandidates2[idx]]
                troughIVTime = inflTimes[troughIVCandidates2[idx]]
            end

        elseif nTroughIVCandidates == 1
            troughIVPnt = troughPnts[troughIVCandidates[1]]
            troughIVTime = troughTimes[troughIVCandidates[1]]
        else
            idx = find(sig[troughPnts[troughIVCandidates]] .== minimum(sig[troughPnts[troughIVCandidates]]))[1]
            troughIVPnt = troughPnts[troughIVCandidates[idx]]
            troughIVTime = troughTimes[troughIVCandidates[idx]]
        end
    else
        troughIVPnt = NaN
        troughIVTime = NaN
    end



    peakPoints = [peakIPnt; peakIIPnt; peakIIIPnt; peakIVPnt; peakVPnt]
    peakLatencies = [peakITime; peakIITime; peakIIITime; peakIVTime; peakVTime]
    troughPoints = [troughIPnt; troughIIPnt; troughIIIPnt; troughIVPnt; troughVPnt]
    troughLatencies = [troughITime; troughIITime; troughIIITime; troughIVTime; troughVTime]

    peakAmps = zeros(length(peakPoints))
    troughAmps = zeros(length(peakPoints))
    for i=1:length(peakAmps)
        if isnan(peakPoints[i]) == false
            peakAmps[i] = sig[1, round(Int, peakPoints[i])]
        else
            peakAmps[i] = NaN
        end

        if isnan(troughPoints[i]) == false
            troughAmps[i] = sig[1, round(Int, troughPoints[i])]
        else
            troughAmps[i] = NaN
        end
    end

    peakTroughAmps = peakAmps-troughAmps

    prms = Dict("minPeakLat" => [minPeakILat, minPeakIILat, minPeakIIILat, minPeakIVLat, minPeakVLat],
                "maxPeakLat" => [maxPeakILat, maxPeakIILat, maxPeakIIILat, maxPeakIVLat, maxPeakVLat],
                "minTroughLat" => [minTroughILat, minTroughIILat, minTroughIIILat, minTroughIVLat, minTroughVLat],
                "maxTroughLat" => [maxTroughILat, maxTroughIILat, maxTroughIIILat, maxTroughIVLat, maxTroughVLat])

    return peakPoints, troughPoints, peakLatencies, troughLatencies, peakAmps, troughAmps, peakTroughAmps, prms
end

