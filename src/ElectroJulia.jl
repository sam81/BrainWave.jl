module ElectroJulia

export averageAverages, averageEpochs, baselineCorrect, chainSegments, deleteSlice3D, filterContinuous, _centered, fftconvolve, findArtefactThresh, getFRatios, getNoiseSidebands, getSpectrum, mergeEventTableCodes, nextPowTwo, removeEpochs, removeSpuriousTriggers, rerefCnt, segment
#segment
using DataFrames
using Distributions
using PyCall
#using Devectorize
#pyinitialize("python3")
@pyimport scipy.signal as scisig
function averageAverages(aveList, nSegments)
    ## """
    ## Perform a weighted average of a list of averages. The weight of
    ## each average in the list is determined by the number of segments
    ## from which it was obtained.
    
    ## Parameters
    ## ----------
    ## aveList : dict of list of 2D numpy arrays
    ##     The list of averages for each experimental condition.
    ## nSegments : dict of ints
    ##     The number of epochs on which each average is based 

    ## Returns
    ## ----------

    ## Examples
    ## ----------
    ## """
    eventList = collect(keys(aveList[1]))
    weightedAve = (String => Array{eltype(aveList[1][eventList[1]]),2})[]
    nSegsSum = (String => Int)[]
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

function averageEpochs(rec)
    ## """
    ## Average the epochs of a segmented recording.

    ## Parameters
    ## ----------
    ## rec : dict of 3D numpy arrays with dimensions (n_channels x n_samples x n_epochs)
    ##     Recording

    ## Returns
    ## ----------
    ## ave : dict of 2D numpy arrays with dimensions (n_channels x n_samples)
    ##     The average epochs for each condition.
    ## n_segs : dict of ints
    ##     The number of epochs averaged for each condition.
        
    ## Examples
    ## ----------
    ## >>> ave, n_segs = average_epochs(rec=rec)
    ## """
    
    eventList = collect(keys(rec))
    ave = (String => Array{eltype(rec[eventList[1]]),2})[]
    nSegs = (String => Int)[]
    for i=1:length(eventList)
        nSegs[eventList[i]] = size(rec[eventList[i]])[3]
        ave[eventList[i]] = mean(rec[eventList[i]], 3)[:,:,1]
    end
    return ave, nSegs
end

function baselineCorrect(rec, baselineStart::Real, preDur::Real, sampRate::Integer)
    
    ## Perform baseline correction by subtracting the average pre-event
    ## voltage from each channel of a segmented recording.

    ## Parameters
    ## ----------
    ## rec: dict of 3D arrays
    ##     The segmented recording.
    ## baseline_start: float
    ##     Start time of the baseline window relative to the event onset, in seconds.
    ##     The absolute value of baseline_start cannot be greater than pre_dur.
    ##     In practice baseline_start allows you to define a baseline window shorter
    ##     than the time window before the experimental event (pre_dur).
    ## pre_dur: float
    ##     Duration of recording before the experimental event, in seconds.
    ## samp_rate: int
    ##     The samplig rate of the EEG recording.
    
    ## Examples
    ## ----------
    ## #baseline window has the same duration of pre_dur
    ## >>> baseline_correct(rec=rec, baseline_start=-0.2, pre_dur=0.2, samp_rate=512)
    ## #now with a baseline shorter than pre_dur
    ## >>> baseline_correct(rec=rec, baseline_start=-0.15, pre_dur=0.2, samp_rate=512)

    eventList = collect(keys(rec))
    epochStartSample = int(round(preDur*sampRate))
    baselineStartSample = int((epochStartSample+1) - abs(round(baselineStart*sampRate)))
    
    for i=1:length(eventList) #for each event
        for j=1:size(rec[eventList[i]])[3] #for each epoch
            for k=1: size(rec[eventList[i]])[1] #for each electrode
                thisBaseline = mean(rec[eventList[i]][k,baselineStartSample:epochStartSample,j])
                rec[eventList[i]][k,:,j] = rec[eventList[i]][k,:,j] - thisBaseline
            end
        end
    end

    
end

function baselineCorrectloop(rec, baselineStart::Real, preDur::Real, sampRate::Integer)
    
    ## Perform baseline correction by subtracting the average pre-event
    ## voltage from each channel of a segmented recording.

    ## Parameters
    ## ----------
    ## rec: dict of 3D arrays
    ##     The segmented recording.
    ## baseline_start: float
    ##     Start time of the baseline window relative to the event onset, in seconds.
    ##     The absolute value of baseline_start cannot be greater than pre_dur.
    ##     In practice baseline_start allows you to define a baseline window shorter
    ##     than the time window before the experimental event (pre_dur).
    ## pre_dur: float
    ##     Duration of recording before the experimental event, in seconds.
    ## samp_rate: int
    ##     The samplig rate of the EEG recording.
    
    ## Examples
    ## ----------
    ## #baseline window has the same duration of pre_dur
    ## >>> baseline_correct(rec=rec, baseline_start=-0.2, pre_dur=0.2, samp_rate=512)
    ## #now with a baseline shorter than pre_dur
    ## >>> baseline_correct(rec=rec, baseline_start=-0.15, pre_dur=0.2, samp_rate=512)

    eventList = collect(keys(rec))
    epochStartSample = int(round(preDur*sampRate))
    baselineStartSample = int((epochStartSample+1) - abs(round(baselineStart*sampRate)))
 
    
    for i=1:length(eventList) #for each event
        for j=1:size(rec[eventList[i]])[3] #for each epoch
            for k=1: size(rec[eventList[i]])[1] #for each electrode
                thisBaseline = mean(rec[eventList[i]][k,baselineStartSample:epochStartSample,j])
                for s=1:size(rec[eventList[i]])[2] #for each sample     
                    rec[eventList[i]][k,s,j] = rec[eventList[i]][k,s,j] - thisBaseline
                end
            end
        end
    end
end

function chainSegments(rec, nChunks::Integer, sampRate::Integer, startTime::Real, endTime::Real, baselineDur::Real, window)
    ## """
    ## Take a dictionary containing in each key a list of segments, and chain these segments
    ## into chunks of length nChunks
    ## baselineDur is for determining what is the zero point
    ## startTime and endTime are given with reference to the zero point
    ## """
    baselinePnts = round(baselineDur * sampRate)
    startPnt = int(round(startTime*sampRate) + baselinePnts) +1
    endPnt = int(round(endTime*sampRate) + baselinePnts) 
    chunkSize = ((endPnt - startPnt)+1)
    sweepSize = chunkSize * nChunks

    if window != None
        n = chunkSize
        w = zeros(1, n)
        if window == "hamming"
            w[1,:] = scisig.hamming(n)
        elseif window == "hanning"
            w[1,:] = scisig.hanning(n)
        elseif window == "blackman"
            w[1,:] = scisig.blackman(n)
        elseif window == "bartlett"
            w[1,:] = scisig.bartlett(n)
        end
    end
    
    nReps = (String => Array{Int,1})[]
    eventList = collect(keys(rec))
    eegChained = (String => Array{eltype(rec[eventList[1]]),2})[]
    fromeegChainedAve = (String => Array{eltype(rec[eventList[1]]),2})[]
    for i=1:length(eventList)
        currCode = eventList[i]
        eegChained[currCode] = zeros(eltype(rec[eventList[1]]), size(rec[currCode])[1], sweepSize)  #two-dimensional array of zeros
        #fromeegChainedAve[currCode] = zeros(size(rec[currCode])[1], chunkSize)
        nReps[currCode] = zeros(Int, nChunks)
        p = 1
        k = 1
        while k <= size(rec[currCode])[3]
            if p > (nChunks)
                p = 1
            end
            
            idxChunkStart = ((p-1)*chunkSize)+1
            idxChunkEnd = (idxChunkStart + chunkSize)-1
            eegChained[currCode][:,idxChunkStart:idxChunkEnd] = eegChained[currCode][:,idxChunkStart:idxChunkEnd] + rec[currCode][:,startPnt:endPnt, k]
            nReps[currCode][p] = nReps[currCode][p] + 1
            #fromeegChainedAve[currCode] = fromeegChainedAve[currCode] + rec[currCode][:,startPnt:endPnt, k]
            p = p+1 #p is the chunk counter
            k = k+1 #k is the epoch counter
        end
    end
    for i=1:length(eventList)
        currCode = eventList[i]
        for p=1:nChunks
            idxChunkStart = ((p-1)*chunkSize)+1
            idxChunkEnd = (idxChunkStart + chunkSize)-1
            if window == None
                eegChained[currCode][:,idxChunkStart:idxChunkEnd] = eegChained[currCode][:,idxChunkStart:idxChunkEnd] / nReps[currCode][p]
            else
                for chn=1:size(eegChained[currCode])[1]
                    eegChained[currCode][chn,idxChunkStart:idxChunkEnd] = eegChained[currCode][chn,idxChunkStart:idxChunkEnd].*w / nReps[currCode][p]
               
                end
            end
        end
        #fromeegChainedAve[currCode] = fromeegChainedAve[currCode] / sum(nReps[currCode])
    end
        
    return eegChained
end

## function deleteSlice3D(x, toRemove)
##     y = similar(x, size(x)[1], size(x)[2], size(x)[3]-length(toRemove))
##     idx = 1
##     for i=1:size(x)[3]
##         if ~contains(toRemove, i)
##             y[:,:,idx] = x[:,:,i]
##             idx = idx+1
##         end
##     end
##     return y
## end

function deleteSlice3D(x, toRemove, axis)
    toKeep = (Int)[]
    for i=1:size(x)[axis]
        if ~in(i, toRemove)
            push!(toKeep, i)
        end
    end
    if axis == 1
        y = x[toKeep,:,:]
    elseif axis == 2
        y = x[:,toKeep,:]
    elseif axis == 3
        y = x[:,:,toKeep]
    end
    
    return y
end

## function deleteSliceND(x, toRemove, axis)
##     toKeep = (Integer)[]
##     for i=1:size(x)[axis]
##         if ~contains(toRemove, i)
##             push!(toKeep, i)
##         end
##     end
##     idx = (Any)[]
##     for i=1:length(size(x))
##         if i == axis
##             push!(idx, toKeep)
##         else
##             push!(idx, [1:size(x)[i]])
##         end
##     end
                
##     y = x[
##     y = similar(x, size(x)[1], size(x)[2], size(x)[3]-length(toRemove))
##     idx = 1
##     for i=1:size(x)[3]
##         if ~contains(toRemove, i)
##             y[:,:,idx] = x[:,:,i]
##             idx = idx+1
##         end
##     end
##     return y
## end

function filterContinuous(rec, channels, sampRate, filterType::String, nTaps::Integer, cutoffs, transitionWidth::Real)
    ## """
    
    ## Parameters
    ## ----------

    ## Returns
    ## ----------

    ## Examples
    ## ----------
    ## """
       
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

    b = convert(Array{eltype(rec),1}, scisig.firwin2(nTaps,f,m))
    #b = ones(Float32,nTaps)
    nChannels = size(rec)[1]
    if channels == None
        channels = [1:nChannels]
    end
   
    for i=1:nChannels
        if in(i, channels) == true
            rec[i,:] = fftconvolve(reshape(rec[i,:], size(rec[i,:])[2]), b, "same")
            rec[i,:] = flipdim(fftconvolve(flipdim(reshape(rec[i,:], size(rec[i,:])[2]),1), b, "same"), 1)
        end
    end
    return rec
end

function _centered(arr, newsize)
    # Return the center newsize portion of the array.
    currsize = size(arr)[1]
    startind = div((currsize - newsize), 2) + 1
    endind = startind + newsize -1 #check indexing is the same pyth julia?
    return arr[startind:endind]
end

function fftconvolve(x, y, mode)
    s1 = size(x)[1]#check if array has two dim?
    s2 = size(y)[1]
    
    convArray = conv(x,y)
    if mode == "full"
        return convArray
    elseif mode == "same"
        return _centered(convArray, s1)
    elseif mode == "valid"
        return _centered(convArray, abs(s1 - s2) + 1)
    end
end
    

function findArtefactThresh(rec, thresh, channels)
    ## """
    ## Find epochs with voltage values exceeding a given threshold.
    
    ## Parameters
    ## ----------
    ## rec : dict of 3D arrays
    ##     The segmented recording
    ## thresh :
    ##     The threshold value.
    ## channels = array or list of ints
    ##     The indexes of the channels on which to find artefacts.
        
    ## Returns
    ## ----------

    ## Examples
    ## ----------
    ## """
    if length(channels) != length(thresh)
        print("The number of thresholds must be equal to the number of channels \n")
        return
    end
    eventList = collect(keys(rec))
    segsToReject = (String => Array{Int,1})[]
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

function getFRatios(ffts, freqs, nSideComp, nExcludedComp, otherExclude)
    ##"""
    ##
    ##"""

    cnds = collect(keys(ffts))
    compIdx = (Int)[]
    for freq in freqs
        thisIdx = find(abs(ffts[cnds[1]]["freq"] - freq) .== minimum(abs(ffts[cnds[1]]["freq"] - freq)))
        append!(compIdx, thisIdx)
    end
    sideBandsIdx = (Int)[]
    idxProtect = (Int)[]
    fftVals = (String => Any)[]
    fRatio = (String => Any)[]
    dfNum = 2
    dfDenom = 2*(nSideComp*2) -1
    for cnd in cnds
        fRatio[cnd] = (String => Array{Float64, 1})[]
        fftVals[cnd] = (String => Array{Float64, 1})[]
        fRatio[cnd]["F"] = []
        fRatio[cnd]["pval"] = []
        fftVals[cnd]["sigPow"] = []
        fftVals[cnd]["noisePow"] = []
        sideBands, sideBandsIdx, idxProtect = getNoiseSidebands(freqs, nSideComp, nExcludedComp, ffts[cnd], otherExclude)
        for c=1:length(compIdx)
            noisePow = mean(sideBands[c])
            sigPow = ffts[cnd]["mag"][compIdx[c]]
            thisF =  sigPow/ noisePow
            fftVals[cnd]["sigPow"] = vcat(fftVals[cnd]["sigPow"], sigPow)
            fftVals[cnd]["noisePow"] = vcat(fftVals[cnd]["noisePow"], noisePow)
            fRatio[cnd]["F"] = vcat(fRatio[cnd]["F"], thisF)
            fRatio[cnd]["pval"] = vcat(fRatio[cnd]["pval"], pdf(FDist(dfNum, dfDenom), thisF))
        end
    end
    minSideFreq = (FloatingPoint)[]
    maxSideFreq = (FloatingPoint)[]
    for c=1:length(compIdx)
        push!(minSideFreq, ffts[cnds[1]]["freq"][minimum(sideBandsIdx[c])])
        push!(maxSideFreq, ffts[cnds[1]]["freq"][maximum(sideBandsIdx[c])])
    end
    res = (String => Any)["fftVals" => fftVals,
                          "fRatio" => fRatio,
                          "compIdx" => compIdx,
                          "sideBandsIdx" => sideBandsIdx,
                          "excludedIdx"  => idxProtect,
                          "minSideFreq" => minSideFreq,
                          "maxSideFreq" => maxSideFreq]
    return res
end

function getNoiseSidebands(freqs, nCompSide, nExcludedComp, fftDict, otherExclude)
    #"""
    #the 2 has the possibility to exclude extra components, useful for distortion products
    #"""
    #components: a list containing the indexes of the target components
    #nCompSide: number of components used for each side band
    #n_exclude_side: number of components adjacent to to the target components to exclude
    #fft_array: array containing the fft values
    compIdx = (Int)[]
    for freq in freqs
        thisIdx = find(abs(fftDict["freq"] - freq) .== minimum(abs(fftDict["freq"] - freq)))
        append!(compIdx, thisIdx)
    end
    
    idxProtect = []; idxProtect = vcat(idxProtect, compIdx)
    if otherExclude != None
        otherExcludeIdx = (Int)[]
        for i=1:length(otherExclude)
            append!(otherExcludeIdx, find(abs(fftDict["freq"] - otherExclude[i]) .== minimum(abs(fftDict["freq"] - otherExclude[i]))))
        end
        idxProtect = vcat(idxProtect, otherExcludeIdx)
    end
    
    for i=1:nExcludedComp
        idxProtect = vcat(idxProtect, compIdx + i)
        idxProtect = vcat(idxProtect, compIdx - i)
        for j=1:length(otherExclude)
            push!(idxProtect, otherExcludeIdx[j] - i)
            push!(idxProtect, otherExcludeIdx[j] + i)
        end
    end
    idxProtect = sort(idxProtect)

    noiseBands = (Any)[]
    noiseBandsIdx = (Any)[]
    for i=1:length(compIdx)
        loSide = []; hiSide = []
        loSideIdx = (Int)[]; hiSideIdx = (Int)[]
        counter = 1
        while length(hiSide) < nCompSide
            currIdx = compIdx[i] + nExcludedComp + counter
            if in(currIdx, idxProtect) == false
                hiSide = vcat(hiSide, fftDict["mag"][currIdx])
                push!(hiSideIdx, currIdx)
            end
            counter = counter + 1
        end
        counter = 1
        while length(loSide) < nCompSide
            currIdx = compIdx[i] - nExcludedComp - counter
            if in(currIdx, idxProtect) == false
                loSide = vcat(loSide, fftDict["mag"][currIdx])
                push!(loSideIdx, currIdx)
            end
            counter = counter + 1
        end
        push!(noiseBands, vcat(loSide, hiSide))
        push!(noiseBandsIdx, vcat(loSideIdx, hiSideIdx))
        #noiseBands = vcat(noiseBands, loSide+hiSide)
    end
    return noiseBands, noiseBandsIdx, idxProtect
end
                              

function getSpectrum(sig, sampRate::Integer, window::String, powerOfTwo::Bool)
    ## """
    
    ## Parameters
    ## ----------

    ## Returns
    ## ----------

    ## Examples
    ## ----------
    ## """
    n = length(sig)
    if powerOfTwo == true
        nfft = 2^nextPowTwo(n)
    else
        nfft = n
    end
    if window != "none"
        if window == "hamming"
             w = scisig.hamming(n)
        elseif window == "hanning"
             w = scisig.hanning(n)
        elseif window == "blackman"
             w = scisig.blackman(n)
        elseif window == "bartlett"
             w = scisig.bartlett(n)
        end
        sig = sig*w
    end
    p = fft(sig)#, nfft) # take the fourier transform
    
    nUniquePts = ceil((nfft+1)/2)
    p = p[1:nUniquePts]
    p = abs(p)
    p = p / sampRate  # scale by the number of points so that
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

    freq_array = [0:(nUniquePts-1)] * (sampRate / nfft)
    x = (String => Array{Float64,1})[]
    x["freq"] = freq_array; x["mag"] = p
    return x
end


function mergeEventTableCodes(eventTable::Dict{String,Any}, trigList, newTrig::Integer)
    ## """
    ## Substitute the event table triggers listed in trig_list
    ## with new_trig

    ## Parameters
    ## ----------
    ## event_table : dict of int arrays
    ##     The event table
    ## trig_list : array of ints
    ##     The list of triggers to substitute
    ## new_trig : int
    ##     The new trigger used to substitute the triggers
    ##     in trig_list
    ## Returns
    ## ----------

    ## Examples
    ## ----------
    ## """
    ## for i=1:length(trigList)
    ##     eventTable["code"][eventTable["code"] .== trigList[i]] = newTrig
    ## end
    eventTable["code"][findin(eventTable["code"], trigList)] = newTrig
    return
end

function nextPowTwo(x::Real)
    ## """
    
    ## Parameters
    ## ----------

    ## Returns
    ## ----------

    ## Examples
    ## ----------
    ## """
    out = int(ceil(log2(x)))
    return out
end

function removeEpochs(rec, toRemove)
    ## """
    ## Remove epochs from a segmented recording.
    
    ## Parameters
    ## ----------
    ## rec : dict of 3D arrays
    ##     The segmented recording
    ## to_remove : dict of 1D arrays
    ##     List of epochs to remove for each condition

    ## Returns
    ## ----------

    ## Examples
    ## ----------
    ## """
    eventList = collect(keys(rec))
    for i=1:length(eventList)
        code = eventList[i]
        rec[code] = deleteSlice3D(rec[code], toRemove[code], 3)
    end
    return rec
end

function removeSpuriousTriggers(eventTable::Dict{String, Any}, sentTrigs::Array{Int}, minTrigDur::Real)
    recTrigs = eventTable["code"]
    recTrigsStart = eventTable["idx"]
    recTrigsDur = eventTable["dur"]

    allowedTrigs = int16(unique(sentTrigs))
    allowedIdx = findin(recTrigs, allowedTrigs)
 
    recTrigsDur = recTrigsDur[allowedIdx]
    recTrigsStart = recTrigsStart[allowedIdx]
    recTrigs = recTrigs[allowedIdx]

    durCondition = recTrigsDur .>= minTrigDur
    recTrigs = recTrigs[durCondition]
    recTrigsStart = recTrigsStart[durCondition]
    recTrigsDur = recTrigsDur[durCondition]

    foo = zeros(Int16, size(sentTrigs)[1])
    for i=1:size(sentTrigs)[1]
        foo[i] = int16(sentTrigs[i])
    end

    if recTrigs == foo
        match_found = true
    else
        match_found = false
    end


    eventTable["code"] = recTrigs
    eventTable["dur"] = recTrigsDur
    eventTable["idx"] = recTrigsStart

    resInfo = (String => Any)[]
    resInfo["match"] = match_found
    resInfo["lenSent"] = length(sentTrigs)
    resInfo["lenFound"] = length(recTrigs)

    return resInfo
end

function rerefCnt(rec, refChan::Integer, channels)
    ## """
    ## Rereference channels in a continuous recording.

    ## Parameters
    ## ----------
    ## rec : 
    ##     Recording
    ## ref_channel: int
    ##     The reference channel (indexing starts from zero).
    ## channels : list of ints
    ##     List of channels to be rereferenced (indexing starts from zero).
  
    ## Returns
    ## -------
    ## rec : an array of floats with dimenions nChannels X nDataPoints
        
    ## Examples
    ## --------
    ## >>> reref_cnt(rec=dats, channels=[1, 2, 3], ref_channel=4)
    ## """

    if channels == None
        nChannels = size(rec)[1]
        channels = [1:nChannels]
    end

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



function segment(rec, eventTable::Dict{String, Any}, epochStart::Real, epochEnd::Real, sampRate::Integer, eventsList=None, eventsLabelsList=None)

    trigs = eventTable["code"]
    trigs_pos = eventTable["idx"]
    if eventsList == None
        eventsList = unique(trigs)
    end

    if eventsLabelsList == None
        eventsLabelsList = String[string(eventsList[i]) for i=1:length(eventsList)]
        end
        

    epochStartSample = int(round(epochStart*sampRate))
    epochEndSample = int(round(epochEnd*sampRate)) - 1

    nSamples = epochEndSample - epochStartSample + 1
    segs = (String => Array{eltype(rec),3})[]
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
                if thisStopPnt > size(rec)[2]#rec.shape[1]
                    print(idx[j], " Epoch ends after end of recording. \n")
                end
            
            else
                segs[eventsLabelsList[i]][:,:,j] = rec[:, thisStartPnt:thisStopPnt]
            
            end
        end
    end
    nSegs = (String => Int)[]
    for i=1:length(eventsList) #count
        nSegs[eventsLabelsList[i]] = size(segs[eventsLabelsList[i]])[3]
    end

    return segs, nSegs
        
end




end #Module
