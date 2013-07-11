#module ElectroJulia

#export segment
using DataFrames
using Distributions
using PyCall
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
    weightedAve = (UTF8String => Array{eltype(aveList[1][eventList[1]]),2})[]
    nSegsSum = (UTF8String => Int)[]
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
    ave = (UTF8String => Array{eltype(rec[eventList[1]]),2})[]
    nSegs = (UTF8String => Int)[]
    for i=1:length(eventList)
        nSegs[eventList[i]] = size(rec[eventList[i]])[3]
        ave[eventList[i]] = mean(rec[eventList[i]], 3)[:,:,1]
    end
    return ave, nSegs
end

function baselineCorrect(rec, baselineStart::Real, preDur::Real, sampRate::Int)
    
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

function chainSegments(rec, nChunks::Int, sampRate::Int, startTime::Real, endTime::Real, baselineDur::Real)
    """
    Take a dictionary containing in each key a list of segments, and chain these segments
    into chunks of length nChunks
    baselineDur is for determining what is the zero point
    startTime and endTime are given with reference to the zero point
    """
    baselinePnts = round(baselineDur * sampRate)
    startPnt = int(round(startTime*sampRate) + baselinePnts) 
    endPnt = int(round(endTime*sampRate) + baselinePnts) - 1
    chunkSize = ((endPnt - startPnt)+1)
    sweepSize = chunkSize * nChunks
    nReps = (UTF8String => Array{Int,1})[]
    eventList = collect(keys(rec))
    eegChained = (UTF8String => Array{eltype(rec[eventList[1]]),2})[]
    fromeegChainedAve = (UTF8String => Array{eltype(rec[eventList[1]]),2})[]
    for i=1:length(eventList)
        currCode = eventList[i]
        eegChained[currCode] = zeros(eltype(rec[eventList[1]]), size(rec[currCode])[1], sweepSize)  #two-dimensional array of zeros
        #fromeegChainedAve[currCode] = zeros(size(rec[currCode])[1], chunkSize)
        nReps[currCode] = zeros(Int, nChunks)
        p = 1
        k = 1
        while k < size(rec[currCode])[3]
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
            eegChained[currCode][:,idxChunkStart:idxChunkEnd] = eegChained[currCode][:,idxChunkStart:idxChunkEnd] / nReps[currCode][p]
        end
        #fromeegChainedAve[currCode] = fromeegChainedAve[currCode] / sum(nReps[currCode])
    end
        
    return eegChained
end

function deleteSlice3D(x, toRemove)
    y = similar(x, size(x)[1], size(x)[2], size(x)[3]-length(toRemove))
    idx = 1
    for i=1:size(x)[3]
        if ~contains(toRemove, i)
            y[:,:,idx] = x[:,:,i]
            idx = idx+1
            end
    end
    return y
end

function filterContinuous(rec, channels, sampRate, filterType::String, nTaps::Int, cutoffs, transitionWidth::Real)
    ## """
    
    ## Parameters
    ## ----------

    ## Returns
    ## ----------

    ## Examples
    ## ----------
    ## """
       
    if filterType == "lowpass"
        f1 = cutoffs[1] * (1-transitionWidth)
        f2 = cutoffs[1]
        f1 = (f1*2) / sampRate
        f2 = (f2*2) / sampRate
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

    b = float32(scisig.firwin2(nTaps,f,m))
    #println(typeof(rec[1,:]))
    #println(typeof(b))
    nChannels = size(rec)[1]
    if channels == None
        channels = [1:nChannels]
    end
   
    for i=1:nChannels
        if contains(channels,i) == true
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
    segsToReject = (UTF8String => Array{Int,1})[]
    for i=1:length(eventList)
        segsToReject[eventList[i]] = []
        for j=1:size(rec[eventList[i]])[3]
            for k=1:length(channels) #size(rec[eventList[i]])[1]
                thisChan = channels[k]
                if (max(rec[eventList[i]][thisChan,:,j]) > thresh[k] || min(rec[eventList[i]][thisChan,:,j]) < -thresh[k]) == true
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

function getFRatios(ffts, compIdx, nSideComp, nExcludedComp, otherExclude)
    ##"""
    ##derived from get_F_Ratios2
    ##"""
    cnds = collect(keys(ffts))
    fftVals = (UTF8String => Any)[]
    fRatio = (UTF8String => Any)[]
    dfNum = 2
    dfDenom = 2*(nSideComp*2) -1
    for cnd in cnds
        fRatio[cnd] = (UTF8String => Array{Float64, 1})[]
        fftVals[cnd] = (UTF8String => Array{Float64, 1})[]
        fRatio[cnd]["F"] = []
        fRatio[cnd]["pval"] = []
        fftVals[cnd]["sigPow"] = []
        fftVals[cnd]["noisePow"] = []
        for c=1:length(compIdx)
            sideBands = getNoiseSidebands(compIdx, nSideComp, nExcludedComp, ffts[cnd]["mag"], otherExclude)
            noisePow = mean(sideBands[c])
            sigPow = ffts[cnd]["mag"][compIdx[c]]
            thisF =  sigPow/ noisePow
            fftVals[cnd]["sigPow"] = vcat(fftVals[cnd]["sigPow"], sigPow)
            fftVals[cnd]["noisePow"] = vcat(fftVals[cnd]["noisePow"], noisePow)
            fRatio[cnd]["F"] = vcat(fRatio[cnd]["F"], thisF)
            fRatio[cnd]["pval"] = vcat(fRatio[cnd]["pval"], pdf(FDist(dfNum, dfDenom), thisF))
        end
    end
    return fftVals, fRatio
end

function getNoiseSidebands(components, nCompSide, nExcludeSide, fftArray, otherExclude)
    #"""
    #the 2 has the possibility to exclude extra components, useful for distortion products
    #"""
    #components: a list containing the indexes of the target components
    #nCompSide: number of components used for each side band
    #n_exclude_side: number of components adjacent to to the target components to exclude
    #fft_array: array containing the fft values
    idxProtect = []; idxProtect = vcat(idxProtect, components)
    if otherExclude != None
        idxProtect = vcat(idxProtect, otherExclude)
    end
    for i=1:length(nExcludeSide)
        idxProtect = vcat(idxProtect, components + (i+1))
        idxProtect = vcat(idxProtect, components - (i+1))
    end
    #idxProtect = sorted(idxProtect)
    #print(idxProtect)

    noiseBands = []
    for i=1:length(components)
        loSide = []
        hiSide = []
        counter = 1
        while length(hiSide) < nCompSide
            currIdx = components[i] + nExcludeSide + counter
            if contains(idxProtect, currIdx) == false
                hiSide = vcat(hiSide, fftArray[currIdx])
            end
            counter = counter + 1
        end
        counter = 1
        while length(loSide) < nCompSide
            currIdx = components[i] - nExcludeSide - counter
            if contains(idxProtect, currIdx) == false
                loSide = vcat(loSide, fftArray[currIdx])
            end
            counter = counter + 1
        end
        noiseBands = vcat(noiseBands, loSide+hiSide)
    end
        
    return noiseBands
end
                              

function getSpectrum(sig, sampRate::Int, window::String, powerOfTwo::Bool)
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
         p[1:length(p)] = p[1:length(p)] * 2
    else
        p[1:(length(p)-1)] = p[1:length(p) - 1] * 2 # we"ve got even number of points fft
    end

    freq_array = [0:(nUniquePts-1)] * (sampRate / nfft)
    x = (ASCIIString => Array{Float64,1})[]
    x["freq"] = freq_array; x["mag"] = p
    return x
end


function mergeEventTableCodes(eventTable::Dict{String,Any}, trigList, newTrig::Int)
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
        rec[code] = deleteSlice3D(rec[code], toRemove[code])
    end
    return rec
end

function removeSpuriousTriggers(eventTable::Dict{String, Any}, sentTrigs::Array{Int}, minTrigDur::Real)
    recTrigs = eventTable["code"]
    recTrigsStart = eventTable["idx"]
    recTrigsDur = eventTable["dur"]

    allowedTrigs = unique(sentTrigs)
    allowedIdx = findin(recTrigs, allowedTrigs)
    recTrigsDur = recTrigsDur[allowedIdx]
    recTrigsStart = recTrigsStart[allowedIdx]
    recTrigs = recTrigs[allowedIdx]

    durCondition = recTrigsDur .>= minTrigDur
    recTrigs = recTrigs[durCondition]
    recTrigsStart = recTrigsStart[durCondition]
    recTrigsDur = recTrigsDur[durCondition]

    if recTrigs == sentTrigs
        match_found = true
    else
        match_found = false
    end


    eventTable["code"] = recTrigs
    eventTable["dur"] = recTrigsDur
    eventTable["idx"] = recTrigsStart

    resInfo = (UTF8String => Any)[]
    resInfo["match"] = match_found
    resInfo["lenSent"] = length(sentTrigs)
    resInfo["lenFound"] = length(recTrigs)

    return resInfo
end

function rerefCnt(rec, refChan::Int, channels)
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

    if contains(channels, refChan)
     rec[refChan,:] = rec[refChan,:] - rec[refChan,:]
    end
    
    return
end

function saveFRatios(fileName::String, subj::String, FRatio, fftValues, cndsTrigs, cndsLabels, nCleanByBlock, nRawByBlock)
    ## """
    
    ## Parameters
    ## ----------

    ## Returns
    ## ----------

    ## Examples
    ## ----------
    ## """
    ## #cnds = list(FRatio.keys())
    
    nRaw = (UTF8String => Int64)[]
    nClean = (UTF8String => Int64)[]
    for cnd in cndsTrigs
        nRaw[cnd] = 0
        nClean[cnd] = 0
        for blk=1:length(nCleanByBlock)
            nRaw[cnd] = nRaw[cnd] + nRawByBlock[blk][cnd]
            nClean[cnd] = nClean[cnd] + nCleanByBlock[blk][cnd]
        end
    end
               
    subjVec = (UTF8String)[]
    compVec = (Int64)[]
    conditionVec = (UTF8String)[]
    nRawVec = (Int64)[]
    nCleanVec = (Int64)[]
    FRatioVec = (Float64)[]
    sigPowVec = (Float64)[]
    noisePowVec = (Float64)[]
    pValVec = (Float64)[]
    #percRej = (Array{Float64,1})[]
            
    for i=1:length(cndsTrigs)
        thisN = length(FRatio[cndsTrigs[i]]["F"])
        subjVec = vcat(subjVec, [subj for j=1:thisN])
        conditionVec = vcat(conditionVec, [cndsLabels[i] for j=1:thisN])
        compVec = vcat(compVec, [1:thisN])
        sigPowVec = vcat(sigPowVec, fftValues[cndsTrigs[i]]["sigPow"][:])
        noisePowVec = vcat(noisePowVec, fftValues[cndsTrigs[i]]["noisePow"][:])
        pValVec = vcat(pValVec, FRatio[cndsTrigs[i]]["pval"][:])
        FRatioVec = vcat(FRatioVec, FRatio[cndsTrigs[i]]["F"][:])
        nRawVec = vcat(nRawVec, [nRaw[cndsTrigs[i]] for j=1:thisN])
        nCleanVec = vcat(nCleanVec, [nClean[cndsTrigs[i]] for j=1:thisN])
    end
                
    datsFrame = DataFrame({"subj" => subjVec,
            "condition" => conditionVec,
            "comp" => compVec,
            "fRatio"=> FRatioVec,
            "pval" => pValVec,
            "sigPow" => sigPowVec,
            "noisePow" => noisePowVec,
            "nRaw" => nRawVec,
            "nClean" => nCleanVec})
            
            datsFrame["percRej"] = 100-((float(datsFrame["nClean"]) ./ float(datsFrame["nRaw"])) * 100)
            #percRej = (datsFrame["nClean"] ./ datsFrame["nRaw"]) 
            writetable(fileName, datsFrame, separator=';')
end

function segment(rec, eventTable::Dict{String, Any}, epochStart::Real, epochEnd::Real, sampRate::Int, eventsList=None, eventsLabelsList=None)

    trigs = eventTable["code"]
    trigs_pos = eventTable["idx"]
    if eventsList == None
        eventsList = unique(trigs)
    end

    if eventsLabelsList == None
        eventsLabelsList = UTF8String[string(eventsList[i]) for i=1:length(eventsList)]
        end
        

    epochStartSample = int(round(epochStart*sampRate))
    epochEndSample = int(round(epochEnd*sampRate)) - 1

    nSamples = epochEndSample - epochStartSample + 1
    segs = (UTF8String => Array{eltype(rec),3})[]
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
    nSegs = (UTF8String => Int64)[]
    for i=1:length(eventsList) #count
        nSegs[eventsLabelsList[i]] = size(segs[eventsLabelsList[i]])[3]
    end

    return segs, nSegs
        
end


function combineChained(dList)

    cnds = collect(keys(dList[1])) 
    cmb = (String => Array{Float32, 2})[]
    for cnd in cnds
        for i=1:length(dList)
            if i == 1
                cmb[cnd] = deepcopy(dList[1][cnd])
            else
                cmb[cnd] = cmb[cnd] + dList[i][cnd]
            end
        end
        cmb[cnd] = cmb[cnd] ./ length(dList)
    end
            
    return cmb
end

#end
