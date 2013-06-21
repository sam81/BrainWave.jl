#module electrojulia

#export segment
using PyCall
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
    weightedAve = (UTF8String => Array{Float32,2})[]
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
        weightedAve[event] = zeros(size(aveList[1][event]))
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
    ave = (UTF8String => Array{Float32,2})[]
    nSegs = (UTF8String => Int)[]
    for i=1:length(eventList)
        nSegs[eventList[i]] = size(rec[eventList[i]])[3]
        ave[eventList[i]] = mean(rec[eventList[i]], 3)[:,:,1]
    end
    return ave, nSegs
end

function baselineCorrect(rec, baselineStart, preDur, sampRate)
    
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
    epochStartSample = int(round(preDur*sampRate))+1
    baselineStartSample = int(epochStartSample - abs(round(baselineStart*sampRate)))
    #print(baselineStartSample, " ", epochStartSample, "\n")
    
    for i=1:length(eventList) #for each event
        for j=1:size(rec[eventList[i]])[3] #for each epoch
            for k=1: size(rec[eventList[i]])[1] #for each electrode
                thisBaseline = mean(rec[eventList[i]][k,baselineStartSample:epochStartSample,j])
                rec[eventList[i]][k,:,j] = rec[eventList[i]][k,:,j] - thisBaseline
            end
        end
    end
    
end

function chainSegments(rec, nChunks, sampRate, startTime, endTime, baselineDur)
    """
    Take a dictionary containing in each key a list of segments, and chain these segments
    into chunks of length nChunks
    baselineDur is for determining what is the zero point
    startTime and endTime are given with reference to the zero point
    """
    baselinePnts = round(baselineDur * sampRate)
    startPnt = int(round(startTime*sampRate) + baselinePnts) 
    endPnt = int(round(endTime*sampRate) + baselinePnts) 
    chunkSize = ((endPnt - startPnt)+1)
    sweepSize = chunkSize * nChunks
    nReps = (UTF8String => Array{Int,1})[]
    eventList = collect(keys(rec))
    eegChained = (UTF8String => Array{Float32,2})[]
    fromeegChainedAve = (UTF8String => Array{Float32,2})[]
    for i=1:length(eventList)
        currCode = eventList[i]
        eegChained[currCode] = zeros(size(rec[currCode])[1], sweepSize)  #two-dimensional array of zeros
        #fromeegChainedAve[currCode] = zeros(size(rec[currCode])[1], chunkSize)
        nReps[currCode] = zeros(nChunks)
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
                if (max(rec[eventList[i]][thisChan,:,j]) > thresh[k] || min(rec[eventList[i]][k,:,j]) < -thresh[k])
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

function getSpectrum(sig, sampRate, window, powerOfTwo)
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
    if nfft % 2 > 0 # we've got odd number of points fft
         p[1:length(p)] = p[1:length(p)] * 2
    else
        p[1:(length(p)-1)] = p[1:length(p) - 1] * 2 # we've got even number of points fft
    end

    freq_array = [0:(nUniquePts-1)] * (sampRate / nfft)
    x = (ASCIIString => Array{Float64,1})[]
    x["freq"] = freq_array; x["mag"] = p
    return x
end


function mergeEventTableCodes(eventTable, trigList, newTrig)
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

function nextPowTwo(x)
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

function removeSpuriousTriggers(eventTable, sentTrigs, minTrigDur)
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

function rerefCnt(rec, refChan, channels=None)
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

function segment(rec, eventTable, epochStart, epochEnd, sampRate, eventsList=None, eventsLabelsList=None)

    trigs = eventTable["code"]
    trigs_pos = eventTable["idx"]
    if eventsList == None
        eventsList = unique(trigs)
    end

    if eventsLabelsList == None
        eventsLabelsList = UTF8String[string(eventsList[i]) for i=1:length(eventsList)]
        end
        

    epochStartSample = int(round(epochStart*sampRate))
    epochEndSample = int(round(epochEnd*sampRate))

    nSamples = epochEndSample - epochStartSample + 1
    segs = (UTF8String => Array{Float32,3})[]
    for i=1:length(eventsList)
        idx = trigs_pos[trigs .== eventsList[i]]
        
        segs[eventsLabelsList[i]] = zeros(Float32, size(rec)[1], nSamples, length(idx))
        for j=1:length(idx)
            thisStartPnt = (idx[j]+epochStartSample)
            #print(thisStartPnt)
            thisStopPnt = (idx[j]+epochEndSample)
            if thisStartPnt < 0 || thisStopPnt > size(rec)[2]
                if thisStartPnt < 0
                    print(idx[j], " Epoch starts before start of recording. Skipping")
                end
                if thisStopPnt > size(rec)[2]#rec.shape[1]
                    print(idx[j], " Epoch ends after end of recording. Skipping")
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


#end
