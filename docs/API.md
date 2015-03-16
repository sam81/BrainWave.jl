# ElectroJulia

## Exported
---

#### baselineCorrect!{T<:Real}(rec::Dict{String, Array{T<:Real, 3}}, baselineStart::Real, preDur::Real, sampRate::Integer)
Perform baseline correction by subtracting the average pre-event
voltage from each channel of a segmented recording.

#### Parameters
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
baselineCorrect(rec=rec, baselineStart=-0.2, preDur=0.2, sampRate=512)
#now with a baseline shorter than pre_dur
baselineCorrect(rec, -0.15, 0.2, 512)
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:106](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

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
[ElectroJulia/src/ElectroJulia.jl:989](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)

## Internal
---

#### baselineCorrectloop!{T<:Real}(rec::Dict{String, Array{T<:Real, 3}}, baselineStart::Real, preDur::Real, sampRate::Integer)
Perform baseline correction by subtracting the average pre-event
voltage from each channel of a segmented recording.

#### Parameters
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
baselineCorrect(rec=rec, baselineStart=-0.2, preDur=0.2, sampRate=512)
#now with a baseline shorter than pre_dur
baselineCorrect(rec, -0.15, 0.2, 512)
```


**source:**
[ElectroJulia/src/ElectroJulia.jl:145](file:///home/sam/.julia/v0.3/ElectroJulia/src/ElectroJulia.jl)


