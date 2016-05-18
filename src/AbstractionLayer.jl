using Compat
import Compat.String

type RawEEG{S<:Real, P<:String}
    data::AbstractMatrix{S}
    eventTable::Any
    sampRate::Real
    nChannels::Int
    nSamples::Int
    chanLabels::AbstractVector{P}
end

function toRawEEG{S<:Real}(data::AbstractMatrix{S}, eventTable, sampRate::Real, ; chanLabels=["" for i=1:size(data)[1]])
    thisRawEEG = RawEEG(data, eventTable, sampRate, size(data)[1], size(data)[2], chanLabels)

    return thisRawEEG
end

