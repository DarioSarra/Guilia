using Survival: KaplanMeier
"""
`_survival(times,events)`
times is a vector of the time of the event
events is a vector of booleans wheter the event occured or was censored
"""
function _survival(times, events; axis = discrete_axis(times), kwargs...)
    if eltype(events) <: AbstractString
        try
            events = parse.(Bool,events)
        catch e
            println(e)
        end
    end
    km = Distributions.fit(KaplanMeier, times, events)
    surv = zeros(length(axis))
    for (i, ax) in enumerate(axis)
        if ax < km.times[1]
            surv[i] = 1
        elseif ax > km.times[end]
            surv[i] = km.survival[end]
        else
            surv[i] = km.survival[searchsortedfirst(km.times, ax)]
        end
    end
    return zip(axis, surv)
end

const Survival = Recombinase.Analysis((; discrete = _survival))
