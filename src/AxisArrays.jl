module AxisArrays

using Requires, Tuples, RangeArrays

export AxisArray, Axis, Interval, axisnames, axisvalues, axisdim, axes, .., atindex,
        shiftaxis, lead, lag, window, TumblingWindow, SlidingWindow, moving, collapse, resample, percentchange, mapaxis, per, cat, permutedims, merge

include("core.jl")
include("intervals.jl")
include("search.jl")
include("indexing.jl")
include("sortedvector.jl")
include("utils.jl")
include("timeseries.jl")

end
