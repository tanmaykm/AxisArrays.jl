# shift, lead and lag
# shift given axis by given intervals
# throw away first n if n > 0
# throw aray last n if n < 0

function shiftaxis{T,N,D,Ax,name,B}(data::AxisArray{T,N,D,Ax}, ax::Axis{name,B})
    by = ax.val
    d = axisdim(data, Axis{name})
    oldvals = axes(data, Axis{name}).val
    newvals = oldvals .+ by
    r = searchsortedfirst(newvals, oldvals[1]):searchsortedlast(newvals, oldvals[end])

    allaxes = [axes(data)...]
    allaxes[d] = Axis{name}(newvals[r])

    allr = repmat(Any[:], N)
    allr[d] = r

    AxisArray(data.data[allr...], allaxes...) #::AxisArray{T,N,D,Ax}
end

# lag: take yesterday's value for today
lag{name}(data, by::Axis{name}) = shiftaxis(data, by)

# lead: take tomorrow's value for today
lead{name}(data, by::Axis{name}) = shiftaxis(data, Axis{name}(-(by.val)))
