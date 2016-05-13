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

# window, percentchange, moving
abstract AbstractWindow{N,D}
immutable TumblingWindow{N,D} <: AbstractWindow{N,D}
    wspec::Axis{N,D}
end
immutable SlidingWindow{N,D} <: AbstractWindow{N,D}
    wspec::Axis{N,D}
end

wd{N,D}(w::AbstractWindow{N,D}) = w.wspec.val
st{N,D}(w::TumblingWindow{N,D}) = w.wspec.val
st{N,D}(w::SlidingWindow{N,D}) = D(1)
wp{N,D}(w::AbstractWindow{N,D}, v) = v .. (v + wd(w) - D(1))

ax{N,D}(data::AxisArray, w::AbstractWindow{N,D}) = ax(axes(data, Axis{N}), w)
function ax{N,T,D}(ax::Axis{N,T}, w::AbstractWindow{N,D})
    r = ax.val[1]:st(w):ax.val[end]
    r[1]:st(w):(r[end] - wd(w))
end

Base.length{N,D}(data::AxisArray, w::AbstractWindow{N,D}) = count(axes(data, Axis{N}), w)
function Base.length{N,T,D}(ax::Axis{N,T}, w::AbstractWindow{N,D})
    ax1 = ax.val[1]
    ax2 = ax.val[end] - wd(w)
    s = st(w)
    n = 0
    while (ax1+s) <= ax2
        n += 1
    end
    n
end

@generated function window{T,N,D,Ax,name,B}(data::AxisArray{T,N,D,Ax}, wspec::AbstractWindow{name,B}, fn)
    axdim = axisdim(AxisArray{T,N,D,Ax}, Axis{name})
    iterdims = deleteat!([1:N;], axdim)
    axnames = axisnames(AxisArray{T,N,D,Ax})
    # ref: https://github.com/JuliaLang/julia/issues/13359
    axargs = [:(axvals[$d]) for d = 1:N]

    X = quote
        wwidth = wd(wspec)
        wstep = st(wspec)
        Aw = Axis{name}

        ax1, ax2 = extrema(axes(data, Aw).val)
        newax = Aw(ax1:wstep:ax2)
        res = similar(data, newax)
        axvals = Array(Any, $N)
    end

    v = gensym("i_$axdim")
    L = quote
        axvals[$axdim] = Colon()
        dtmp = data[$(axargs...)]
        dvec = AxisArray(reshape(dtmp.data, length(dtmp.data)), axes(data)[$axdim])
        Atype = Axis{$axnames[$axdim]}
        for $v in 1:length(newax)
            wdata = dvec[Atype(wp(wspec, newax.val[$v]))].data
            axvals[$axdim] = $v
            res[$(axargs...)] = fn(wdata)
        end
    end

    for idx in reverse(iterdims)
        v = gensym("i_$idx")
        L = :(for $v in 1:length(axes(data, Axis{$axnames[$idx]}))
                axvals[$idx] = $v
                $L
            end)
    end
    push!(X.args, L)
    push!(X.args, :(return res))
   
    X 
end

moving(data, ax::Axis, fn) = window(data, SlidingWindow(ax), fn)
