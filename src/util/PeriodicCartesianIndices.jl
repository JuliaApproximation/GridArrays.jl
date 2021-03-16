module PeriodicIndexing

export PeriodicCartesianIndices

using Base: tail

import Base: getindex, length, size, IndexStyle, IndexCartesian, in, iterate, mod,
    first, last

mod_tuple(c::NTuple{N,Int}, size::NTuple{N,Int}) where N  =
    ntuple(k->mod(c[k]-1,size[k])+1, Val(N))

mod(c::CartesianIndex{N}, size::NTuple{N,Int}) where N =
    CartesianIndex(mod_tuple(c.I, size))

add_offset_mod(i::NTuple{N,Int}, start::NTuple{N,Int}, size::NTuple{N,Int}, periodic::NTuple{N,Bool}) where N =
    ntuple(k->( periodic[k]  ? add_offset_mod(i[k], start[k], size[k]) : i[k] ), Val(N))

add_offset_mod(i::Int, start::Int, size::Int) =
    start + mod(i-start, size)

struct PeriodicCartesianIndices{N,K} <:AbstractArray{CartesianIndex{N},N}
    size::NTuple{N,Int}
    iter::CartesianIndices{N,K}
    start::CartesianIndex{N}
    stop::CartesianIndex{N}
    periodic::NTuple{N,Bool}
end

getindex(A::PeriodicCartesianIndices{N}, I::Vararg{Int,N}) where N = mod(getindex(A.iter,I...),A.size)
IndexStyle(A::PeriodicCartesianIndices) = IndexCartesian()
function PeriodicCartesianIndices(size::NTuple{N,Int}, iranges::UnitRange...) where {N}
    start = CartesianIndex(ntuple(k->first(iranges[k]), Val(N)))
    stop  = CartesianIndex(ntuple(k->last(iranges[k]), Val(N)))
    PeriodicCartesianIndices(size, CartesianIndices(iranges), start, stop, ntuple(k->true,Val(N)))
end

function PeriodicCartesianIndices(size::NTuple{N,Int}, start::CartesianIndex{N}, stop::CartesianIndex,  periodic::NTuple{N,Bool}=ntuple(k->true,Val(N))) where {N}
    iranges = ntuple(k-> ((start[k] <= stop[k]) ? (start[k]:stop[k]) : (start[k]:(stop[k]+size[k]))), Val(N))
    PeriodicCartesianIndices(size, CartesianIndices(iranges), start, stop, periodic)
end

size(m::PeriodicCartesianIndices) = size(m.iter)
length(m::PeriodicCartesianIndices) = length(m.iter)


struct PeriodicUnitRange{I} <:AbstractUnitRange{I}
    size::I
    iter::UnitRange{I}
    periodic::Bool

    function PeriodicUnitRange(size::I, iter::UnitRange{I}, periodic=true) where I<:Integer
        if periodic
            e = last(iter)
            l = length(iter)
            e_aligned = mod(e-1,size)+1
            f_aligned = e_aligned - l + 1
            new{I}(size, f_aligned:e_aligned, periodic)
        else
            new{I}(size, max(1,start(iter)):min(size,last(iter)), periodic)
        end
    end
end

first(r::PeriodicUnitRange) = mod(first(r.iter)-1,r.size)+1
last(r::PeriodicUnitRange) = mod(last(r.iter)-1,r.size)+1

function in(i::Integer, r::PeriodicUnitRange{<:Integer})
    if r.periodic
        i_aligned =  mod(i-1,r.size)+1
        return (i_aligned ∈ r.iter) | (i_aligned ∈ r.iter .+ r.size)
    else
        return in(i, r.iter)
    end
end

function iterate(m::PeriodicUnitRange)
    tuple = iterate(m.iter)
    if tuple != nothing
        index, state = tuple
        mod(index-1, m.size)+1, state
    end
end

function iterate(m::PeriodicUnitRange, state)
    tuple = iterate(m.iter, state)
    if tuple != nothing
        index, state = tuple
        mod(index-1, m.size)+1, state
    end
end

length(m::PeriodicUnitRange) = length(m.iter)

nbindexlist(index::Base.CartesianIndex{N}, size::NTuple{N,Int}, periodic=ntuple(k->false, Val(N))) where {N}  =
    PeriodicCartesianIndices(size, index+CartesianIndex(ntuple(k->-1, Val(N))), index+CartesianIndex(ntuple(k->1, Val(N))), periodic)
nbindexlist(index::Int, size::NTuple{1,Int}, periodic=false)  =
    PeriodicUnitRange(size[1], index-1:index+1, periodic)


end # module
