
"""
Make a grid with the given size, for which the given domain is the covering domain.
"""
function discretize_togrid end

discretize_togrid(d::AbstractInterval, n::Int) = EquispacedGrid(n, d)
discretize_togrid(d::Interval{:closed,:open}, n::Int) = PeriodicEquispacedGrid(n, d)
discretize_togrid(d::Interval{:open,:open}, n::Int) = MidpointEquispacedGrid(n, d)
discretize_togrid(d::UnitInterval{T}, n::Int) where {T} = UnitEquispacedGrid{T}(n)
discretize_togrid(d::ChebyshevInterval{T}, n::Int) where {T} = ChebyshevNodes{T}(n)

discretize_togrid(d::ProductDomain, dims) = productgrid(map(discretize_togrid, components(d), dims)...)

discretize_togrid(d::MappedDomain, dims) =
    MappedGrid(discretize_togrid(superdomain(d), dims), forward_map(d))

discretize_togrid(d::UnitCircle{T}, n::Int) where {T} =
    MappedGrid(UnitPeriodicEquispacedGrid{T}(n), fromcanonical(d, DomainSets.Parameterization()))
