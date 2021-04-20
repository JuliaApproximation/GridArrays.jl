
################
# Bounding boxes
#################

# A bounding box is an Interval or ProductDomain of intervals that encompasses the domain.

# If the boundingbox is not a product of intervals, something has gone wrong.

import DomainSets: boundingbox

@deprecate boundingbox(a::Number, b::Number) boundingbox(a..b)
@deprecate boundingbox(a::SVector{1}, b::SVector{1}) boundingbox(a[1]..b[1])
@deprecate boundingbox(a, b) boundingbox(Rectangle(a, b))
