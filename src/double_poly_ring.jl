
#==============================================================================================================================+
	
	double_poly_ring.jl
	David Anderson, June 2024.

+==============================================================================================================================#

export DoublePolyRing, xy_ring

#==============================================================================================================================#

# Polynomial constructors

######
struct DoublePolyRing
           ring::ZZMPolyRing
           x_vars::Vector{ZZMPolyRingElem}
           y_vars::Vector{ZZMPolyRingElem}
       end

#####



"""
    xy_ring(xx,yy)

Form a DoublePolyRing object

## Arguments
- `xx::Vector{String}`: the names of x variables
- `yy::Vector{String}`: the names of y variables

If xx, yy are given as Int, the names x1, x2, y1, y2, etc., are chosen.

## Returns
- `DoublePolyRing`, `Vector{ZZMPolyRingElem}`, `Vector{ZZMPolyRingElem}`: The double polynomial ring, with specified variable sets.

# Examples
```julia-repl
# Define a DoublePolyRing with x-variables [a,b,c] and y-variables [\u03B1,\u03B2,\u03B3]
julia> xx = ["a","b","c"]; yy = ["\u03B1","\u03B2","\u03B3"];

julia> R,x,y = xy_ring(xx,yy)
(DoublePolyRing(Multivariate polynomial ring in 6 variables over ZZ, ZZMPolyRingElem[a, b, c], ZZMPolyRingElem[\u03B1,\u03B2,\u03B3]), ZZMPolyRingElem[a, b, c], ZZMPolyRingElem[\u03B1,\u03B2,\u03B3])

# Define a DoublePolyRing with three x-variables and two y-variables
julia> R = xy_ring(3,2)[1]
DoublePolyRing(Multivariate polynomial ring in 5 variables over ZZ, ZZMPolyRingElem[x1, x2, x3], ZZMPolyRingElem[y1, y2])

# Extract the x variables
julia> x = R.x_vars
3-element Vector{ZZMPolyRingElem}:
 x1
 x2
 x3

# Extract the y variables
julia> y = R.y_vars
2-element Vector{ZZMPolyRingElem}:
 y1
 y2

# If only one argument is supplied, the result has empty y-variable set
julia> R,x,y = xy_ring(3)
(DoublePolyRing(Multivariate polynomial ring in 3 variables over ZZ, ZZMPolyRingElem[x1, x2, x3], ZZMPolyRingElem[]), ZZMPolyRingElem[x1, x2, x3], ZZMPolyRingElem[])

```

# Notes

The ring is constructed from a `ZZMPolyRing` object in Nemo, so elements are of type `ZZMPolyRingElem`.

"""
function xy_ring(xx::Vector{String},yy::Vector{String})

    n=length(xx)
    m=length(yy)

    R,all_vars = polynomial_ring(ZZ,vcat(xx,yy))

    x = all_vars[1:n]
    y = all_vars[n+1:n+m]

    return DoublePolyRing(R,x,y), x, y

end


function xy_ring(xx::Vector{String})
    return xy_ring(xx,String[])
end


function xy_ring(n::Int,m::Int)

  local xvars = ["x$(i)" for i=1:n]
  local yvars = ["y$(i)" for i=1:m]

  R,all_vars = polynomial_ring(ZZ,vcat(xvars,yvars))

  x = all_vars[1:n]
  y = all_vars[n+1:n+m]

  return DoublePolyRing(R,x,y), x, y

end


function xy_ring(n::Int)
  return xy_ring(n,0)
end


