# Polynomial ring helpers for SchubertPolynomials.jl
#
# The package now works with a bare `MPolyRing` (mirroring `ssyt_ring` in
# SemistandardTableaux): variables are named x1..xn then y1..ym, and the x- and
# y-families are recovered by name with `extract_vars` (re-exported from
# SemistandardTableaux).  The old `DoublePolyRing` type and `xy_ring` constructor
# are kept for one version as deprecated shims.


export schub_ring


#################
# Polynomial ring constructor
#################

"""
    schub_ring(n, m=0; coeff=ZZ, xname=:x, yname=:y) -> MPolyRing

Polynomial ring with variables `x1..xn` then `y1..ym` over the coefficient ring
`coeff` (default Nemo `ZZ`).  Mirrors `ssyt_ring` from SemistandardTableaux.

Returns the ring itself (an `MPolyRing`), not a tuple of variable families; use
`extract_vars(R; varname=:x)` / `extract_vars(R; varname=:y)` to recover them.

A single (ordinary) polynomial lives in any such ring — pass `double=false` to a
constructor and only the x-variables are used.  For a double/factorial
polynomial, build a ring with y-variables (`m > 0`) and pass `double=true`.

# Examples
```julia-repl
julia> R = schub_ring(3, 2)
Multivariate polynomial ring in 5 variables x1, x2, x3, y1, y2 over integers

julia> extract_vars(R; varname=:x)
3-element Vector{ZZMPolyRingElem}:
 x1
 x2
 x3

# Over the rationals
julia> R = schub_ring(3; coeff=QQ);
```
"""
function schub_ring(n::Int, m::Int=0; coeff=ZZ, xname::Symbol=:x, yname::Symbol=:y)
    names = vcat(["$(xname)$(i)" for i in 1:n], ["$(yname)$(j)" for j in 1:m])
    R, _ = polynomial_ring(coeff, names)
    return R
end


################################################################################
# Deprecated: DoublePolyRing + xy_ring
#
# Retained for one version so existing scripts keep working.  Migrate
#   R = xy_ring(n, m)[1]      ->  R = schub_ring(n, m)
#   R.x_vars                  ->  extract_vars(R; varname=:x)
#   R.y_vars                  ->  extract_vars(R; varname=:y)
# and pass the ring as the `ring` keyword instead of positionally.
################################################################################

export DoublePolyRing, xy_ring

struct DoublePolyRing
           ring::ZZMPolyRing
           x_vars::Vector{ZZMPolyRingElem}
           y_vars::Vector{ZZMPolyRingElem}
       end


function xy_ring(xx::Vector{String}, yy::Vector{String})
    Base.depwarn("`xy_ring`/`DoublePolyRing` are deprecated; use `schub_ring` " *
                 "(returns an MPolyRing) with `extract_vars`.", :xy_ring)
    n=length(xx)
    m=length(yy)
    R, all_vars = polynomial_ring(ZZ, vcat(xx,yy))
    x = all_vars[1:n]
    y = all_vars[n+1:n+m]
    return DoublePolyRing(R,x,y), x, y
end

xy_ring(xx::Vector{String}) = xy_ring(xx, String[])

function xy_ring(n::Int, m::Int)
    xvars = ["x$(i)" for i=1:n]
    yvars = ["y$(i)" for i=1:m]
    xy_ring(xvars, yvars)
end

xy_ring(n::Int) = xy_ring(n, 0)
