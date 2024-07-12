# Functions for specializing Schubert polynomials
# David Anderson, July 2024.


export localization, principal_specialization

# TO DO: docs, tests

#########


function localization( u::Vector{Int}, v::Vector{Int}, R::DoublePolyRing=xy_ring( length(v),length(v) )[1] )
# localize Schubert class u at fixed point v
  x=R.x_vars
  y=R.y_vars

  uu=trimw(u) 
  vv=trimw(v)

  if length(x)<length(vv) || length(y)<length(vv)
    throw(ArgumentError("not enough variables"))
  end

  if length(x)!=length(y)
    throw(ArgumentError("incompatible x and y variables"))
  end

  yv = [ -y[vi] for vi in vv ]

  sp = schub_poly(uu,R)

  spv = evaluate( sp, x[1:length(vv)], yv )

  return spv

end


function principal_specialization( w::Vector{Int}, R::DoublePolyRing=xy_ring(length(w)-1)[1] )
# return principal specialization of Schubert polynomial for w

  x=R.x_vars
  y=R.y_vars

  S,q = polynomial_ring( R.ring, "q" )

  spw = schub_poly( w, R )

  spq = evaluate( spw, vcat(x,y), vcat( [q^(i-1) for i in 1:length(x)], [S(0) for i in 1:length(y)] ) )

  return spq

end
