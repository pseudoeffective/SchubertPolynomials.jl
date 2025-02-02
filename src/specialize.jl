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


function principal_specialization( pol::ZZMPolyRingElem, R::Union{ZZMPolyRing,DoublePolyRing}=parent(pol) )

   if isa(R,DoublePolyRing)
     RR = R.ring
     xx = R.x_vars
   else
     RR = R
     xx = gens(R)
   end

   S,qq = polynomial_ring(RR,[:q])

   n = maxvar(pol)

   if n==0 return(pol) end

   return evaluate( pol, [xx[i] for i=1:n], [ qq[1]^(i-1) for i=1:n ] )

end


function principal_specialization( w::Vector{Int}, R::DoublePolyRing=xy_ring(length(w)-1)[1] )
# return principal specialization of Schubert polynomial for w

  if length(R.y_vars)>0
    throw(ArgumentError("no y variables allowed"))
  end

  spw = schub_poly( w, R )

  spq = principal_specialization( spw, R )

  return spq

end
