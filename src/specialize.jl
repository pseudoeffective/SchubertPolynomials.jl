# Functions for specializing Schubert polynomials
# David Anderson, July 2024.


export localization, principal_specialization

# TO DO: docs, tests

#########


function localization( u::Vector{Int}, v::Vector{Int};
                       double::Bool=true, ring::MPolyRing=schub_ring( length(v), length(v) ) )
# localize Schubert class u at fixed point v
  x=extract_vars(ring; varname=:x)
  y=extract_vars(ring; varname=:y)

  uu=trimw(u)
  vv=trimw(v)

  if length(x)<length(vv) || length(y)<length(vv)
    throw(ArgumentError("not enough variables"))
  end

  if length(x)!=length(y)
    throw(ArgumentError("incompatible x and y variables"))
  end

  yv = [ -y[vi] for vi in vv ]

  sp = schub_poly(uu; ring=ring, double=double)

  spv = evaluate( sp, x[1:length(vv)], yv )

  return spv

end


function principal_specialization( pol::ZZMPolyRingElem; ring::MPolyRing=parent(pol) )

   RR = ring
   xx = extract_vars(ring; varname=:x)

   S,qq = polynomial_ring(RR,[:q])

   n = maxvar(pol)

   if n==0 return(pol) end

   return evaluate( pol, [xx[i] for i=1:n], [ qq[1]^(i-1) for i=1:n ] )

end


function principal_specialization( w::Vector{Int}; double::Bool=false,
                                   ring::MPolyRing=schub_ring(length(w)-1, 0) )
# return principal specialization of Schubert polynomial for w

  if length(extract_vars(ring; varname=:y))>0
    throw(ArgumentError("no y variables allowed"))
  end

  spw = schub_poly( w; ring=ring, double=double )

  spq = principal_specialization( spw; ring=ring )

  return spq

end
