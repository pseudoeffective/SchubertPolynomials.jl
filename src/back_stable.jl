# Tools for computing back stable Schubert polynomials in Julia, using BPDs
# David Anderson, June 2024


export back_schub_poly, acoeff, schub2dom


#############

"""
    back_schub_poly(w; double=false, ring=schub_ring(length(w)-1, length(w)-1))

Return the back stable Schubert polynomial for the permutation `w`

## Arguments
- `w::Vector{Int}`: A permutation

## Keywords
- `double::Bool`: when `true`, coefficients are the double/equivariant polynomials in the y-variables of `ring`.  Default `false`.
- `ring::MPolyRing`: The ambient polynomial ring (build with `schub_ring`).  Defaults to `schub_ring(length(w)-1, length(w)-1)`.

## Returns
`DominantSum`: the back stable Schubert polynomial, written in the basis of dominant polynomials, with coefficients in `ring`

## Example

```julia-repl
julia> R = schub_ring(4,4);

julia> back_schub_poly( [1,3,4,2]; ring=R, double=true )
(x2*x3 + x2*y2 + x3*y2 + y2^2)*ss[] + ss[1, 1] + (x3 + y2)*ss[1]
```
Here `ss[lambda]` stands for the dominant polynomial corresponding to a partition `lambda`.
"""
function back_schub_poly(w; double::Bool=false,
                         ring::MPolyRing=schub_ring( max(length(w)-1,1), max(length(w)-1,1) ) )

  w=trimw(w)

  pars = Vector{Int}[]
  cfs = Union{Int,ZZMPolyRingElem}[]

  bpds = all_bpds(w)

  for b in bpds
    la = dominant_part(b)
    if !(la in pars)
      aa = acoeff(w,la; double=double, ring=ring)
      push!(pars, la)
      push!(cfs, aa)
    end
  end

  return sort_dominant_sum( DominantSum(cfs,pars) )
end



# get the coefficient of the (dominant) s[lambda] in the (back stable) Schubert polynomial S[w]
function acoeff( w::Vector{Int}, lambda::Vector{Int}; double::Bool=false,
                 ring::MPolyRing=schub_ring(length(w)-1, length(w)-1) )

  lam = trimp(lambda)

  apol = ring(0)


  for b in all_bpds(w)
    if dominant_part(b)==lam
      apol = apol + bpd2bin_star( b, ring; double=double )
    end
  end

  return apol
end


function bpd2bin_star( bpd::BPD, R::MPolyRing=schub_ring( size(bpd.mtx)[1]-1, size(bpd.mtx)[2]-1 ); double::Bool=false )
# product of binomials for bpd
# only schub version now
# `double=false` (default) ignores any y-variables, giving the single polynomial
  local n=size(bpd.mtx)[1]-1
  bin = R(1)

  x = extract_vars(R; varname=:x)
  y = extract_vars(R; varname=:y)

  local aa=length(x)
  local bb= double ? length(y) : 0

  la = dominant_part(bpd)
  la = vcat( la, zeros(Int, n-length(la)) )

  for i=1:n
    for j=la[i]+2:n

      if bpd.mtx[i,j]==0
        p=R(0)
        if i<=aa
          p=p+x[i]
        end
        if j<=bb
          p=p+y[j]
        end

        bin = bin*p

      end

    end
  end

  return bin

end

"""
    schub2dom(SS; ring=schub_ring(maxlen-1, maxlen-1), double=false)

Convert a `SchubertSum` to a `DominantSum` by expanding each Schubert polynomial S[w]
as a back stable Schubert polynomial in the dominant basis.

## Arguments
- `SS::SchubertSum`: A sum of Schubert polynomials

## Keywords
- `ring::MPolyRing`: The ambient polynomial ring (build with `schub_ring`).  Defaults to a ring sized by the largest permutation in `SS`.
- `double::Bool`: when `true`, coefficients are equivariant (in the y-variables of `ring`).  Default `false`.

## Returns
`DominantSum`: the sum expressed in the basis of dominant polynomials

## Example

```julia-repl
julia> SS = SchubertSum([1, 1], [[2,1,3], [1,3,2]]);

julia> schub2dom(SS; ring=schub_ring(4,4), double=true)
```
"""
function schub2dom(SS::SchubertSum;
                   ring::MPolyRing=schub_ring( max(maximum(length.(SS.schubs))-1, 1),
                                               max(maximum(length.(SS.schubs))-1, 1) ),
                   double::Bool=false)

    # Start with the zero DominantSum in the polynomial ring
    result = DominantSum([ring(0)], [[0]])

    # For each term c*S[w] in the SchubertSum
    for (c, w) in zip(SS.coeffs, SS.schubs)
        # Convert S[w] to dominant basis
        dom_poly = back_schub_poly(w; ring=ring, double=double)

        # Multiply by the coefficient and add to result
        result = result + c * dom_poly
    end

    return condense(result)
end
