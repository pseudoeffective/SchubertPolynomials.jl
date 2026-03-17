# Tools for computing back stable Schubert polynomials in Julia, using BPDs
# David Anderson, June 2024


export back_schub_poly, acoeff, schub2dom


#############

"""
    back_schub_poly(w, R)

Return the back stable Schubert polynomial for the permutation `w`

## Arguments
- `w::Vector{Int}`: A permutation
- `R::DoublePolyRing`: The ambient polynomial ring, with underlying ring `R.ring` and variables `R.x_vars` and `R.y_vars`


## Returns
`DominantSum`: the back stable Schubert polynomial, written in the basis of dominant polynomials, with coefficients in R

## Example

```julia-repl
julia> R = xy_ring(4,4)[1];

julia> back_schub_poly( [1,3,4,2], R )
(x2*x3 + x2*y2 + x3*y2 + y2^2)*ss[] + ss[1, 1] + (x3 + y2)*ss[1]
```
Here `ss[lambda]` stands for the dominant polynomial corresponding to a partition `lambda`.
"""
function back_schub_poly(w, R::DoublePolyRing=xy_ring( max(length(w)-1,1))[1] )

  w=trimw(w)

  pars = Vector{Int}[]
  cfs = Union{Int,ZZMPolyRingElem}[]

  bpds = all_bpds(w)

  for b in bpds
    la = dominant_part(b)
    if !(la in pars)
      aa = acoeff(w,la,R)
      push!(pars, la)
      push!(cfs, aa)
    end
  end

  return sort_dominant_sum( DominantSum(cfs,pars) )
end



# get the coefficient of the (dominant) s[lambda] in the (back stable) Schubert polynomial S[w]
function acoeff( w::Vector{Int}, lambda::Vector{Int}, R::DoublePolyRing=xy_ring(length(w)-1)[1] )

  lam = trimp(lambda)

  apol = R.ring(0)


  for b in all_bpds(w)
    if dominant_part(b)==lam
      apol = apol + bpd2bin_star( b, R )
    end
  end     

  return apol
end


function bpd2bin_star( bpd::BPD, R::DoublePolyRing=xy_ring( size(bpd.mtx)[1]-1, size(bpd.mtx)[2]-1 )[1] )
# product of binomials for bpd
# only schub version now
# can get single polyn by using no y_vars
  local n=size(bpd.mtx)[1]-1
  bin = R.ring(1)

  x = R.x_vars
  y = R.y_vars

  local aa=length(x)
  local bb=length(y)

  la = dominant_part(bpd)
  la = vcat( la, zeros(Int, n-length(la)) )

  for i=1:n
    for j=la[i]+2:n

      if bpd.mtx[i,j]==0
        p=R.ring(0)
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
    schub2dom(SS, R)

Convert a `SchubertSum` to a `DominantSum` by expanding each Schubert polynomial S[w]
as a back stable Schubert polynomial in the dominant basis.

## Arguments
- `SS::SchubertSum`: A sum of Schubert polynomials
- `R::DoublePolyRing`: The ambient polynomial ring (optional, will be inferred from the largest permutation if not provided)

## Returns
`DominantSum`: the sum expressed in the basis of dominant polynomials

## Example

```julia-repl
julia> R = xy_ring(4,4)[1];

julia> SS = SchubertSum([1, 1], [[2,1,3], [1,3,2]]);

julia> schub2dom(SS, R)
```
"""
function schub2dom(SS::SchubertSum, R::DoublePolyRing)

    # Start with the zero DominantSum in the polynomial ring
    result = DominantSum([R.ring(0)], [[0]])

    # For each term c*S[w] in the SchubertSum
    for (c, w) in zip(SS.coeffs, SS.schubs)
        # Convert S[w] to dominant basis
        dom_poly = back_schub_poly(w, R)

        # Multiply by the coefficient and add to result
        result = result + c * dom_poly
    end

    return condense(result)
end

# Version without explicit R argument - infer the ring size
function schub2dom(SS::SchubertSum)

    # Find the maximum length among all permutations in SS
    max_len = maximum(length.(SS.schubs))

    # Create appropriate ring
    R = xy_ring(max(max_len - 1, 1))[1]

    return schub2dom(SS, R)
end
