# Tools for computing back stable Schubert polynomials in Julia, using BPDs
# David Anderson, June 2024


export back_schub_poly, acoeff


#############

"""
    back_schub_poly(w, R)

Return the back stable Schubert polynomial for the permutation `w`

## Arguments
- `w::Vector{Int}`: A permutation
- `R::DoublePolyRing`: The ambient polynomial ring, with underlying ring `R.ring` and variables `R.x_vars` and `R.y_vars`


## Returns
`ZZMPolyRingElem`: the Schubert polynomial in the ring R

```
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

  return DominantSum(cfs,pars)
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


function bpd2bin_star( bpd::BPD, R::DoublePolyRing=xy_ring( size(bpd.m)[1]-1, size(bpd.m)[2]-1 )[1] )
# product of binomials for bpd
# only schub version now
# can get single polyn by using no y_vars
  local n=size(bpd.m)[1]-1
  bin = R.ring(1)

  x = R.x_vars
  y = R.y_vars

  local aa=length(x)
  local bb=length(y)

  la = dominant_part(bpd)
  la = vcat( la, zeros(Int, n-length(la)) )

  for i=1:n
    for j=la[i]+2:n

      if bpd.m[i,j]==0
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

