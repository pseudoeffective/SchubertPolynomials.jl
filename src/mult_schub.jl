# Tools for multiplying Schubert polynomials in Julia
# David Anderson, June 2024.


export mult_2schub, lrc, mult_schub, expand_schub


#########

"""
    mult_2schub( uu, vv, rnk, R )

Multiply two Schubert classes

## Arguments
- `uu::Vector{Int}`: a permutation
- `uu::Vector{Int}`: a permutation
- `rnk::Int`: an integer specifying the ambient flag variety, Fl(rnk).  Defaults to the larger of the lengths of uu and vv.
- `R::DoublePolyRing`: the ambient polynomial ring, with underlying ring `R.ring` and variables `R.x_vars` and `R.y_vars`.  Defaults to a single-variable ring in rnk-1 x variables.

## Returns
`SchubertSum`: a sum of Schubert classes

## Examples
```julia-repl
# Choose two permutations
julia> u = [2,1];
julia> v = [1,3,2];

# Compute the (non-equivariant) product
julia> mult_2schub( u, v, 3 )
S[3, 1, 2] + S[2, 3, 1]

# Define a double-variable DoublePolyRing
julia> R = xy_ring(3,3)[1];

julia> mult_2schub( [2,3,1], [3,1,2], 4, R )
(y1 - y3)*S[3, 2, 1] + S[4, 2, 1, 3]
```
"""
function mult_2schub( uu::Vector{Int}, vv::Vector{Int}, rnk::Int=maximum([length(uu),length(vv)]), R=xy_ring(rnk-1)[1] )

  len(uu)>len(vv) && return mult_2schub( vv, uu, rnk, R )

  f = schub_poly(uu,R)

  return mult_schub( f, vv, rnk, R )

end



# get the schubert structure constant, default non-equivariant
function lrc( uu::Vector{Int}, vv::Vector{Int}, ww::Vector{Int}, rnk::Int=maximum([length(uu),length(vv),length(ww)]), R=xy_ring(rnk-1)[1] )

  ss = mult_2schub( uu, vv, rnk, R )

  return coeff( ss, ww )

end



# multiply a polynomial by a schubert sum, default non-equivariant
function mult_schub(f::Union{Int,ZZMPolyRingElem}, ss::SchubertSum, rnk::Int=maximum( length.(ss.schubs) ), R=xy_ring(rnk-1)[1] )

  i = maxvar(f)
  cfs = copy(ss.coeffs)

  if i==0
    cfs = [f*c for c in cfs]
    return SchubertSum( cfs, ss.schubs )
  end

  xx = R.x_vars

  ii = findfirst( isequal("x$(i)"), string.(xx) )

  f1 = evaluate( f, [xx[ii]], [0] )

  f2 = divrem( f-f1, xx[ii] )[1]

  ss1 = mult_schub( f1, ss, rnk, R )  # maxvar less than i

  ss2 = mult_schub( f2, ss, rnk, R )  # lower degree in xi

  ssmk = monk(ss2,i,rnk,R) # put the xi back in

  return condense( ss1+ssmk )

end

# to multiply by a single w, default non-equivariant
function mult_schub(f::Union{Int,ZZMPolyRingElem}, w::Vector{Int}, rnk::Int=length(w), R=xy_ring(rnk-1)[1] )
  ss = SchubertSum(w)
  return mult_schub(f,ss,rnk,R)
end



# expand a polynomial in Schubert basis, default non-equivariant
function expand_schub( f::Union{Int,ZZMPolyRingElem}, rnk::Int=maxvar( f ), R=xy_ring(rnk)[1] )

  return mult_schub( f, [1], rnk, R )

end

#############



function maxvar(f::ZZMPolyRingElem)  # get max x-variable index

  vrs = vars(f)
  maxv = 0

  for v in vrs
     st = string(v)
     if st[1]=='x'
       ndx = parse(Int, st[2:length(st)])

       maxv = max(maxv,ndx)
     end
  end

  return maxv

end

function maxvar(f::Int)
  return 0
end


function monk(ss::SchubertSum, i::Int, rnk::Int, R=xy_ring(rnk)[1])

  if i>rnk || i<0
    return(ss)
  end

  ll = length(ss)
  yy = R.y_vars

  cfs2 = Vector{Union{Int,ZZMPolyRingElem}}([])  # new list of coeffs
  sch2 = Vector{Vector{Int}}([])  # new list of schubs

  for k in 1:ll
     w=ss.schubs[k]
     n=length(w)

     # compute equivariant coefficient
     i <= n ? wi=w[i] : wi=i
   
     wi <= length(yy) ? eq_coeff = -yy[wi] : eq_coeff = 0

     if eq_coeff !=0
        vvi = findfirst(isequal(w),sch2)
        if vvi==nothing
          push!(sch2, w)
          push!(cfs2, eq_coeff*ss.coeffs[k])
        else
          cfs2[vvi]=cfs2[vvi]+eq_coeff*ss.coeffs[k]
        end
     end


     # now do classical coeff    

     if i<=n+1

       last=0
       i<n ? tl = w[i+1:n] : tl = []

       for j in i-1:-1:1

         if last < w[j] && w[j] < wi  # ensures length goes up by 1
            last = w[j]

            vv = vcat( w[1:j-1], wi, w[j+1:i-1], last, tl )  # multiply by transposition for monk

            vvi = findfirst(isequal(vv),sch2)  # see if vv is in the list

            if vvi==nothing
              push!(sch2,vv)
              push!(cfs2, -ss.coeffs[k])
            else  
              cfs2[vvi]=cfs2[vvi]-ss.coeffs[k]
            end
         end
       end
        

     else
       vv = vcat( w, n+1:i-2, i,i-1 )

       vvi = findfirst(isequal(vv),sch2)  # see if vv is in the list

       if vvi==nothing
          push!(sch2,vv)
          push!(cfs2, -ss.coeffs[k])
       else  
          cfs2[vvi]=cfs2[vvi]-ss.coeffs[k]
       end

     end

     if i>=n+1

        i==rnk && continue

        vv = vcat( w, n+1:i-1, i+1,i )
        vvi = findfirst(isequal(vv),sch2)  # see if vv is in the list

        if vvi==nothing
          push!(sch2,vv)
          push!(cfs2, ss.coeffs[k])
        else  
          cfs2[vvi]=cfs2[vvi]+ss.coeffs[k]
        end

        continue

     end

     last = Inf

     for j in i+1:n

       if wi<w[j] && w[j]<last
           last = w[j]

           vv = vcat( w[1:i-1], last, w[i+1:j-1], wi, w[j+1:n] )  # multiply by transposition for monk

           vvi = findfirst(isequal(vv),sch2)  # see if vv is in the list

           if vvi==nothing
             push!(sch2,vv)
             push!(cfs2, ss.coeffs[k])
           else  
             cfs2[vvi]=cfs2[vvi]+ss.coeffs[k]
           end
       end
     end

     if last==Inf && n<rnk

        vv = vcat( w[1:i-1], n+1, w[i+1:n], wi )

        vvi = findfirst(isequal(vv),sch2)  # see if vv is in the list

        if vvi==nothing
          push!(sch2,vv)
          push!(cfs2, ss.coeffs[k])
        else  
          cfs2[vvi]=cfs2[vvi]+ss.coeffs[k]
        end
     end

  end

  return SchubertSum( cfs2, sch2 )

end



