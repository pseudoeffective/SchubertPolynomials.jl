# Tools for multiplying Schubert polynomials in Julia
# David Anderson, March 2025.


export mult_2schub, lrc, mult_schub, expand_schub, maxvar



#########

"""
    mult_2schub( uu, vv, rnk; double=false, ring=schub_ring(rnk-1, rnk-1) )

Multiply two Schubert classes

## Arguments
- `uu::Vector{Int}`: a permutation
- `vv::Vector{Int}`: a permutation
- `rnk::Int`: an integer specifying the ambient flag variety, Fl(rnk).  Defaults to the larger of the lengths of uu and vv.

## Keywords
- `double::Bool`: when `true`, computes the equivariant (double) product, using the y-variables of `ring`.  Default `false`.
- `ring::MPolyRing`: the ambient polynomial ring (build with `schub_ring`).  Defaults to `schub_ring(rnk-1, rnk-1)`.

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

# Equivariant product in a ring with y-variables
julia> R = schub_ring(3,3);

julia> mult_2schub( [2,3,1], [3,1,2], 4; ring=R, double=true )
(y1 - y3)*S[3, 2, 1] + S[4, 2, 1, 3]
```
"""
function mult_2schub( uu::Vector{Int}, vv::Vector{Int}, rnk::Int=maximum([length(uu),length(vv)]);
                      double::Bool=false, ring::MPolyRing=schub_ring(rnk, rnk) )

  len(uu)>len(vv) && return mult_2schub( vv, uu, rnk; double=double, ring=ring )

  f = schub_poly(uu; ring=ring, double=double)

  return mult_schub( f, vv, rnk; double=double, ring=ring )

end



# get the schubert structure constant, default non-equivariant
function lrc( uu::Vector{Int}, vv::Vector{Int}, ww::Vector{Int}, rnk::Int=maximum([length(uu),length(vv),length(ww)]);
              double::Bool=false, ring::MPolyRing=schub_ring(rnk, rnk) )

  ss = mult_2schub( uu, vv, rnk; double=double, ring=ring )

  return coeff( ss, ww )

end



# multiply a polynomial by a schubert sum, default non-equivariant
function mult_schub(f::Union{Int,ZZMPolyRingElem}, ss::SchubertSum, rnk::Int=maximum( length.(ss.schubs) );
                    double::Bool=false, ring::MPolyRing=schub_ring(rnk, rnk) )

  i = maxvar(f)
  cfs = copy(ss.coeffs)

  if i==0
    cfs = [f*c for c in cfs]
    return SchubertSum( cfs, ss.schubs )
  end

  xx = extract_vars(ring; varname=:x)

  ii = findfirst( isequal("x$(i)"), string.(xx) )

  f1 = evaluate( f, [xx[ii]], [0] )

  f2 = divrem( f-f1, xx[ii] )[1]

  ss1 = mult_schub( f1, ss, rnk; double=double, ring=ring )  # maxvar less than i

  ss2 = mult_schub( f2, ss, rnk; double=double, ring=ring )  # lower degree in xi

  ssmk = monk(ss2,i,rnk,ring; double=double) # put the xi back in

  return condense( ss1+ssmk )

end

# to multiply by a single w, default non-equivariant
function mult_schub(f::Union{Int,ZZMPolyRingElem}, w::Vector{Int}, rnk::Int=length(w);
                    double::Bool=false, ring::MPolyRing=schub_ring(rnk, rnk) )
  ss = SchubertSum(w)
  return mult_schub(f,ss,rnk; double=double, ring=ring)
end



# expand a polynomial in Schubert basis, default non-equivariant
function expand_schub( f::Union{Int,ZZMPolyRingElem}, rnk::Int=maxvar( f );
                       double::Bool=false, ring::MPolyRing=parent(f) )

  return mult_schub( f, [1], rnk; double=double, ring=ring )

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


function monk(ss::SchubertSum, i::Int, rnk::Int, R::MPolyRing=schub_ring(rnk, rnk); double::Bool=false)

  if i>rnk || i<0
    return(ss)
  end

  ll = length(ss)

  if double
    yy = extract_vars(R; varname=:y)
  else
    yy = ZZMPolyRingElem[]
  end


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
        if vvi===nothing
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

            if vvi===nothing
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

       if vvi===nothing
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

        if vvi===nothing
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

           if vvi===nothing
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

        if vvi===nothing
          push!(sch2,vv)
          push!(cfs2, ss.coeffs[k])
        else  
          cfs2[vvi]=cfs2[vvi]+ss.coeffs[k]
        end
     end

  end

  return SchubertSum( cfs2, sch2 )

end



