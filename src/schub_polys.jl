# Tools for computing Schubert polynomials in Julia, using BPDs
# David Anderson, April 2024


export schub_poly, groth_poly, nschub, ngroth, ddx, pdx

# TO DO: simplify calling methods for schub_poly
#############

"""
    schub_poly(w, R; method)

Return the Schubert polynomial for the permutation `w`

## Arguments
- `w::Vector{Int}`: A permutation
- `R::DoublePolyRing`: The ambient polynomial ring, with underlying ring `R.ring` and variables `R.x_vars` and `R.y_vars`
- `method`: An optional argument specifying the algorithm for computing.
- `memo`: An optional argument for memoization (default=false)

The options for `method` are:
- `method="transition"`, computes via transition formula
- `method="drift"`, computes from flat BPDs by drift class formula
- `method="bpd"`, computes by summing over all BPDs
- `method="dd"`, computes by divided differences from longest permutation 
- `method="auto" (default), computes by transition if there are y variables, or if w has relatively small length, otherwise computes by divided difference

## Returns
`ZZMPolyRingElem`: the Schubert polynomial in the ring R

## Examples
```julia-repl
# Define a single-variable DoublePolyRing
julia> R = xy_ring(5,0)[1];

# Choose a permutation
julia> w = [3,5,1,4,2];

# Compute via transition
julia> pol1 = schub_poly(w, R, method="transition");

# Compute via divided differences
julia> pol2 = schub_poly(w, R, method="dd");

julia> pol1==pol2
true
```

### If there are not enough x-variables to compute by descending induction, the "dd" method throws an error
```julia-repl
julia> R = xy_ring(3,0)[1];

julia> schub_poly(w, R, method="dd")
ERROR: ArgumentError: Not enough x variables for this method

```

### If R is not specified, a double polynomial ring containing the single Schubert polynomial is chosen

```julia-repl
julia> pol3 = schub_poly(w, method="bpd");

julia> R,x,y = xy_ring(4,0);

julia> pol4 = schub_poly(w, R, method="drift");

julia> pol3==pol4
true
```
"""
function schub_poly(w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1) )[1]; method="auto", memo::Bool=false )

  w=trimw(w)
  ws = cutw(w)

  if length(ws)>1
    return( schub_poly(ws[1], R; method=method, memo=memo )*schub_poly(ws[2], R; method=method, memo=memo ) ) 
  end

  if memo
    return schub_trans_memo(w,R)
  end

  if method=="dd"
    return schub_dd(w,R)
  elseif method=="bpd"
    return schub_bpd(w,R)
  elseif method=="drift"
    return schub_drifts(w,R)
  elseif method=="transition"
    return schub_trans(w,R)
  end

  if length(R.y_vars)>0 || len(w)<.3*length(w)^2
    return schub_trans(w,R)
  end

  return schub_dd(w,R)

end


"""
    groth_poly(w, R; method)

Return the Grothendieck polynomial for the permutation `w`

## Arguments
- `w::Vector{Int}`: A permutation
- `R::DoublePolyRing`: The ambient polynomial ring, with underlying ring `R.ring` and variables `R.x_vars` and `R.y_vars`
- `method`: An optional argument specifying the algorithm for computing.

The options for `method` are:
- `method="bpd"` (default), computes by summing over all BPDs
- `method="dd"`, computes by divided differences from longest permutation 

## Returns
`ZZMPolyRingElem`: the Grothendieck polynomial in the ring R

## Examples
```julia-repl
# Define a single-variable DoublePolyRing
julia> R,x,y = xy_ring(5,0);

# Choose a permutation
julia> w = [3,5,1,4,2];

# Compute via divided differences (default)
julia> pol1 = groth_poly(w, R);

# Compute via bpds
julia> pol2 = groth_poly(w, R, method="bpd");

julia> pol1==pol2
true
```
### If R is not specified, a polynomial ring containing the single Grothendieck polynomial is chosen

```julia-repl
julia> pol3 = groth_poly(w, method="bpd");

julia> R,x,y = xy_ring(4,0);

julia> pol4 = groth_poly(w, R, method="dd");

julia> pol3==pol4
true
```
"""
function groth_poly(w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1), 0 )[1]; method="dd" )

  if method=="dd"
    return groth_dd(w,R)
  elseif method=="drift"
    return groth_drifts(w,R)
  end

  return groth_bpd(w,R)

end


## TO DO: add "transition" method for groth_poly
## TO DO: figure out "drift" method for groth_poly


######
# just counting terms

function nschub(w::Vector{Int}, cache::Dict{Int, Int} = Dict{Int,Int}())
# count the number of terms in the Schubert polynomial

    idx = nthperm(w)

    # check cache first
    if haskey(cache, idx)
        return cache[idx]
    end

    if w == collect(1:length(w))  # identity permutation
        cache[idx] = 1
        return 1
    end

    mtx = max_transition(w)

    if isempty(mtx)
        la, ff = vex2flag(w)
        result = round(Int, vex_det(la, ff))
        cache[idx] = result
        return result
    end

    vv = mtx[2]
    wws = mtx[3]

    ss = nschub(vv, cache)

    for ww in wws
        ss += nschub(ww, cache)
    end

    cache[idx] = ss
    return ss
end



function ngroth(w::Vector{Int})
# count the number of terms in the Grothendieck polynomial
# need to implement transition

  w = trimw(w)

  n = length(w)-1

  R = xy_ring(n,0)[1]

  pw = groth_poly(w,R,method="dd")

  return sum(abs.(coefficients(pw)))  

end

######


function bpd2bin( bpd::BPD, R::DoublePolyRing=xy_ring( size(bpd.m)[1]-1, size(bpd.m)[2]-1 )[1]; version="schub"  )
# product of binomials for bpd
# requires DoublePolyRing
# can get single polyn by using no y_vars
  local n=size(bpd.m)[1]-1
  bin = R.ring(1)

  x = R.x_vars
  y = R.y_vars

  local aa=length(x)
  local bb=length(y)

  for i=1:n
    for j=1:n

      if bpd.m[i,j]==0
        p=R.ring(0)
        if i<=aa
          p=p+x[i]
        end
        if j<=bb
          p=p+y[j]
        end
        if version=="groth" && i<=aa && j<=bb
          p=p-x[i]*y[j]
        end     
        bin = bin*p
        if version=="groth"
          bin = -bin
        end
      end

      if version=="groth" && bpd.m[i,j]==3
        p=R.ring(1)
        if i<=aa
          p=p-x[i]
        end
        if j<=bb
          p=p-y[j]
        end
        if i<=aa && j<=bb
          p=p+x[i]*y[j]
        end     
        bin = bin*p
      end

    end
  end

  return bin

end



function schub_bpd( w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1]  )
# compute schubert pol by bpd formula
  bpds=all_bpds(w)

  pol=R.ring(0)

  for bp in bpds
    pol = pol+bpd2bin(bp,R)
  end

  return(pol)

end



#@memoize 
function schub_trans( w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by transition formula

  mxt = max_transition(w)

  if length(mxt)==0
    return schub_drifts(w,R)
  end

  (r,s) = mxt[1]
  vv = mxt[2]
  wws = mxt[3]

  x = R.x_vars
  y = R.y_vars

  pv = schub_trans(vv, R)

  p1 = R.ring(0)

  if r<=length(x)
    p1 = p1+x[r]
  end

  if w[s]<=length(y)
    p1 = p1 + y[w[s]]
  end

  p1 = p1*pv

  for ww in wws
    pw = schub_trans(ww,R)
    p1 = p1+pw
  end

  return p1

end

@memoize function schub_trans_memo( w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by transition formula, memoized

  mxt = max_transition(w)

  if length(mxt)==0
    return schub_drifts(w,R)
  end

  (r,s) = mxt[1]
  vv = mxt[2]
  wws = mxt[3]

  x = R.x_vars
  y = R.y_vars

  pv = schub_trans_memo(vv, R)

  p1 = R.ring(0)

  if r<=length(x)
    p1 = p1+x[r]
  end

  if w[s]<=length(y)
    p1 = p1 + y[w[s]]
  end

  p1 = p1*pv

  for ww in wws
    pw = schub_trans_memo(ww,R)
    p1 = p1+pw
  end

  return p1

end


function schub_dd( w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by divided differences

  w=trimw(w)

  n=length(w)

  if n==0
    return(R.ring(1))
  end

  if length(R.x_vars)<n-1
    throw(ArgumentError("Not enough x variables for this method"))
  end

  i=n-1
  while i>0 && w[i]>w[i+1]
    i=i-1
  end

  if i==0
    return sp0( n, R )
  end

  w = vcat( w[1:i-1], w[i+1], w[i], w[i+2:n] )

  pol1 = schub_dd( w, R )

  return ddx( pol1, i, R )

end



function groth_bpd( w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1]  )
# compute grothendieck pol by bpd formula
  bpds=all_Kbpds(w)

  pol=R.ring(0)

  for bp in bpds
    pol = pol+bpd2bin(bp,R, version="groth")
  end

  return((-1)^len(w)*pol)

end


@memoize function groth_dd( w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute grothendieck pol by divided differences

  n=length(w)
  while n>0 && w[n]==n
    pop!(w)
    n=length(w)
  end

  if n==0
    return(R.ring(1))
  end

  if length(R.x_vars)<n-1
    throw(ArgumentError("Not enough x variables for this method"))
  end

  i=n-1
  while i>0 && w[i]>w[i+1]
    i=i-1
  end

  if i==0
    return gp0( n, R )
  end

  w = vcat( w[1:i-1], w[i+1], w[i], w[i+2:n] )

  pol1 = groth_dd( w, R )

  return pdx( pol1, i, R )

end





#############

function sp0( n::Int, R::DoublePolyRing=xy_ring( n-1, n-1 )[1] )

  x=R.x_vars
  y=R.y_vars

  aa=length(x)
  bb=length(y)

  pol=1

  for i in 1:n-1
    for j in 1:n-i
      p=R.ring(0)
      if i<=aa p=p+x[i] end
      if j<=bb p=p+y[j] end
      pol = pol*p
    end
  end

  return(pol)

end


function gp0( n::Int, R::DoublePolyRing=xy_ring( n-1, n-1 )[1] )

  x=R.x_vars
  y=R.y_vars

  aa=length(x)
  bb=length(y)

  pol=1

  for i in 1:n-1
    for j in 1:n-i
      p=R.ring(0)
      if i<=aa p=p+x[i] end
      if j<=bb p=p+y[j] end
      if i<=aa && j<=bb p=p-x[i]*y[j] end
      pol = pol*p
    end
  end

  return(pol)

end




#####################
# difference operator

function ddx(p::ZZMPolyRingElem, i::Int, R::Union{ZZMPolyRing,DoublePolyRing}=parent(p))
# p is a polynomial in R.x_vars

  if isa(R,ZZMPolyRing)
    RR = R
    x = gens(R)
  elseif isa(R,DoublePolyRing)
    RR = R.ring
    x = R.x_vars
  end

  if i>length(x)
    return(RR(0))
  end

  if i==length(x)
    p1 = evaluate( p, [x[i]], [RR(0)] )
    q=divrem( p-p1, x[i] )[1]
    return q
  end

  p1 = evaluate( p, [x[i],x[i+1]], [x[i+1],x[i]] )

  q=divrem( p-p1, x[i]-x[i+1] )[1]

  return q

end


# isobaric difference operator
function pdx(p::ZZMPolyRingElem, i::Int, R::Union{ZZMPolyRing,DoublePolyRing}=parent(p))
# p is a polynomial in R.x_vars

  if isa(R,ZZMPolyRing)
    RR = R
    x = gens(R)
  elseif isa(R,DoublePolyRing)
    RR = R.ring
    x = R.x_vars
  end

  if i>length(x)
    return(p)
  end

  if i==length(x)
    p1 = evaluate( p, [x[i]], [RR(0)] )
    q=divrem( p-(1-x[i])*p1, x[i] )[1]
    return q
  end

  p1 = evaluate( p, [x[i],x[i+1]], [x[i+1],x[i]] )

  q=divrem( (1-x[i+1])*p-(1-x[i])*p1, x[i]-x[i+1] )[1]

  return q

end


# for counting vexillary terms
function vex_det( la::Vector{Int}, ff::Vector{Int} )

  n=length(la)

  if n==0 return 1 end

  A = zeros(BigFloat,n,n) # may need this depending on float tolerance
#  A = zeros(n,n)

  for i in 1:n
    for j in 1:n
      if la[i]-i+j>=0
        A[i,j]=binomial( la[i]-i+j+ff[i]-1, ff[i]-1 )
      end
    end
  end

  return det(A)

end




#=
# temporary alternative, not used
function schub_2drifts( w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by drift class formula

  fbpds = flat_bpds(w)

  pol = R.ring(0)

  for b in fbpds
    d=bpd2drift(b)
    pol = pol+drift_poly(d,R)
  end

  return pol

end

=#