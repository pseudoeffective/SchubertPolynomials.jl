# Tools for computing Schubert polynomials in Julia, using BPDs
# David Anderson, April 2024


export schub_poly, groth_poly, nschub, nschub_legacy, ngroth, ddx, pdx

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

"""
    nschub_legacy(w, cache)

Legacy implementation of nschub using max_transition and vex_det.
"""
function nschub_legacy(w::Vector{Int}, cache::Dict{BigInt, BigInt} = Dict{BigInt,BigInt}())

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
        result = round(BigInt, vex_det(la, ff))
        cache[idx] = result
        return result
    end

    vv = mtx[2]
    wws = mtx[3]

    ss = nschub_legacy(vv, cache)

    for ww in wws
        ss += nschub_legacy(ww, cache)
    end

    cache[idx] = ss
    return ss
end


######
# Transition DFS for nschub (bit-packed)

function _nschub_transition(code::T, n::Int, cache::Dict{T, BigInt}) where T <: PackedPerm
    haskey(cache, code) && return cache[code]

    r, s = find_transition_rs(code, n)
    if r == 0
        # Dominant permutation: return 1
        cache[code] = BigInt(1)
        return BigInt(1)
    end

    v = swap_positions(code, r, s)
    result = _nschub_transition(v, n, cache)

    vr = get_val(v, r)
    for i in 1:r-1
        if get_val(v, i) < vr && is_bruhat_cover_up(v, i, r)
            wprime = swap_positions(v, i, r)
            result += _nschub_transition(wprime, n, cache)
        end
    end

    cache[code] = result
    return result
end


######
# Cotransition DFS for nschub (bit-packed)

function _nschub_cotrans(code::T, n::Int, cache::Dict{T, BigInt}) where T <: PackedPerm
    haskey(cache, code) && return cache[code]

    ci = find_cotrans_index(code, n)
    if ci == 0
        # w0: return 1
        cache[code] = BigInt(1)
        return BigInt(1)
    end

    pi_ci = get_val(code, ci)
    result = BigInt(0)

    for a in 1:n
        for b in a+1:n
            if is_bruhat_cover_up(code, a, b)
                sigma_ci = ci == a ? get_val(code, b) : (ci == b ? get_val(code, a) : pi_ci)
                if sigma_ci != pi_ci
                    cover = swap_positions(code, a, b)
                    result += _nschub_cotrans(cover, n, cache)
                end
            end
        end
    end

    cache[code] = result
    return result
end


######
# CRT-based modular arithmetic versions
# Run transition/cotransition mod a prime, storing UInt64 residues instead of BigInt

# Large primes that fit in UInt64 (each ~63 bits, ~19 decimal digits)
const _CRT_PRIMES = UInt64[
    9223372036854775783,  # largest prime < 2^63
    9223372036854775643,
    9223372036854775549,
    9223372036854775507,
    9223372036854775433,
    9223372036854775421,
]

"""
Reconstruct a BigInt from residues mod the CRT primes using the Chinese Remainder Theorem.
"""
function _crt_reconstruct(residues::Vector{UInt64}, primes::Vector{UInt64})
    # Standard CRT reconstruction
    M = prod(BigInt.(primes))
    result = BigInt(0)
    for i in eachindex(primes)
        p = BigInt(primes[i])
        Mi = div(M, p)
        # Compute modular inverse of Mi mod p using extended gcd
        yi = invmod(Mi, p)
        result += BigInt(residues[i]) * Mi * yi
    end
    return mod(result, M)
end

function _nschub_transition_mod(code::T, n::Int, cache::Dict{T, UInt64}, p::UInt64) where T <: PackedPerm
    haskey(cache, code) && return cache[code]

    r, s = find_transition_rs(code, n)
    if r == 0
        cache[code] = UInt64(1)
        return UInt64(1)
    end

    v = swap_positions(code, r, s)
    result = _nschub_transition_mod(v, n, cache, p)

    vr = get_val(v, r)
    for i in 1:r-1
        if get_val(v, i) < vr && is_bruhat_cover_up(v, i, r)
            wprime = swap_positions(v, i, r)
            result = mod(result + _nschub_transition_mod(wprime, n, cache, p), p)
        end
    end

    cache[code] = result
    return result
end

function _nschub_cotrans_mod(code::T, n::Int, cache::Dict{T, UInt64}, p::UInt64) where T <: PackedPerm
    haskey(cache, code) && return cache[code]

    ci = find_cotrans_index(code, n)
    if ci == 0
        cache[code] = UInt64(1)
        return UInt64(1)
    end

    pi_ci = get_val(code, ci)
    result = UInt64(0)

    for a in 1:n
        for b in a+1:n
            if is_bruhat_cover_up(code, a, b)
                sigma_ci = ci == a ? get_val(code, b) : (ci == b ? get_val(code, a) : pi_ci)
                if sigma_ci != pi_ci
                    cover = swap_positions(code, a, b)
                    result = mod(result + _nschub_cotrans_mod(cover, n, cache, p), p)
                end
            end
        end
    end

    cache[code] = result
    return result
end

function _nschub_cotrans_bfs_mod(code::T, n::Int, p::UInt64) where T <: PackedPerm
    max_ell = n * (n - 1) ÷ 2
    ell = perm_length(code, n)
    ell == max_ell && return UInt64(1)

    layer = Dict{T, UInt64}(code => UInt64(1))

    for level in ell:max_ell-1
        next = Dict{T, UInt64}()
        for (pi_code, pi_val) in layer
            ci = find_cotrans_index(pi_code, n)
            ci == 0 && continue
            pi_ci = get_val(pi_code, ci)
            for a in 1:n
                for b in a+1:n
                    if is_bruhat_cover_up(pi_code, a, b)
                        sigma_ci = ci == a ? get_val(pi_code, b) : (ci == b ? get_val(pi_code, a) : pi_ci)
                        if sigma_ci != pi_ci
                            sigma_code = swap_positions(pi_code, a, b)
                            next[sigma_code] = mod(get(next, sigma_code, UInt64(0)) + pi_val, p)
                        end
                    end
                end
            end
        end
        layer = next
    end

    result = UInt64(0)
    for v in values(layer)
        result = mod(result + v, p)
    end
    return result
end

"""
    _nschub_crt(w, n, code, method)

Run the nschub computation using CRT: compute residues mod several large primes,
then reconstruct the BigInt result. Uses ~6 primes for ~114 digits of range.
"""
function _nschub_crt(code::T, n::Int; method::Symbol=:transition) where T <: PackedPerm
    primes = _CRT_PRIMES
    residues = Vector{UInt64}(undef, length(primes))

    for (idx, p) in enumerate(primes)
        if method == :cotransition
            residues[idx] = _nschub_cotrans_mod(code, n, Dict{T, UInt64}(), p)
        elseif method == :cotransition_bfs
            residues[idx] = _nschub_cotrans_bfs_mod(code, n, p)
        else
            residues[idx] = _nschub_transition_mod(code, n, Dict{T, UInt64}(), p)
        end
    end

    return _crt_reconstruct(residues, primes)
end


######
# Cotransition BFS (sort-reduce) for nschub (bit-packed)

function _nschub_cotrans_bfs(code::T, n::Int) where T <: PackedPerm
    max_ell = n * (n - 1) ÷ 2
    ell = perm_length(code, n)
    ell == max_ell && return BigInt(1)

    layer = Dict{T, BigInt}(code => BigInt(1))

    for level in ell:max_ell-1
        next = Dict{T, BigInt}()
        for (pi_code, pi_val) in layer
            ci = find_cotrans_index(pi_code, n)
            ci == 0 && continue
            pi_ci = get_val(pi_code, ci)
            for a in 1:n
                for b in a+1:n
                    if is_bruhat_cover_up(pi_code, a, b)
                        sigma_ci = ci == a ? get_val(pi_code, b) : (ci == b ? get_val(pi_code, a) : pi_ci)
                        if sigma_ci != pi_ci
                            sigma_code = swap_positions(pi_code, a, b)
                            next[sigma_code] = get(next, sigma_code, BigInt(0)) + pi_val
                        end
                    end
                end
            end
        end
        layer = next
    end

    return sum(values(layer); init=BigInt(0))
end


######
# Public API

"""
    nschub(w; method=:auto)

Count the number of monomials in the Schubert polynomial S_w.

## Arguments
- `w::Vector{Int}`: A permutation
- `method::Symbol`: Algorithm to use. Options:
  - `:auto` or `:transition` — transition DFS with bit-packed memoization (default)
  - `:cotransition` — cotransition DFS with bit-packed memoization
  - `:cotransition_bfs` — cotransition BFS (sort-reduce, no recursion)
  - `:crt` — transition DFS with CRT-based modular arithmetic (lower memory)
  - `:crt_cotransition` — cotransition DFS with CRT-based modular arithmetic
  - `:crt_cotransition_bfs` — cotransition BFS with CRT-based modular arithmetic
  - `:legacy` — original implementation using max_transition/vex_det

The CRT methods store `UInt64` residues mod several large primes instead of `BigInt`,
reducing memory usage. Results are reconstructed via Chinese Remainder Theorem,
supporting values up to ~114 decimal digits.

## Returns
`BigInt`: the number of monomials in the Schubert polynomial
"""
function nschub(w::Vector{Int}; method::Symbol=:auto)
    w = trimw(w)
    n = length(w)
    n == 0 && return BigInt(1)

    # Factor through cutw for disconnected permutations
    ws = cutw(w)
    if length(ws) > 1
        return prod(nschub(wi; method=method) for wi in ws)
    end

    if method == :legacy
        return nschub_legacy(w)
    end

    # Fall back to legacy for n > 25 (exceeds UInt128 packing)
    if n > 25
        return nschub_legacy(w)
    end

    # CRT methods
    if method in (:crt, :crt_cotransition, :crt_cotransition_bfs)
        inner = method == :crt ? :transition :
                method == :crt_cotransition ? :cotransition : :cotransition_bfs
        if n <= 16
            return _nschub_crt(pack_perm64(w), n; method=inner)
        else
            return _nschub_crt(pack_perm128(w), n; method=inner)
        end
    end

    if n <= 16
        code = pack_perm64(w)
        if method == :cotransition
            return _nschub_cotrans(code, n, Dict{UInt64, BigInt}())
        elseif method == :cotransition_bfs
            return _nschub_cotrans_bfs(code, n)
        else
            return _nschub_transition(code, n, Dict{UInt64, BigInt}())
        end
    else
        code = pack_perm128(w)
        if method == :cotransition
            return _nschub_cotrans(code, n, Dict{UInt128, BigInt}())
        elseif method == :cotransition_bfs
            return _nschub_cotrans_bfs(code, n)
        else
            return _nschub_transition(code, n, Dict{UInt128, BigInt}())
        end
    end
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


function bpd2bin( bpd::BPD, R::DoublePolyRing=xy_ring( size(bpd.mtx)[1]-1, size(bpd.mtx)[2]-1 )[1]; version="schub"  )
# product of binomials for bpd
# requires DoublePolyRing
# can get single polyn by using no y_vars
  local n=size(bpd.mtx)[1]-1
  bin = R.ring(1)

  x = R.x_vars
  y = R.y_vars

  local aa=length(x)
  local bb=length(y)

  for i=1:n
    for j=1:n

      if bpd.mtx[i,j]==0
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

      if version=="groth" && bpd.mtx[i,j]==3
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