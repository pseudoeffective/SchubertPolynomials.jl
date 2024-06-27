# Constructing semistandard tableaux and Schur polynomials
# David Anderson, February 2024.
# This version uses Memoization

# TO DO -- harmonize schur_poly and ssyt functions with Oscar versions
# at least when there are no y variables



export Tableau, schur_poly, ssyt


struct Tableau
  t::Vector{Vector{Int}}
end



function sstab( la::Vector{Int}, mu=[]; rowmin=false )
# make the superstandard tableau from partition shape la(/mu)
  tab=[]

  if rowmin
    for i in 1:length(la)
      push!(tab, fill(i, la[i]) )
    end

    for i=1:length(mu)
      for j=1:mu[i]
        tab[i][j]=0
      end
    end
    return Tableau(tab)
  end  


  for i in 1:length(la)
    tr=[]
    if i<=length(mu)
      push!(tr,fill(0,mu[i])...)
      j=mu[i]+1
    else
      j=1
    end
    while j<=la[i]
      k=1
      while k<i && tab[i-k][j]>0
        k+=1
      end
      push!(tr,k)
    j+=1
    end
    push!(tab,tr)
  end

  return Tableau(tab)
end



# extract the (skew) shape of a tableau
function shape(tab::Tableau)
  la=[]
  mu=[]

  len=length(tab.t)

  for i in 1:len
    push!(la,length(tab.t[i]))
    j=1
    while tab.t[i][j]==0
      j+=1
    end
    if j>1
      push!(mu,j-1)
    end  
  end

  return la,mu
end


# overload show to display Tableau
function Base.show(io::IO, tab::Tableau)
    println(io)
    for i in 1:length(tab.t)
        for j in 1:length(tab.t[i])
            tt=tab.t[i][j]
            if tt==0
              print(io, 'x', " ")
            else
              print(io, tt, " ")
            end
        end
        println(io)
    end
end

# overload identity for Tableau type
Base.:(==)(tab1::Tableau,tab2::Tableau) = tab1.t==tab2.t


# return the product of binomials
function tab2bin( tab::Tableau, RR::DoublePolyRing; xoffset = 0, yoffset = 0 )
  len = length(tab.t)

  x = RR.x_vars
  y = RR.y_vars

  n = length(x)
  m = length(y)

  bin = RR.ring(1)

  for i=1:len
    for j=1:length( tab.t[i] )
      tt = tab.t[i][j]
      if tt>0
        p=RR.ring(0)
        if tt+xoffset<=n
          p=p+x[tt+xoffset]
        end
        if tt+j-i+yoffset<=m && tt+j-i+yoffset>0
          p=p+y[tt+j-i+yoffset]
        end
        bin = bin*p
      end
    end
  end

  return bin

end


# sum of binomials for a set of tableaux
function ssyt2pol( tabs, RR::DoublePolyRing; xoffset=0, yoffset=0 )

  pol=RR.ring(0)

  for tab in tabs
    pol = pol + tab2bin( tab, RR; xoffset=xoffset, yoffset=yoffset )
  end

  return pol

end



"""
    schur_poly(la, ff, RR=xy_ring(length(la), length(la)+la[1])[1]; mu=[], xoffset=0, yoffset=0, rowmin=false)

Compute the Schur polynomial corresponding to a given (skew) partition `la/mu` and a flag `ff`, in an optionally specified ring `RR`. The polynomial is constructed as an enumerator of semistandard Young tableaux of (skew) shape `la/mu` and bounded by the flagging condition `ff`.

## Arguments
- `la::Vector{Int}`: A partition represented as a vector of integers, specifying the shape of the Young diagram.
- `ff::Union{Int,Vector{Int},Vector{Vector{Int}}}`: A flag specifying bounds on the tableaux. If `ff` is given as a single integer, it bounds the entries of the tableaux.  If `ff` is a vector of integers, it must be of length at least that of `la`; then `ff[i]` bounds the entries in the `i`th row of the tableaux.  If `ff` is a vector of vectors, it is interpreted as a tableaux whose shape is assumed to contain `la`; then the entries of `ff` bound the tableaux entrywise.
- `RR::DoublePolyRing`: An optional argument specifying the double polynomial ring to use for constructing the Schur polynomial. Defaults to a ring constructed based on the size of `la`.
- `mu::Vector{Int}`: An optional argument specifying a subpartition of `la`, for skew Schur polynomials. Defaults to an empty vector, for the straight shape `la`.
- `xoffset::Int`: An optional argument specifying an offset value for the x-variable indices in the polynomial. Defaults to 0.
- `yoffset::Int`: An optional argument specifying an offset value for the y-variable indices in the polynomial. Defaults to 0.
- `rowmin::Bool`: An optional argument specifying whether to use row-minimal tableau, i.e., to require that entries in row `i` be at least `i`.  (This is a nontrivial condition only for skew shapes.) Defaults to `false`.

## Returns
- `ZZMPolyRingElem`: The Schur polynomial as an element of the specified polynomial ring `RR`.

# Examples
```julia-repl
# Specify a partition
julia> la = [2, 1]

# Specify a bound for the x-variables
julia> ff = 3

# Compute the Schur polynomial
julia> poly = schur_poly(la, ff)


### To get the single Schur polynomial, change the coefficient ring
julia> R = xy_ring(3,0)[1];

julia> poly1 = schur_poly(la,ff,R)
x1^2*x2 + x1^2*x3 + x1*x2^2 + 2*x1*x2*x3 + x1*x3^2 + x2^2*x3 + x2*x3^2

```
"""
function schur_poly( la, ff::Vector{Vector{Int}}, RR::DoublePolyRing=xy_ring( length(la) , length(la)+la[1] )[1]; mu = Int[], xoffset=0, yoffset=0, rowmin=false )
  if length(la)==0
    return RR.ring(1)
  end

  if length(RR.y_vars)==0
     return schur_polynomial1_combinat( la, ff, RR, mu=mu, xoffset=xoffset, rowmin=rowmin )
  end

  tbs = ssyt( la, ff, mu=mu, rowmin=rowmin )

  pol = ssyt2pol( tbs, RR; xoffset=xoffset, yoffset=yoffset )

#=
  pol = RR.ring(0)

  for tab in tbs
    pol = pol + tab2bin( tab, RR; xoffset=xoffset, yoffset=yoffset )
  end
=#

  return pol

end

function schur_poly( la, ff::Vector{Int}, RR::DoublePolyRing=xy_ring( length(la) , length(la)+la[1] )[1]; mu = Int[], xoffset=0, yoffset=0, rowmin=false )
  if length(la)==0
    return RR.ring(1)
  end

  return schur_poly( la, Vector{Vector{Int}}([fill(ff[i],la[i]) for i=1:length(la)]), RR; mu = mu, xoffset=xoffset, yoffset=yoffset, rowmin=rowmin )

end

function schur_poly( la, ff::Int, RR::DoublePolyRing=xy_ring( length(la) , length(la)+la[1] )[1]; mu = Int[], xoffset=0, yoffset=0, rowmin=false )
  if length(la)==0
    return RR.ring(1)
  end

  return schur_poly( la, Vector{Int}(fill(ff,length(la))), RR; mu = mu, xoffset=xoffset, yoffset=yoffset, rowmin=rowmin )


end



"""
    ssyt(la::Vector{Int}, ff::Union{Int, Vector{Int}, Vector{Vector{Int}}}; mu::Vector{Int}=[], rowmin::Bool=false) -> Vector{Tableau}

Constructs semistandard Young tableaux on a shape `la`, with given flagging conditions `ff`, optionally skewed by a subshape `mu`.  Method based on Oscar version by Ulrich Thiel and collaborators.

## Arguments
- `la::Vector{Int}`: A partition represented as a nonincreasing vector of integers.

- `ff::Union{Int, Vector{Int}, Vector{Vector{Int}}}`: Flagging conditions for the tableaux. If provided as an integer, it specifies the largest value allowed for an entry. If given as a single vector of integers, it specifies a row flagging condition. If provided as a vector of vectors, it specifies a filling which bounds the tableaux entrywise.

- `mu::Vector{Int}`: An optional keyword argument, giving a subshape of the partition `la`. Defaults to an empty vector for a straight shape.

- `rowmin::Bool`: An optional boolean keyword. When set to `true`, indicates the entries in row `i` must be at least `i`. Condition is redundant for straight-shape tableaux. Defaults to `false`.

## Returns
- `Vector{Tableau}`: A vector of `Tableau` objects, each representing a semistandard Young tableau that satisfies the given shape and flagging conditions.

# Examples
```julia-repl
# Generate SSYTs for the partition [3, 2] with largest entry 3
julia> tabs = ssyt([3, 2], 3)

# Generate SSYTs for the shape [3, 2, 1]/[1] and row-minimal condition
julia> skewtabs = ssyt([3, 2, 1], [3, 3, 3], mu=[1], rowmin=true)
"""
function ssyt( lambda::Vector{Int}, ff::Vector{Vector{Int}}; mu::Vector{Int} = fill(0, length(lambda)), rowmin::Bool=false)

  len = length(lambda)
  tab = sstab(lambda, mu, rowmin=rowmin).t
  m = len
  n = lambda[m]
  mu = vcat(mu, [0 for s=length(mu)+1:len])  # extend mu by 0 if necessary


  tabs_list = Tableau[]
  valid = true
  while true

    for i = 1:len
      for j = mu[i]+1:lambda[i]
        if tab[i][j] > ff[i][j]
          valid = false
          break
        end
      end
    end


    if valid
       push!(tabs_list,Tableau(deepcopy(tab)) )
    end

    #raise one element by 1
    while !(tab[m][n] < ff[m][n] &&
            (n == lambda[m] || tab[m][n] < tab[m][n + 1]) &&
            (m == len || lambda[m + 1] < n || tab[m][n] + 1 < tab[m + 1][n]))
      if n > mu[m]+1
        n -= 1
      elseif m > 1
        m -= 1
        n = lambda[m]
      else
        return tabs_list
      end
    end

    if tab[m][n]>0
       tab[m][n] += 1
    end

    #minimize trailing elements
    if n < lambda[m]
      i = m
      j = n + 1
    else
      i = m + 1
      i<=len ? j = mu[i]+1 : j=1
    end
    while (i <= len && j <= lambda[i])
      if i == 1
        tab[1][j] = tab[1][j - 1]  #if i==1 then j!=1 by initialization
      elseif j == 1
        if rowmin  # ensure tab entry is at least row index
           tab[i][1] = max(i,tab[i - 1][1] + 1) #likewise if j==1 then i!=1
        else
           tab[i][1] = tab[i - 1][1] + 1 #likewise if j==1 then i!=1
        end
      else
        if rowmin  # ensure tab entry is at least row index
           tab[i][j] = max(tab[i][j - 1], tab[i - 1][j] + 1, i) #likewise if j==1 then i!=1
        else
           tab[i][j] = max(tab[i][j - 1], tab[i - 1][j] + 1) #likewise if j==1 then i!=1
        end
      end
      if j < lambda[i]
        j += 1
      else
        i += 1
        i<=len ? j = mu[i]+1 : j=1
      end
    end
    m = len
    n = lambda[len]
  end #while true


end


###
function ssyt( lambda, ff::Vector{Int}; mu = fill(0, length(lambda)), rowmin=false)
# case of flagged tableaux
  bd = Vector{Vector{Int}}([])
  for i=1:length(lambda)
    push!( bd, fill( ff[i], lambda[i] ) )
  end

  ssyt( lambda, bd; mu=mu, rowmin=rowmin )
end


###
function ssyt( lambda, ff::Int=length(la); mu = fill(0, length(lambda)), rowmin=false)
# uniformly bounded, default to number of rows of lambda
  bd = Vector{Vector{Int}}([])
  for i=1:length(lambda)
    push!( bd, fill( ff, lambda[i] ) )
  end

  ssyt( lambda, bd; mu=mu, rowmin=rowmin )
end


################


# this function does much better for single polynomials

function schur_polynomial1_combinat(lambda::Vector{Int}, ff::Vector{Vector{Int}}, R::DoublePolyRing=xy_ring(max(max(ff...)...),0)[1]; mu::Vector{Int}=Int[], xoffset::Int=0, rowmin::Bool=false)
  if isempty(lambda)
    return one(R.ring)
  end

  S = base_ring(R.ring)
  sf = MPolyBuildCtx(R.ring)

  xx = R.x_vars

  #version of the function semistandard_tableaux(shape::Vector{T}, max_val = sum(shape))
  len = length(lambda)
  tab = sstab(lambda, mu, rowmin=rowmin).t
  m = len
  n = lambda[m]
  mu = vcat(mu, [0 for s=length(mu)+1:len])  # extend mu by 0 if necessary


  count = zeros(Int,length(xx))
  valid = true
  while true
    count .= 0
    for i = 1:len
      for j = mu[i]+1:lambda[i]
        if tab[i][j] <= ff[i][j]
          if tab[i][j]+xoffset <= length(xx)
             count[ tab[i][j]+xoffset ] +=1
          end

        else
          valid = false
          break
        end
      end
    end

    if valid
       push_term!(sf, S(1), count )
    end
    #raise one element by 1
    while !(tab[m][n] < ff[m][n] &&
            (n == lambda[m] || tab[m][n] < tab[m][n + 1]) &&
            (m == len || lambda[m + 1] < n || tab[m][n] + 1 < tab[m + 1][n]))
      if n > mu[m]+1
        n -= 1
      elseif m > 1
        m -= 1
        n = lambda[m]
      else
        return finish(sf)
      end
    end

    if tab[m][n]>0
       tab[m][n] += 1
    end

    #minimize trailing elements
    if n < lambda[m]
      i = m
      j = n + 1
    else
      i = m + 1
      i<=len ? j = mu[i]+1 : j=1
    end
    while (i <= len && j <= lambda[i])
      if i == 1
        tab[1][j] = tab[1][j - 1]  #if i==1 then j!=1 by initialization
      elseif j == 1
        if rowmin
           tab[i][1] = max(i,tab[i - 1][1] + 1) #likewise if j==1 then i!=1
        else
           tab[i][1] = tab[i - 1][1] + 1 #likewise if j==1 then i!=1
        end
      else
        if rowmin
           tab[i][j] = max(tab[i][j - 1], tab[i - 1][j] + 1, i) #likewise if j==1 then i!=1
        else
           tab[i][j] = max(tab[i][j - 1], tab[i - 1][j] + 1) #likewise if j==1 then i!=1
        end
      end
      if j < lambda[i]
        j += 1
      else
        i += 1
        i<=len ? j = mu[i]+1 : j=1
      end
    end
    m = len
    n = lambda[len]
  end #while true
end



####### not used
#=
@memoize function ssyt_old(la, bd::Vector{Vector{Int}}; mu = fill(0, length(la)), rowmin=false)
    len = length(la)
    mmu = copy(mu)

    ns = collect( 1:maximum(vcat(bd...)) )
    
    while length(mmu) < len
        push!(mmu, 0)
    end
    
    if length(mmu) > len || mmu[end] > la[end]
        println("error: mu not contained in lambda")
        return nothing
    end
    
    tabs = Vector{Tableau}([])
    
    if len == 0
        return tabs
    end
    
    if len == 1

        if la[1] == mmu[1]
            push!(tabs, Tableau( [fill(0, mmu[1])] ) )
        end

        if la[1] == mmu[1] + 1
            for i in 1:length(ns)
                if ns[i] <= bd[1][la[1]]
                    push!(tabs, Tableau( [append!( copy(fill(0, mmu[1])), ns[i])] ) )
                end
            end
            return tabs
        end

        if la[1] > mmu[1] + 1
            tabs1 = ssyt([la[1] - 1], bd, mu=mmu)
            for tt in tabs1
                for k in 1:length(ns)
                    if ns[k] <= bd[1][la[1]]
                        t2=add_tabi(tt,1,ns[k])
                        if t2!=nothing
                          push!(tabs, t2 )
                        end
                    end
                end
            end
        end

        return tabs
    end
    
    if len > 1
        la1 = la[1:len-1]
        mu1 = mmu[1:len-1]
        
        tabs1 = ssyt(la1, bd, mu=mu1, rowmin=rowmin)
        
        if mmu[end] > 0
            tabs2 = Vector{Tableau}([])
            for tt in tabs1
                push!(tabs2, Tableau( vcat( tt.t, [fill(0, mmu[len])] ) ) )
            end
            tabs1 = tabs2
        end
        
        for i in 1:(la[len] - mmu[len])
            tabs2 = Vector{Tableau}([])
            for tt in tabs1
                for k in 1:length(ns)
                    if ns[k] <= bd[len][mmu[len]+i] && (!rowmin || ns[k]>=len )
                        t2 = add_tabi(tt, len, ns[k])
                        if t2!=nothing
                          push!(tabs2, t2)
                        end
                    end
                end
            end
            tabs1 = tabs2
        end
        tabs = tabs1
    end
        
    return tabs
end


################
# Build tableaux as paths in Young lattice, with memoization


###
# add 'a' to row k of tab, return nothing if fails
function add_tabi(tab::Tableau, k, a)
    len = length(tab.t)

    if k > len + 1
        return nothing
    end

    if k == 1
        if a >= last(tab.t[1])
            temp = [[tab.t[1]; a], tab.t[2:end]...]
            return Tableau(temp)
        else
            return nothing
        end
    end

    if k > 1 && k < len + 1
        lam = length(tab.t[k])
        if lam < length(tab.t[k-1]) && 
           a >= last(tab.t[k]) && 
           a > tab.t[k-1][lam+1]
            temp = [tab.t[1:k-1]..., [tab.t[k]; a], tab.t[k+1:end]...]
            return Tableau(temp)
        else
            return nothing
        end
    end

    if k == len + 1
        if a > first(tab.t[len])
            temp = [tab.t..., [a]]
            return Tableau(temp)
        else
            return nothing
        end
    end
end


=#