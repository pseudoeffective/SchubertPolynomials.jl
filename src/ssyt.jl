# Constructing semistandard tableaux and Schur polynomials
# David Anderson, June 2024.
# This version based on OSCAR


	#- - - - - - - - - - - - - - - - - - - -#
	#	Type, Constructors, & Iterators		#
	#- - - - - - - - - - - - - - - - - - - -#

export Tableau, ssyt

struct Tableau
  t::Vector{Vector{Int}}
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



		# * * Base Methods * * #

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



		# * * SSYT Iterator* * * #
	
# *Not yet implemented as an interator.
	

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






	#- - - - - - - - - - - -#
	#	Tableau Methods		#
	#- - - - - - - - - - - -#
	
		# * * Schur Polynomials * * #

export schur_poly

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

  return pol

end

###
function schur_poly( la, ff::Vector{Int}, RR::DoublePolyRing=xy_ring( length(la) , length(la)+la[1] )[1]; mu = Int[], xoffset=0, yoffset=0, rowmin=false )
  if length(la)==0
    return RR.ring(1)
  end

  return schur_poly( la, Vector{Vector{Int}}([fill(ff[i],la[i]) for i=1:length(la)]), RR; mu = mu, xoffset=xoffset, yoffset=yoffset, rowmin=rowmin )

end

###
function schur_poly( la, ff::Int, RR::DoublePolyRing=xy_ring( length(la) , length(la)+la[1] )[1]; mu = Int[], xoffset=0, yoffset=0, rowmin=false )
  if length(la)==0
    return RR.ring(1)
  end

  return schur_poly( la, Vector{Int}(fill(ff,length(la))), RR; mu = mu, xoffset=xoffset, yoffset=yoffset, rowmin=rowmin )
end




# faster polynomial constructor for single polynomials, based on OSCAR version

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





	# * * Tableaux operations * * #

# returns the transpose of the Tableau T
function Base.transpose(T::Tableau)
	t = []
	if length(T.t) > 0
		for i in 1:length(T.t[1]) push!(t,[]) end
		for row in T.t for j in 1:length(row)
			push!(t[j],row[j])
		end end
	end
	return Tableau(t)
end



	# * * Insertion Algorithms * * #

export rsk_insert, rsk

# RSK insert the biletter (q,p) into the bitableau (P,Q).
function rsk_insert((P,Q)::Tuple{Tableau,Tableau},(q,p)::Tuple{Integer,Integer})
	local PP = deepcopy(P)
	local QQ = deepcopy(Q)
	local pp = p
	
	local i = 1
	while true
		if i > length(PP.t)
			push!(PP.t,[pp])
			push!(QQ.t,[q])
			break
		elseif pp >= last(PP.t[i])
			push!(PP.t[i],pp)
			push!(QQ.t[i],q)
			break
		end
		
		for j in 1:length(PP.t[i])
			if pp < PP.t[i][j]
				local temp = PP.t[i][j]
				PP.t[i][j] = pp
				pp = temp
				break
			end
		end
		
		i += 1
	end
	
	return (PP,QQ)
end

# RSK insert the letter p into the tableau P
function rsk_insert(P::Tableau,p::Integer)
	return rsk_insert((P,P),(1,p))[1]
end

# RSK insert a biword into the bitableau (P,Q).
function rsk_insert((P,Q)::Tuple{Tableau,Tableau},biword::Vector{<:Tuple{Vararg{Integer}}})
	local PP = deepcopy(P)
	local QQ = deepcopy(Q)
	for (q,p) in sort(biword)
		(PP,QQ) = rsk_insert((PP,QQ),(q,p))
	end
	return (PP,QQ)
end

# RSK insert a word into the tableau P.
function rsk_insert(P::Tableau,word::Vector{<:Integer})
	local PP = deepcopy(P)
	for p in word
		PP = rsk_insert(PP,p)
	end
	return PP
end

# sends a biword to its RSK tableau (P,Q).
function rsk(biword::Vector{<:Tuple{Vararg{Integer}}})
	local P = Tableau([])
	local Q = Tableau([])
	for (q,p) in sort(biword)
		(P,Q) = rsk_insert((P,Q),(q,p))
	end
	return (P,Q)
end

# sends a word to its RSK recording insertion P
function rsk(word::Vector{<:Integer})
	return rsk([(i,word[i]) for i in 1:length(word)])
end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> T0 = Tableau([[1,1,3],[2,4],[5,5]])
	1 1 3
	2 4
	5 5
	
	julia> T1 = rsk_insert(T0,2)
	1 1 2
	2 3
	4 5
	5
	
	julia> T2 = rsk_insert(T1,[4,3,4])
	1 1 2 3 4
	2 3 4
	4 5
	5
	
	julia> biword = [(1,2),(1,5),(2,3),(2,5),(3,1),(3,1),(3,4),(3,5),(4,3),(5,4)];
	
	julia> (P,Q) = rsk(biword)
	(
	1 1 3 4
	2 3 4 5
	5 5
	,
	1 1 2 3
	2 3 3 5
	3 4
	)
	
	julia> rsk_insert((P,Q),(6,2))
	(
	1 1 2 4
	2 3 3 5
	4 5
	5
	,
	1 1 2 3
	2 3 3 5
	3 4
	6
	)
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



export edelman_greene_insert, edelman_greene

# Edelman-Greene insert the biletter (q,p) into the bitableau (P,Q)
function edelman_greene_insert((P,Q)::Tuple{Tableau,Tableau},(q,p)::Tuple{Integer,Integer})
	local Pt = transpose(P)
	local Qt = transpose(Q)
	local pp = p
	
	local j = 1
	while true
		if j > length(Pt.t)
			push!(Pt.t,[pp])
			push!(Qt.t,[q])
			break
		elseif pp > last(Pt.t[j])
			push!(Pt.t[j],pp)
			push!(Qt.t[j],q)
			break
		end
		
		for i in 1:length(Pt.t[j])
			if pp < Pt.t[j][i]
				if Pt.t[j][i]==pp+1 && i>1 if Pt.t[j][i-1]==pp
					pp += 1
					break
				end end
				
				local temp = Pt.t[j][i]
				Pt.t[j][i] = pp
				pp = temp
				break
			end
		end
		
		j += 1
	end
	
	return (transpose(Pt),transpose(Qt))
end

# Edelman-Greene insert a biword into the bitableau (P,Q)
function edelman_greene_insert((P,Q)::Tuple{Tableau,Tableau},biword::Vector{<:Tuple{Vararg{Integer}}})
	local PP = deepcopy(P)
	local QQ = deepcopy(Q)
	for (q,p) in biword
		(PP,QQ) = edelman_greene_insert((PP,QQ),(q,p))
	end
	return (PP,QQ)
end

function edelman_greene(biword::Vector{<:Tuple{Vararg{Integer}}})
	local P = Tableau([])
	local Q = Tableau([])
	for (q,p) in biword
		(P,Q) = edelman_greene_insert((P,Q),(q,p))
	end
	return (P,Q)
end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> biword = [(1,4),(1,2),(2,5),(2,3),(2,2),(3,4),(3,3)];
	
	julia> (P,Q) = edelman_greene_insert(biword)
	(
	2 3 4
	3 4 5
	4
	,
	1 1 2
	2 2 3
	3
	)
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#
