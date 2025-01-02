
#==============================================================================================================================+
	
	permtools.jl
	David Anderson, April 2024.

	Tools for manipulating permutations in Julia

+==============================================================================================================================#

export len, descents, sij, dominant_transition, max_transition, trimw

#==============================================================================================================================#

# permutation tools

function len( w::Vector{Int} )
# coxeter length of a permutation (or word)
  n = length(w)
  a = 0

  for i in 1:n-1
    for j in i+1:n
      if w[i]>w[j]
        a=a+1
      end
    end
  end

  return a
end



function descents(w::Vector{Int})
# compute the descent set
  n=length(w)

  ds=[]

  for i=1:n-1
    w[i]>w[i+1] && push!(ds,i)
  end

  return ds

end


function cutw( w::Vector{Int} )
# split a permutation w = w1 w2 into two which permute distinct indices, if possible

  n=length(w)

  if n<3
    return [w]
  end

  for i in 1:n-3
    if Set(w[n-i:n])==Set(n-i:n)
      w1 = w[1:n-i-1]
      w2 = vcat(1:n-i-1,w[n-i:n])
      if w1==collect(1:n-i-1)
        return [w2]
      end
      return [w1,w2]
    end
  end

  return [w]

end


function findlast132( w::Vector{Int} )
# find the lex-last subsequence of pattern 132, return its indices or [] if none

  n=length(w)

  if n<3
    return []
  end

  if n==3
    w==[1,3,2] && return [1,2,3]
    return []
  end

  for i in n-2:-1:1
    for j in n-1:-1:i+1
      for k in n:-1:j+1
        if invperm(sortperm( w[[i,j,k]] ) )==[1,3,2]
          return [i,j,k]
        end
      end
    end
  end

  return []

end


function findlast2143( w::Vector{Int} )
# find the lex-last subsequence of pattern 2143, return its indices or [] if none

  n=length(w)

  if n<4
    return []
  end

  if n==4
    w==[2,1,4,3] && return [1,2,3,4]
    return []
  end

  for i in n-3:-1:1
    for j in n-2:-1:i+1
      for k in n-1:-1:j+1
        for l in n:-1:k+1
          if invperm(sortperm( w[[i,j,k,l]] ) )==[2,1,4,3]
            return [i,j,k,l]
          end
        end
      end
    end
  end

  return []

end


function maxinversion(w::Vector{Int})

  n=length(w)

  if n<2 return [] end
  if n==2
    w==[2,1] && return [1,2]
    return []
  end

  for i in n-1:-1:1
    for j in n:-1:i+1
      if w[i]>w[j]
        return [i,j]
      end
    end
  end

  return []

end



function sij( w::Vector{Int}, i::Int, j::Int )
# transposition w t_ij

  n=length(w)
  if i==j return w end

  if i>j
    a=j
    b=i
    return sij(w,a,b)
  end

  wt = vcat( w[1:i-1], w[j], w[i+1:j-1], w[i], w[j+1:n] )

  return wt

end



function dominant_transition( w::Vector{Int} )
# do step in Fan-Guo-Sun transition tree

  pos = findlast132(w)

  if length(pos)==0
    return []
  end

  r,s = pos[2],pos[3]

  vv = sij(w,r,s)

  ws = Vector{Int}[]
  ll = len(w)

  for i in 1:r-1
    ww = sij(vv, i, r )
    if len(ww)==ll
      push!(ws,ww)
    end
  end

  return [(r,s),vv, ws]

end


function max_transition( w::Vector{Int} )
# do step in Lascoux-Schutzenberger transition tree

  pos = findlast2143(w)

  if length(pos)==0
    return []
  end

  r,s = pos[3],pos[4]

  vv = sij(w,r,s)

  ws = Vector{Int}[]
  ll = len(w)

  for i in 1:r-1
    ww = sij(vv, i, r )
    if len(ww)==ll
      push!(ws,ww)
    end
  end

  return [(r,s),vv, ws]

end


function trimw( w::Vector{Int} )
# remove trailing w(i)=i
  n=length(w)

  ww = deepcopy(w)
  while n>0 && ww[n]==n
    pop!(ww)
    n=length(ww)
  end

  return ww

end
