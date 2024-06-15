# Tools for working with drift configurations in Julia
# David Anderson, February 2024.


export Drift, drift_class, drift_poly, nw_reset, se_reset, markconfig, random_drift, partition2drift, bpd2drift


# TO DO: review function drift2bin(d::Drift, R::DoublePolyRing) which acts like bpd2bin, records a product of binomials x[i]+y[j]
# TO DO: review function drift_poly(d::Drift, R::DoublePolyRing) which produces the drift polynomial from the iterator
# TO DO: clarify the logic in markbox

#########

struct Drift
    m::Matrix{Union{Int8,Tuple}}
end



# Symbol to integer mapping
const DRIFT_TO_INT = Dict(
    "O" => Int8(0),
    "+" => Int8(1),
    "." => Int8(6),
    "*" => Int8(7),
    "" => Int8(8)
)

function Drift(matrix::Matrix{String})
    int_matrix = map(x -> DRIFT_TO_INT[x], matrix)
    return Drift(int_matrix)
end

#= moved to bpds.jl

# convert integers back to symbols for display
function int_to_symbol(i::Int8)
    symbols = ['O', '+', '/', '%', '|', '-', '.', '*', ' ', 'o']
    return symbols[i+1]  
end

function int_to_symbol(t::Tuple)
    return t[1]
end

=#


# add method to display Drift
function Base.show(io::IO, dc::Drift)
    println(io)
    for i in 1:size(dc.m, 1)
        for j in 1:size(dc.m, 2)
            print(io, int_to_symbol(dc.m[i, j]), " ")
        end
        println(io)
    end
end

# overload identity for Drift type
Base.:(==)(d1::Drift, d2::Drift) = d1.m == d2.m



function can_drift(dc,i1,j1)

 # check corners
    if dc.m[i1,j1] != 0 && !isa( dc.m[i1,j1], Tuple )
      return(false)
    end

    if dc.m[i1+1,j1] == 0 || isa( dc.m[i1+1,j1], Tuple ) || dc.m[i1+1,j1] == 1 || dc.m[i1+1,j1]==7
      return(false)
    end

    if dc.m[i1,j1+1] == 0 || isa( dc.m[i1,j1+1], Tuple ) || dc.m[i1,j1+1] == 1 || dc.m[i1,j1+1]==7
      return(false)
    end

    if dc.m[i1+1,j1+1] == 0 || isa( dc.m[i1+1,j1+1], Tuple ) || dc.m[i1+1,j1+1] == 1 || dc.m[i1+1,j1+1]==7
      return(false)
    end


  return(true)

end

function drift( dc, i, j)
# perform drift move of dc at i,j if possible

  if can_drift(dc,i,j)

    dc2=copy(dc.m)

    if dc2[i+1,j+1]==6
      dc2[i+1,j+1]=Int8(7)
    elseif isa( dc2[i,j], Tuple )
      dc2[i+1,j+1]=(dc2[i,j][1]-1, dc2[i,j][2])
    else
      dc2[i+1,j+1]=Int8(0)
    end

    dc2[i,j]=Int8(8)

    return Drift(dc2)

  end

end 


function can_undrift(dc,i1,j1)

 # check corners
    if dc.m[i1,j1] != 0 && dc.m[i1,j1] != 7
      return(false)
    end

    if dc.m[i1-1,j1] == 0 || dc.m[i1-1,j1] == 1 || dc.m[i1-1,j1]==7
      return(false)
    end

    if dc.m[i1,j1-1] == 0 || dc.m[i1,j1-1] == 1 || dc.m[i1,j1-1]==7
      return(false)
    end

    if dc.m[i1-1,j1-1] == 0 || dc.m[i1-1,j1-1] == 1 || dc.m[i1-1,j1-1]==7 || dc.m[i1-1,j1-1]==6
      return(false)
    end

  return(true)

end

function undrift( dc, i, j)
# perform drift move of dc at i,j if possible

  if can_undrift(dc,i,j)

    dc2=copy(dc.m)

    if dc2[i,j]==7
      dc2[i,j]=Int8(6)
    else
      dc2[i,j]=Int8(8)
    end

    dc2[i-1,j-1]=Int8(0)

    return Drift(dc2)

  end

end 


function step_drifts(dc)
# produce all one-step drifts of dc
   local n=size(dc.m)[1]

   local dfts = []

   for i=1:n-1
     for j=1:n-1
            if can_drift(dc,i,j)
              local dc2=drift(dc,i,j)
              push!(dfts,dc2)
       end
     end
   end

   return(dfts)
end



struct DriftIterator
  stack::Vector{Any}
  seen::Set{Matrix}
end


function DriftIterator(dc::Drift)
  # Initialize
  seen = Set([dc.m])
  dfts = step_drifts(dc)
  stack = [(dc, dfts)]
  return DriftIterator(stack,seen)
end


Base.eltype(::Type{DriftIterator}) = Drift


Base.IteratorSize(::Type{<:DriftIterator}) = Base.SizeUnknown()


function Base.iterate(iter::DriftIterator, state=nothing)

    while !isempty(iter.stack)
        current, drfts = pop!(iter.stack)

        unseen_drfts = filter( b -> !(b.m in iter.seen), drfts )

        for b in unseen_drfts
          push!(iter.seen, b.m)  # mark new drift as seen
          push!( iter.stack, (b, step_drifts(b)) )
        end

        return( current, isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # End of iteration
end


function drift_class(dc::Drift)
# iterator of all diagrams in drift class of dc

  dc2 = nw_reset(dc)

  iter = DriftIterator(dc2)

#  empty!(hash_all_drifts_SE)  # clear the lookup table

  return iter

end



function nw_reset(dc)
# returns flat diagram in drift class of dc
   local n=size(dc.m)[1]

   for i=2:n
     for j=2:n
       if can_undrift(dc,i,j)
         local dc2=undrift(dc,i,j)
         return nw_reset(dc2)
       end
     end
   end

   return(dc)   

end   


function se_reset(dc)
# returns sharp diagram in drift class of dc
   local n=size(dc.m)[1]

   for i=1:n-1
     for j=1:n-1

       if can_drift(dc,i,j)
         dc2=copy(dc)
         dc2=drift(dc2,i,j)
         return se_reset(dc2)
       end
     end
   end

   return(dc)   

end



#####
# Drift polynomials
#####

function drift2bin( d::Drift, R::DoublePolyRing=xy_ring( size(d.m)[1]-1, size(d.m)[2]-1 )[1] )
# product of binomials for d
# requires DoublePolyRing
# can get single polyn by using no y_vars
  local n=size(d.m)[1]-1
  bin = R.ring(1)

  x = R.x_vars
  y = R.y_vars

  local aa=length(x)
  local bb=length(y)

  for i=1:n
    for j=1:n

      if d.m[i,j]==0 || d.m[i,j]==7 || isa(d.m[i,j],Tuple)
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


function drift_poly( d::Drift, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1]  )
# compute drift pol by iterator
  dc=drift_class(d)

  pol=R.ring(0)

  for dd in dc
    pol = pol+drift2bin(dd,R)
  end

  return(pol)

end


#####
# Towards drift transition
#####


 
function drift_split( dc, i, j )
# do split of dc at (i,j)

# assumes input dc is marked

dr=dc.m

# only split at colliding box
  if !isa(dr[i,j], Tuple)
    return dc
  end

  if !dr[i,j][2]
    return dc
  end

# find colliding box at end of column or row
  if isa(dr[i+1,j], Tuple) && dr[i+1,j][2]
    return drift_split(dc,i+1,j)
  end

  if isa(dr[i,j+1], Tuple) && dr[i,j+1][2]
    return drift_split(dc,i,j+1)
  end


  local dr1=copy(dr)

  local k=dr[i,j][1]

  if k>0
    dr1[i+k,j+k]=Int8(6)
    dr1[i,j]=(k,false)
  else
    dr1[i,j]=(0,false,0)
  end

  dc1=Drift(dr1)

  local dr2=copy(dr)

  dr2[i+k+1,j+k+1]=Int8(0)
  dr2[i,j]=Int8(8)

# changed 6 to 9
  dr2[i+k,j+k]=Int8(9)

  b=1
  while isa( dr2[i,j+b], Tuple )
    if dr2[i+k+1,j+b+k+1]==6
      dr2[i+k+1,j+b+k+1]=(0,false,0)
    else
      dr2[i+k+1,j+b+k+1]=( dr[i,j+b][1]-k-1, false )
    end
    dr2[i,j+b]=Int8(8)
    b+=1
  end

  b=1
  while isa( dr2[i+b,j], Tuple )
    if dr2[i+b+k+1,j+k+1]==6
      dr2[i+b+k+1,j+k+1]=(0,false,0)
    else
      dr2[i+b+k+1,j+k+1]=( dr[i+b,j][1]-k-1, false )
    end
    dr2[i+b,j]=Int8(8)
    b+=1
  end

  dc2=Drift(dr2)
  dc2=unmarkconfig(dc2)
  dc2=markconfig(dc2)

  return (dc1,dc2)

end


# better: label boxes with pairs (int,boolean), int=distance box can drift, boolean=collides or not
# probably should adapt the marking functions so they can work on marked diagrams, not just the unmarked drift diagrams
# no, this leads to errors -- for now an inefficient "unmarking"
# done?  check!


function markbox( dc::Drift, i, j )

  dr=dc.m

  local n = size(dr)[1]

  if dr[i,j]==7
    return(0,false,0)
  end

  if dr[i,j]!=0 && !isa(dr[i,j],Tuple)
    return dr[i,j]
  end

  if i==n || j==n
    return (0,false)
  end

# this case should never occur for bpds
  if (dr[i+1,j]==8 || dr[i+1,j]==6) && (dr[i,j+1]==8 || dr[i,j+1]==6) && dr[i+1,j+1]==0
      return (0,true)
  end

  if (dr[i,j+1]==0 || isa( dr[i,j+1], Tuple ) ) && ( dr[i+1,j]!=0 && !isa( dr[i+1,j], Tuple ) )
    local k1=markbox( dc, i, j+1 )[1]
    local k2=0
#    while i+k2+1<n && j+k2<n && dr[i+k2+1,j+k2+1]==8
    while k2<k1 && i+k2+1<n && j+k2<n
      if dr[i+k2+2,j+k2+1]==0 || isa( dr[i+k2+2,j+k2+1], Tuple )
        return (k2,true)
      end
      k2+=1
      if dr[i+k2,j+k2]==6
        return (k2,false)
      end
    end
    return (k2,false)
  end

  if ( dr[i+1,j]==0 || isa( dr[i+1,j], Tuple ) ) && dr[i,j+1]!=0 && !isa( dr[i,j+1], Tuple )
    local k1=markbox( dc, i+1, j )[1]
    local k2=0
#    while i+k2<n && j+k2+1<n && dr[i+k2+1,j+k2+1]==8
    while k2<k1 && i+k2<n && j+k2+1<n
      if dr[i+k2+1,j+k2+2]==0 || isa( dr[i+k2+1,j+k2+2], Tuple )
        return (k2,true)
      end
      k2+=1
      if dr[i+k2,j+k2]==6
        return (k2,false)
      end
    end
    return (k2,false)
  end

  if ( dr[i+1,j]==0 || isa( dr[i+1,j], Tuple ) ) && ( dr[i,j+1]==0 || isa( dr[i,j+1], Tuple ) )
    local k1 = markbox( dc, i+1, j )[1]
    local k2 = markbox( dc, i, j+1 )[1]
    return ( minimum([k1,k2]), false )
  end


  if can_drift( dc, i, j)

    if dr[i+1,j+1]!=6

      if i<n-1 && dr[i+2,j+1]==0
        return (0,true)
      end

      if j<n-1 && dr[i+1,j+2]==0
        return (0,true)
      end

    else

      return (1,false)

    end

    local dc1 = drift(dc,i,j)

    local (kk,collides)=markbox( dc1, i+1, j+1 )

    return (kk+1,collides)
  end

  return (0,false)

end


function markconfig( dc::Drift )
# mark interfering boxes in drift config dc

  local n=size(dc.m)[1]

  local mm=Matrix{Any}(undef,n,n)

  for i=1:n
    for j=1:n
      mm[i,j]=markbox( dc, i, j )
    end
  end

  return Drift(mm)

end


function markconfig( bpd::BPD )

  return markconfig( bpd2drift(bpd) )

end


function unmarkbox( dc::Drift, i, j)

  if length(dc.m[i,j])>2
    return Int8(7)
  elseif isa( dc.m[i,j], Tuple )
    return Int8(0)
  else
    return dc.m[i,j]
  end

end

function unmarkconfig(dc::Drift)

  local n=size(dc.m)[1]

  local mm=Matrix{Any}(undef,n,n)

  for i=1:n
    for j=1:n
      mm[i,j]=unmarkbox(dc,i,j)
    end
  end

  return Drift(mm)

end


#####
# Generating drift configurations
#####

function random_drift( n )
# random drift config of size n

  local possible_entries = Int8[0, 1, 8, 6, 7]

  return Drift(rand( possible_entries, n,n ))

end

# make drift config from partition, so that boxes can drift to rows bounded by ff
function partition2drift( lambda::Vector{Int}, ff::Vector{Int}=fill(length(lambda),length(lambda)), n::Int=maximum( [ maximum(ff), length(lambda), lambda[1] ] )  )
# drift config with lambda in NW corner

  local mtx = fill( Int8(8), n, n)

  for i=1:length(lambda)
    for j=1:lambda[i]
      mtx[i,j]=Int8(0)
    end
    if i<=length(ff) && ff[i]+1 <= n && lambda[i]+ff[i]+1-i <= n
      mtx[ff[i]+1,lambda[i]+ff[i]-i+1]=Int8(1)
    end
  end

  return Drift(mtx)

end



function bpd2drift( bpd::BPD )
# generate drift config from BPD

  local n = size(bpd.m)[1]

  local bpd2=copy(bpd.m)


  for i=1:n
    for j=1:n
      if !isa(bpd2[i,j],Tuple) && bpd2[i,j]!=0 && bpd2[i,j]!=1
        bpd2[i,j]=Int8(8)
      end
    end
  end


  return Drift(bpd2)

end





#=
# not used

function tabcomps(bpd)
# return labelled tableaux for a flat bpd

  local n=size(bpd)[1]

  if !isflat(bpd)
    return( tableau_components( flatten(bpd) ) )
  end

  local lyds=Vector{Vector}([])

  local corners=Vector{Tuple{Int,Int}}([])

  for i=1:n
    for j=1:n
      if !( (i,j) in corners) && bpd[i,j]==0 && ((i,j)==(1,1) || (i>1 && j>1 &&bpd[i-1,j-1]==1)) #find a new NW corner
        push!(corners,(i,j))

        local la=Vector{Int}([])
        local mu=Vector{Int}([])
        local rr=Vector{Int}([])

        local s=0
        while i+s<=n && bpd[i+s,j]==0

          local k=0
          while j+k<=n && bpd[i+s,j+k]==0  # find SE boxes
            k +=1
          end          
          push!(la,k)

          local el=1
            while i+s+el<=n && j+k-1+el<=n && bpd[i+s+el,j+k-1+el]=="%" || bpd[i+s+el,j+k-1+el]=="|"  # || bpd[i+s+el,j+k-1+el]=="-" 
              el +=1
            end
          push!(rr,el-1)


          local kk=0
          while j-kk-1>0 && bpd[i+s,j-kk-1]==0  # find skew boxes
            kk +=1
          end

          mu=mu+fill(kk,length(mu))
          la=la+fill(kk,length(la))
          push!(mu,0)

          if s>0 && i+s>1 && j-kk>1 && bpd[i+s-1,j-kk-1]==1
            push!(corners,(i+s,j-kk) ) # record new corner
          end
          j=j-kk
          s +=1
        end

#=
        for k=length(la):-1:2
          if la[k-1]==la[k]
            rr[k-1]=rr[k]
          end
        end
=#

        push!(lyds,[la,rr,mu,[i-1,j-1]])

      end
    end
  end

  return lyds
end

=#