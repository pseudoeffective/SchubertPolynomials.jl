# Tools for working with drift configurations in Julia
# David Anderson, May 2025.


export Drift, drift_class, nw_reset, se_reset, random_drift, partition2drift, bpd2drift, isintegrable, iscancelable, drift2rkmtx, rkmtx2asm, empty_drift


# TO DO: improve Base.show overload, better cross symbol

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



# convert integers back to symbols for display
function int_to_symbol_drift(i::Int8)
    symbols = ["\u25A1",  # "□ "
        "\u002B",    # \u271A "✚ "
        "\u256D\u2500",    # "╭─"
        "\u256F",         # "╯ "
        "\u2502",         # "│ "
        "\u2500\u2500",    # "──"
        "\u2B27",         # \u2022 "• "
        '*', "\u00B7", 'o'
    ]
    return symbols[i+1]  
end

function int_to_symbol_drift(t::Tuple)
    return t[1]
end




# add method to display Drift
function Base.show(io::IO, dc::Drift)
    println(io)
    for i in 1:size(dc.m, 1)
        for j in 1:size(dc.m, 2)
            print(io, int_to_symbol_drift(dc.m[i, j]), " ")
        end
        println(io)
    end
end

# overload identity for Drift type
Base.:(==)(d1::Drift, d2::Drift) = d1.m == d2.m

# add method to Base.size for Drift
function Base.size(d::Drift)
  size(d.m)
end


function can_drift(dc::Drift,i::Int,j::Int)

 # check corners
    if dc.m[i,j] != 0 && !isa( dc.m[i,j], Tuple )
      return(false)
    end

    if dc.m[i+1,j] == 0 || isa( dc.m[i+1,j], Tuple ) || dc.m[i+1,j] == 1 || dc.m[i+1,j]==7
      return(false)
    end

    if dc.m[i,j+1] == 0 || isa( dc.m[i,j+1], Tuple ) || dc.m[i,j+1] == 1 || dc.m[i,j+1]==7
      return(false)
    end

    if dc.m[i+1,j+1] == 0 || isa( dc.m[i+1,j+1], Tuple ) || dc.m[i+1,j+1] == 1 || dc.m[i+1,j+1]==7
      return(false)
    end


  return(true)

end

function drift( dc::Drift, i::Int, j::Int)
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


function can_undrift(dc::Drift,i::Int,j::Int)

 # check corners
    if dc.m[i,j] != 0 && dc.m[i,j] != 7
      return(false)
    end

    if dc.m[i-1,j] == 0 || dc.m[i-1,j] == 1 || dc.m[i-1,j]==7
      return(false)
    end

    if dc.m[i,j-1] == 0 || dc.m[i,j-1] == 1 || dc.m[i,j-1]==7
      return(false)
    end

    if dc.m[i-1,j-1] == 0 || dc.m[i-1,j-1] == 1 || dc.m[i-1,j-1]==7 || dc.m[i-1,j-1]==6
      return(false)
    end

  return(true)

end

function undrift( dc::Drift, i::Int, j::Int)
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


function step_drifts(dc::Drift)
# produce all one-step drifts of dc
   local n=size(dc)[1]

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



function nw_reset(dc::Drift)
# returns flat diagram in drift class of dc
   local n,m=size(dc)

   for i=2:n
     for j=2:m
       if can_undrift(dc,i,j)
         local dc2=undrift(dc,i,j)
         return nw_reset(dc2)
       end
     end
   end

   return(dc)   

end   


function se_reset(dc::Drift)
# returns sharp diagram in drift class of dc
   local n,m=size(dc)

   for i=1:n-1
     for j=1:m-1

       if can_drift(dc,i,j)
         dc2=copy(dc)
         dc2=drift(dc2,i,j)
         return se_reset(dc2)
       end
     end
   end

   return(dc)   

end


function can_cancel(dc::Drift, i::Int,j::Int)

  if i>=size(dc)[1] || j>=size(dc)[2] return false end
 
  if dc.m[i,j]!=0 return false end

  if can_drift(dc,i,j) return can_cancel(drift(dc,i,j),i+1,j+1) end
 
  if dc.m[i+1,j+1]!=1 return false end
  if (dc.m[i+1,j] in [0,1]) || (dc.m[i,j+1] in [0,1]) return false end

  return true

end


function cancel_drift!(dc::Drift, i::Int, j::Int)

  if !can_cancel(dc,i,j) return dc end

  dc.m[i,j]=Int8(8)
  for s=1:min(size(dc)[1]i,size(dc)[2]-j)
    if dc.m[i+s,j+s]==1
      dc.m[i+s,j+s]=Int8(8)
      return dc
    end
  end

end

function cancel_all_drift(dc::Drift)

   n,m=size(dc)
   if n!=m throw(error("Drift configuration must be square")) end
  
  if !all(x -> x == 0 || x == 1 || x==8, dc.m)
    return dc
  end

  dc2=deepcopy(dc)

  while countboxes(dc2)>0
    ctbx=countboxes(dc2)
    for s=-n+1:0
      for t=0:s+n-1
        cancel_drift!(dc2,n-t,n+s-t)
      end
    end
    for s=1:n-1
      for t=0:n-s-1
        cancel_drift!(dc2,n-s-t,n-t)
      end
    end
    if ctbx==countboxes(dc2) return dc2 end
  end

  return dc2
end

function iscancelable(dc::Drift)
  n,m=size(dc)
  if n!=m return false end

  return (cancel_all_drift(dc)==empty_drift(n))
end


function drift2rkmtx(dc::Drift)

  n,m=size(dc)
  rkmtx = fill(0,n,m)
  #first row
  for j=1:m
    if !(dc.m[1,j] in [0,1])
      rkmtx[1,j]=1
    elseif dc.m[1,j]==1
      rkmtx[1,j]=2
    end
  end
  #first column
  for i=2:n
    if !(dc.m[i,1] in [0,1])
      rkmtx[i,1]=1
    elseif dc.m[i,1]==1
      rkmtx[i,1]=2
    end
  end
  #rest
  for i=2:n
    for j=2:m
      if dc.m[i,j]==0
        rkmtx[i,j]=rkmtx[i-1,j-1]
      elseif dc.m[i,j]==1
        rkmtx[i,j]=rkmtx[i-1,j-1]+2
      else
        rkmtx[i,j]=rkmtx[i-1,j-1]+1
      end
    end
  end

  return rkmtx
end

function rkmtx2asm(mtx::Matrix{<:Integer})
  n,m=size(mtx)
  aa = fill(0,n,m)
  if n==0 || m==0 return mtx end
  #first row
  aa[1,1]=mtx[1,1]
  for j=2:m
    aa[1,j]=mtx[1,j]-mtx[1,j-1]
  end
  #first column
  for i=2:n
    aa[i,1]=mtx[i,1]-mtx[i-1,1]
  end
  #rest
  for i=2:n
    for j=2:m
      aa[i,j]=mtx[i-1,j-1]+mtx[i,j]-mtx[i-1,j]-mtx[i,j-1]    
    end
  end
  return aa
end

function isintegrable(dc::Drift)
  aa = rkmtx2asm( drift2rkmtx( dc ) )

  return is_asm(aa)
end




#####
# Generating drift configurations
#####

function empty_drift( n::Int, m::Int=n )

  mtx = fill( Int8(8),n,m )

  return Drift(mtx)

end

# random drift config of size n,m
function random_drift( n::Int, m::Int=n ; extended::Bool=false)
  if extended 
    possible_entries = Int8[0, 1, 8, 6, 7]
  else
    possible_entries = Int8[0, 1, 8]
  end

  return Drift(rand( possible_entries, n,m ))

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

  local n = size(bpd)

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

function countboxes( dc::Drift )

  ct = count(==(0),dc.m)
  return ct

end