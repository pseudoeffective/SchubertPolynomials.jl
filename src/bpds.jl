# Tools for generating BPDs in Julia
# David Anderson, 29 February 2024.



export BPD, Rothe, all_bpds, all_Kbpds, flat_bpds, bpd2asm, asm2bpd

#############
# BPD type and constructors


"""
    BPD

A type for bumpless pipe dreams

## Structure

A BPD `b` has one field, `b.m`, which is a Matrix{Int8}.  The entries are integers 0-5 which encode the six possible tiles as follows:

   `O` <-> `0`

   `+` <-> `1`

   `/` <-> `2` (r-elbow)

   `%` <-> `3` (j-elbow)

   `|` <-> `4`

   `-` <-> `5`

## Constructor

The function `BPD(m)` takes as its argument `m` either a matrix with entries of type Int8 (with values 0-5) or of type String (with the six possible tiles).


## Example
```julia-repl
julia> mtx = Matrix( [ 0 0 2 ; 0 2 1 ; 2 1 1 ] );

julia> b = BPD( mtx )
O O / 
O / + 
/ + + 

julia> b == mtx
false

julia> b.m == mtx
true

julia> mtx2 = Matrix( [ "O" "O" "/"; "O" "/" "+"; "/" "+" "+" ] )
O O / 
O / + 
/ + + 

julia> b2 = BPD(mtx2);

julia> b2 == b
true

julia> b2.m == b.m
true
```
"""
struct BPD
    m::Matrix{Int8}
end

# Symbol to integer mapping
const SIXVTX_TO_INT = Dict(
    "O" => 0,
    "+" => 1,
    "/" => 2,
    "%" => 3,
    "|" => 4,
    "-" => 5
)

function BPD(matrix::Matrix{String})
    int_matrix = map(x -> SIXVTX_TO_INT[x], matrix)
    return BPD(int_matrix)
end


# convert integers back to symbols for display
# integers 0-5 map to bpd symbols, 6-8 for drifts, 9?
function int_to_symbol(i::Int8)
    symbols = ['O', '+', '/', '%', '|', '-', '.', '*', ' ', 'o']
    return symbols[i+1]
end

function int_to_symbol(t::Tuple)
    return t[1]
end


# add method to Base.show for BPD display
function Base.show(io::IO, bpd::BPD)
    println(io)
    for i in 1:size(bpd.m, 1)
        for j in 1:size(bpd.m, 2)
            print(io, int_to_symbol(bpd.m[i, j]), " ")
        end
        println(io)
    end
end

# overload identity for BPD type
Base.:(==)(b1::BPD, b2::BPD) = b1.m == b2.m



"""
    Rothe(w)

Construct the Rothe BPD for a permutation

## Argument
`w::Vector{Int}`: A permutation

## Returns
`BPD`: The Rothe bumpless pipe dream for `w`.

## Example

```julia-repl
# Produce the BPD
julia> w = [3,2,5,1,4];

julia> b = Rothe(w)

O O / - - 
O / + - - 
O | | O / 
/ + + - + 
| | | / + 


# View the integer matrix which is stored
julia> b.m
5×5 Matrix{Int8}:
 0  0  2  5  5
 0  2  1  5  5
 0  4  4  0  2
 2  1  1  5  1
 4  4  4  2  1

```
"""
function Rothe(w)
        local n=length(w)
        local m=minimum(w)
        local r=Matrix{String}(undef,n,n)

             for j = 1:n
                if j+m-1<w[1]
                   r[1,j]="O"
                elseif j+m-1==w[1]
                   r[1,j]="/"
                else
                   r[1,j]="-"
                end
             end

             for i=2:n
               for j=1:n
                 if j+m-1<w[i]
                   if r[i-1,j]=="/" || r[i-1,j]=="|" || r[i-1,j]=="+"
                      r[i,j]="|"
                   else
                      r[i,j]="O"
                   end
                 elseif j+m-1==w[i]
                   r[i,j]="/"
                 else
                   if r[i-1,j]=="|" || r[i-1,j]=="/" || r[i-1,j]=="+"
                     r[i,j]="+"
                   else
                     r[i,j]="-"
                   end
                 end
               end
             end

        return(BPD(r))
end



############
# droop/drip/drop moves

function can_droop(bpd,i1,j1,i2,j2)

 # check rectangle bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # check NW and SE corners
    if (bpd.m[i1,j1],bpd.m[i2,j2])  != (2,0)
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2
       if bpd.m[i1,j] in [2,3] || bpd.m[i2,j] in [2,3]
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2
       if bpd.m[i,j1] in [2,3] || bpd.m[i,j2] in [2,3]
         return(false)
       end
    end

 # check inside of rectangle
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.m[i,j] in [2,3]
          return(false)
        end
      end
    end
  return(true)
end


# to generate bpds by drooping from flats, need to allow larger droops.
# can have 3 on SW or NE corners of rectangle.

function can_flat_drop(bpd,i1,j1,i2,j2)

 # check bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # small droop not allowed
    if (i2,j2)==(i1+1,j1+1)
       return(false)
    end

 # check corners
    if (bpd.m[i1,j1],bpd.m[i2,j2]) != (2,0)
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2-1
       if bpd.m[i1,j] in [2,3,0] || bpd.m[i2,j] in [2,3,0]
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2-1
       if bpd.m[i,j1] in [2,3,0] || bpd.m[i,j2] in [2,3,0]
         return(false)
       end
    end

 # check interior
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.m[i,j] in [2,3,0]
          return(false)
        end
      end
    end
  return(true)

end


# drip is the small droop
function can_drip(bpd,i1,j1)

 # check corners
    if (bpd.m[i1,j1], bpd.m[i1+1,j1+1]) != (2,0)
      return(false)
    end

    if !(bpd.m[i1,j1+1] in [3,5])
      return(false)
    end

    if !(bpd.m[i1+1,j1] in [3,4])
      return(false)
    end

  return(true)

end

# K-theory version
function can_Kdroop(bpd,i1,j1,i2,j2)

 # check rectangle bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # check NW and SE corners
    if (bpd.m[i1,j1],bpd.m[i2,j2]) != (2,3)
      return(false)
    end

 # check NE and SW corners
    if !( (bpd.m[i1,j2],bpd.m[i2,j1]) in [(1,4), (5,1)] )
      return false
    end

 # check N and W borders
    for j=j1+1:j2-1
       if bpd.m[i1,j]!=1 && bpd.m[i1,j]!=5
         return(false)
       end
    end
    for i=i1+1:i2-1
       if bpd.m[i,j1]!=1 && bpd.m[i,j1]!=4
         return(false)
       end
    end

 # check S and E borders
    aa=0
    for j=j1+1:j2-1
       if bpd.m[i2,j]==3
         return(false)
       end
       if bpd.m[i2,j]==2
         aa+=1
         if aa>1 return(false) end
       end
    end

    for i=i1+1:i2-1
       if bpd.m[i,j2]==3
         return(false)
       end
       if bpd.m[i,j2]==2
         aa+=1
         if aa>1 return(false) end
       end
    end

 # check inside of rectangle
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.m[i,j] in [2,3]
          return(false)
        end
      end
    end
  return(true)
end





function droop(bpd,i1,j1,i2,j2)
 # assumes can_[flat_]droop==true
    
    local bpd2=deepcopy(bpd.m)
    # set corners of rectangle
    bpd2[i1,j1]=0
    bpd2[i2,j2]=3
    if bpd2[i1,j2]==5
      bpd2[i1,j2]=2
    elseif bpd2[i1,j2]==3
      bpd2[i1,j2]=4
    end
    if bpd2[i2,j1]==4
      bpd2[i2,j1]=2
    elseif bpd2[i2,j1]==3
      bpd2[i2,j1]=5
    end

    # set west edge
    for i=i1+1:i2-1
       if bpd2[i,j1]==4
          bpd2[i,j1]=0
       elseif bpd2[i,j1]==1
          bpd2[i,j1]=5
       end
    end

    # set north edge
    for j=j1+1:j2-1
       if bpd2[i1,j]==5
          bpd2[i1,j]=0
       elseif bpd2[i1,j]==1
          bpd2[i1,j]=4
       end
    end

    # set east edge
    for i=i1+1:i2-1
       if bpd2[i,j2]==0
          bpd2[i,j2]=4
       elseif bpd2[i,j2]==5
          bpd2[i,j2]=1
       end
    end

    # set south edge
    for j=j1+1:j2-1
       if bpd2[i2,j]==0
          bpd2[i2,j]=5
       elseif bpd2[i2,j]==4
          bpd2[i2,j]=1
       end
    end

    return(BPD(bpd2))
end



function Kdroop(bpd,i1,j1,i2,j2)
 # assumes can_[flat_]Kdroop==true
    
    local bpd2=deepcopy(bpd.m)
    # set corners of rectangle
    bpd2[i1,j1]=0
    bpd2[i2,j2]=3
    if bpd2[i1,j2]==5
      bpd2[i1,j2]=2
    elseif bpd2[i1,j2]==3
      bpd2[i1,j2]=4
    end
    if bpd2[i2,j1]==4
      bpd2[i2,j1]=2
    elseif bpd2[i2,j1]==3
      bpd2[i2,j1]=5
    end

    # set east edge
    ii=i2
    for i=i1+1:i2-1
       if bpd2[i,j2]==0
          bpd2[i,j2]=4
       elseif bpd2[i,j2]==5
          bpd2[i,j2]=1
       elseif bpd2[i,j2]==2
          bpd2[i,j2]=1
          ii=i
       end
    end

    # set south edge
    jj=j2
    for j=j1+1:j2-1
       if bpd2[i2,j]==0
          bpd2[i2,j]=5
       elseif bpd2[i2,j]==4
          bpd2[i2,j]=1
       elseif bpd2[i2,j]==2
          bpd2[i2,j]=1
          jj=j
       end
    end

    # set west edge
    for i=i1+1:min(ii,i2)-1
       if bpd2[i,j1]==4
          bpd2[i,j1]=0
       elseif bpd2[i,j1]==1
          bpd2[i,j1]=5
       end
    end
    # straighten bottom pipe
    if ii<i2
      bpd2[ii,j1]=2
      for j=j1+1:j2-1
        if bpd2[ii,j]==0
          bpd2[ii,j]=5
        elseif bpd2[ii,j]==4
          bpd2[ii,j]=1
        end
      end
    end

    # set north edge
    for j=j1+1:min(jj,j2)-1
       if bpd2[i1,j]==5
          bpd2[i1,j]=0
       elseif bpd2[i1,j]==1
          bpd2[i1,j]=4
       end
    end
    # straighten bottom pipe
    if jj<j2
      bpd2[i1,jj]=2
      for i=i1+1:i2-1
        if bpd2[i,jj]==0
          bpd2[i,jj]=4
        elseif bpd2[i,jj]==5
          bpd2[i,jj]=1
        end
      end
    end

    return(BPD(bpd2))
end



# produce all droops of bpd
function all_droops(bpd)
   local n=size(bpd.m)[1]

   local dps = []

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_droop(bpd,i1,j1,i2,j2)
              local bpd2=droop(bpd,i1,j1,i2,j2)
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end


# produce all (flat) drops of bpd
function flat_drops(bpd)
   local n=size(bpd.m)[1]

   local dps = []

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_flat_drop(bpd,i1,j1,i2,j2)
              bpd2=makeflat(droop(bpd,i1,j1,i2,j2))
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end


# produce all drips of bpd
function all_drips(bpd)
   local n=size(bpd.m)[1]

   local dps = []

   for i1=1:n-1
     for j1=1:n-1
            if can_drip(bpd,i1,j1)
              local bpd2=droop(bpd,i1,j1,i1+1,j1+1)
              push!(dps,bpd2)
       end
     end
   end

   return(dps)
end


# produce all Kdroops of bpd
function all_Kdroops(bpd)
   local n=size(bpd.m)[1]

   local dps = []

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_droop(bpd,i1,j1,i2,j2)
              local bpd2=droop(bpd,i1,j1,i2,j2)
              push!(dps,bpd2.m)
            end
            if can_Kdroop(bpd,i1,j1,i2,j2)
              local bpd2=Kdroop(bpd,i1,j1,i2,j2)
              if !(bpd2.m in dps)
                push!(dps,bpd2.m)
              end
            end
          end
       end
     end
   end

   return(map(BPD,dps))
end





###############
# iterator generating all BPDs for w


struct BPDIterator
    stack::Vector{Tuple{BPD,Vector{BPD}}}
    seen::Set{Matrix}
end


function BPDIterator(bpd::BPD)
    # initialize with the first element
    seen = Set([bpd.m])
    droops = all_droops(bpd)
    stack = [(bpd, droops)]
    return BPDIterator(stack,seen)
end

Base.eltype(::Type{BPDIterator}) = BPD

Base.IteratorSize(::Type{<:BPDIterator}) = Base.SizeUnknown()

function Base.iterate(iter::BPDIterator, state=nothing)

    while !isempty(iter.stack)
        current, droops = pop!(iter.stack)

        unseen_droops = filter( b -> !(b.m in iter.seen), droops )

        for b in unseen_droops
          push!(iter.seen, b.m)  # mark new droop as seen
          push!( iter.stack, (b, all_droops(b)) )
        end

        return( current, isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # end of iteration
end

########

"""
    all_bpds(w)

An iterator generating all reduced BPDs for a permutation `w`

## Argument
`w::Vector{Int}`: a permutation

## Returns
`BPDIterator`: an iterator type generating all reduced BPDs for `w`.

## Examples

```julia-repl
# Define the iterator
julia> w = [3,2,5,1,4];

julia> bps = all_bpds(w);

# Run a loop over the iterator

julia> i=0; for b in bps i+=1; end; i
3

# The iterator is exhausted; to reset it, define it again

julia> bps = all_bpds(w);

# Form a vector of all BPDs for w

julia> bpds = collect(bps)
3-element Vector{Any}:
 
O O / - - 
O / + - - 
O | | O / 
/ + + - + 
| | | / + 

 
O O / - - 
O O | / - 
O / + % / 
/ + + - + 
| | | / + 

 
O O O / - 
O / - + - 
O | / % / 
/ + + - + 
| | | / + 
```
"""
function all_bpds(w)
    local bpd = Rothe(w)

    iter = BPDIterator(bpd)

    return iter
end


##########
# iterator generating all K-BPDs for w

struct KBPDIterator
    stack::Vector{Any}
    seen::Set{Matrix}
end


function KBPDIterator(bpd::BPD)
    # Initialize with the first element
    seen = Set([bpd.m])
    droops = all_Kdroops(bpd)
    stack = [(bpd, droops)]
    return KBPDIterator(stack,seen)
end

Base.eltype(::Type{KBPDIterator}) = BPD

Base.IteratorSize(::Type{<:KBPDIterator}) = Base.SizeUnknown()


function Base.iterate(iter::KBPDIterator, state=nothing)

    while !isempty(iter.stack)
        current, droops = pop!(iter.stack)

        unseen_droops = filter( b -> !(b.m in iter.seen), droops )

        for b in unseen_droops
          push!(iter.seen, b.m)  # mark new droop as seen
          push!( iter.stack, (b, all_Kdroops(b)) )
        end

        return( current, isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # End of iteration
end



"""
    all_Kbpds(w)

An iterator generating all BPDs for a permutation `w`, including non-reduced K-theoretic ones

## Argument
`w::Vector{Int}`: a permutation

## Returns
`KBPDIterator`: an iterator type generating all BPDs for `w`.

## Examples
```julia-repl
# Define the iterator
julia> w = [3,2,5,1,4]

julia> bps = all_Kbpds(w);

# Run a loop over the iterator

julia> i=0; for b in bps i+=1; end; i
4

# To reset the iterator, define it again

julia> bps = all_Kbpds(w);

# Form a vector of all BPDs for w

julia> bpds = collect(bps)
4-element Vector{Any}:
 
O O / - - 
O / + - - 
O | | O / 
/ + + - + 
| | | / + 

 
O O / - - 
O O | / - 
O / + % / 
/ + + - + 
| | | / + 

 
O O O / - 
O O / + - 
O / + % / 
/ + + - + 
| | | / + 

 
O O O / - 
O / - + - 
O | / % / 
/ + + - + 
| | | / + 


```
"""
function all_Kbpds(w)
    local bpd = Rothe(w)
    iter = KBPDIterator(bpd)

    return iter
end

#######


# determine if bpd is flat
function isflat(bpd)
   local n=size(bpd.m)[1]

   for i=2:n-1
     for j=2:n-1
       if bpd.m[i,j]==0 && bpd.m[i-1,j]!=0 && bpd.m[i,j-1]!=0 && bpd.m[i-1,j-1]==2
          return(false)
       end
     end
   end
   return(true)
end


# returns flat bpd in drift class of bpd
function makeflat(bpd)
   local n=size(bpd.m)[1]

   for i=2:n-1
     for j=2:n-1
       if bpd.m[i,j]==0 && bpd.m[i-1,j]!=0 && bpd.m[i,j-1]!=0 && bpd.m[i-1,j-1]==2
         local bpd2=droop(bpd,i-1,j-1,i,j)
         return makeflat(bpd2)
       end
     end
   end

   return(bpd)   

end   
 


#############
# iterator generating all flat BPDs for w

struct FlatBPDIterator
    stack::Vector{Any}
    seen::Set{Matrix}
end

function FlatBPDIterator(bpd::BPD)
    # initialize with the first element
    seen = Set([makeflat(bpd).m])
    drops = flat_drops(makeflat(bpd))
    stack = [(makeflat(bpd), drops)]
    return FlatBPDIterator(stack,seen)
end

Base.eltype(::Type{FlatBPDIterator}) = BPD

Base.IteratorSize(::Type{<:FlatBPDIterator}) = Base.SizeUnknown()


function Base.iterate(iter::FlatBPDIterator, state=nothing)

    while !isempty(iter.stack)
        current, drops = pop!(iter.stack)

        unseen_drops = filter( b -> !(makeflat(b).m in iter.seen), drops )

        for b in unseen_drops
          b=makeflat(b)
          push!(iter.seen, b.m)  # mark new drop as seen
          push!( iter.stack, (b, flat_drops(b)) )
        end

        return( makeflat(current), isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # end of iteration
end


"""
    flat_bpds(w)

An iterator generating all flat reduced BPDs for a permutation `w`

## Argument
`w::Vector{Int}`: a permutation

## Returns
`FlatBPDIterator`: an iterator type generating all flat reduced BPDs for `w`.

## Examples
```julia-repl
# Define the iterator
julia> w = [3,2,5,1,4];

julia> fbps = flat_bpds(w);

# Run a loop over the iterator

julia> i=0; for b in fbps i+=1 end; i
3

# To reset the iterator, define it again

julia> fbps = flat_bpds(w);

# Form a vector of flat BPDs for w

julia> fbpds = collect(fbps)
3-element Vector{Any}:
 
O O / - - 
O / + - - 
O | | O / 
/ + + - + 
| | | / + 

 
O O / - - 
O O | / - 
O / + % / 
/ + + - + 
| | | / + 

 
O O O / - 
O / - + - 
O | / % / 
/ + + - + 
| | | / + 
```
"""
function flat_bpds(w)
    local bpd = Rothe(w)
    iter = FlatBPDIterator(bpd)

    return iter
end


function vec_flat_bpds(w)
    local bpd = Rothe(w)
    iter = FlatBPDIterator(bpd)

    return collect(iter)
end



#############
# Conversions to and from ASMs
# where to put these functions?


"""
    bpd2asm(b)

Convert a bumpless pipe dream to an alternating sign matrix

## Argument
`b::BPD`: a bumpless pipe dream

## Returns
`Matrix`: the alternating sign matrix corresponding to `b`.

## Example
```julia-repl
# Generate all reduced BPDs for a permutation
julia> w = [3,2,5,1,4];

julia> bps = all_Kbpds(w);

# Construct a vector of the corresponding ASMs

julia> asms = [];

julia> for b in bps push!(asms, bpd2asm(b)) end;

julia> asms[3]
5×5 Matrix{Int8}:
 0  0  0   1  0
 0  0  1   0  0
 0  1  0  -1  1
 1  0  0   0  0
 0  0  0   1  0

```
"""
function bpd2asm( b::BPD )

  local n=size(b.m)[1]

  local a=zeros(Int8,n,n)

  for i=1:n
    for j=1:n
      if b.m[i,j]==2
        a[i,j]=1
      elseif b.m[i,j]==3
        a[i,j]=-1
      end
    end
  end

  return a
end


"""
    asm2bpd(a)

Convert an alternating sign matrix to a bumpless pipe dream

## Argument
`a::Matrix`: an alternating sign matrix

## Returns
`BPD`: the corresponding BPD

## Example
```julia-repl
# Generate all reduced BPDs for a permutation
julia> w = [3,2,5,1,4];

julia> bpds = collect(all_Kbpds(w));

# Convert a BPD to an ASM

julia> b = bpds[3]

O O O / - 
O O / + - 
O / + % / 
/ + + - + 
| | | / + 

julia> a = bpd2asm(b)
5×5 Matrix{Int8}:
 0  0  0   1  0
 0  0  1   0  0
 0  1  0  -1  1
 1  0  0   0  0
 0  0  0   1  0

# Convert the ASM back to a BPD

julia> b2 = asm2bpd(a);

julia> b == b2
true

```
"""
function asm2bpd( a )
# improve this, the rules are not local

  local n=size(a)[1]

  local b=Matrix{Int8}(undef,n,n)

  if a[n,1]==1
    b[n,1]=2
  else
    b[n,1]=4
  end

  for j=2:n
    if a[n,j]==1
      b[n,j]=2
    elseif b[n,j-1]==2 || b[n,j-1]==1
      b[n,j]=1
    else
      b[n,j]=4
    end
  end

  for i=n-1:-1:1
    for j=1:n
      if a[i,j]==1
        b[i,j]=2

      elseif a[i,j]==-1
        b[i,j]=3

      elseif a[i,j]==0

        if j==1
          if b[i+1,j]==4
            b[i,j]=4
          else
            b[i,j]=0
          end

        else
          if b[i,j-1]==2 || b[i,j-1]==1 || b[i,j-1]==5
            local bb=5
          else
            bb=0
          end
          if b[i+1,j]==3 || b[i+1,j]==1 || b[i+1,j]==4
            local cc=4
          else
            cc=0
          end
          if bb==5
            if cc==4
              b[i,j]=1
            else
              b[i,j]=5
            end
          else
            b[i,j]=cc
          end
        end
      end
    end
  end

  return BPD(b)
end



####################
#not used
function can_sharp_drop(bpd,i1,j1,i2,j2)

 # check bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # small droop not allowed
    if (i2,j2)==(i1+1,j1+1)
       return(false)
    end

 # check corners
    if bpd.m[i1,j1] != 2 || bpd.m[i2,j2] != 0
      return(false)
    end

 # check active pipe
    if i2==i1+1 && bpd.m[i2-1,j2-1] == 5
      return(false)
    end
    if j2==j1+1 && bpd.m[i2-1,j2-1]==4
      return(false)
    end

 # check NW of destination
    if bpd.m[i2-1,j2-1]==0
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2-1
       if bpd.m[i1,j]==2 || bpd.m[i2,j]==2 || bpd.m[i1,j]==3 || bpd.m[i2,j]==3
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2-1
       if bpd.m[i,j1]==2 || bpd.m[i,j2]==2 || bpd.m[i,j1]==3 || bpd.m[i,j2]==3
         return(false)
       end
    end

 # check interior
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.m[i,j] == 2 || bpd.m[i,j]==3
          return(false)
        end
      end
    end
  return(true)
end



function sharp_drops(bpd)
# produce all (sharp) drops of bpd
   local n=size(bpd.m)[1]

   local dps = []

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_sharp_drop(bpd,i1,j1,i2,j2)
              bpd2=droop(bpd,i1,j1,i2,j2)
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end



function issharp(bpd)
# determine if bpd is sharp
   local n=size(bpd.m)[1]

   for i=1:n-1
     for j=1:n-1
       if bpd.m[i,j]==0 && bpd.m[i+1,j]!=0 && bpd.m[i,j+1]!=0 && bpd.m[i+1,j+1]!=1
          return(false)
       end
     end
   end
   return(true)
end



function makesharp(bpd)
# returns sharp bpd in drift class of bpd
   local n=size(bpd.m)[1]

   for i=1:n-2
     for j=1:n-2
       if bpd.m[i,j]==0 && bpd.m[i+1,j]!=0 && bpd.m[i,j+1]!=0 && bpd.m[i+1,j+1]==3
         local bpd2=deepcopy(bpd.m)
         bpd2[i,j]=2
         bpd2[i+1,j+1]=0

         if bpd2[i+1,j]==2
           bpd2[i+1,j]=4
         else
           bpd2[i+1,j]=3
         end

         if bpd2[i,j+1]==2
           bpd2[i,j+1]=5
         else
           bpd2[i,j+1]=3
         end

         return makesharp(BPD(bpd2))
       end
     end
   end

   return(bpd)   

end

