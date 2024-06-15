# Type for writing sums of Schubert classes
# David Anderson, June 2024.


export SchubertSum

# TO DO: documentation


#########


# type for representing a sum of Schubert classes
mutable struct SchubertSum{T<:Union{Int,ZZMPolyRingElem}}
    coeffs::Vector{T}
    schubs::Vector{Vector{Int}}
end


# WARNING - so far no safety for possibly repeated elements in ss.schubs
# always assume these are distinct


# a single Schubert class
function SchubertSum( w::Vector{Int} )
  ww=trimw(w)
  length(ww)==0 && push!(ww,1)
  return SchubertSum( [1], [ww] )
end


# FIXME meant to ensure length of coeffs and schubs is the same
function SchubertSum(coeffs::Vector{T}, schubs::Vector{Vector{Int}}) where T <: Any
    cfs = Vector{Union{Int,ZZMPolyRingElem}}()
    for c in coeffs
      push!(cfs, c)
    end

    if length(cfs)==length(schubs)
        return SchubertSum(cfs, schubs)

    elseif length(cfs) < length(schubs)
        sch = schubs[1:length(cfs)]
        return SchubertSum(cfs, sch)

    elseif length(cfs) > length(schubs)
        sch = schubs
        for i in length(schubs)+1:length(cfs)
           push!(sch,[1])
        end

        return SchubertSum(cfs, sch)
    end
  
end


function sort_schubert_sum( ss::SchubertSum )

    p = sortperm( ss.schubs )
    sorted_coeffs = ss.coeffs[p]
    sorted_schubs = ss.schubs[p]
    return SchubertSum( sorted_coeffs, sorted_schubs )

end


# overload identity for SchubertSum
import Base: ==

function ==(a::SchubertSum, b::SchubertSum)

  a1 = condense(a)
  b1 = condense(b)

  a1 = sort_schubert_sum(a1)
  b1 = sort_schubert_sum(b1)

  return ( a1.schubs == b1.schubs && a1.coeffs == b1.coeffs )

end


# overload length for SchubertSum
function Base.length(ss::SchubertSum)
    length(ss.coeffs) 
end

# overload coeff for SchubertSum
import Nemo.coeff
function coeff( ss::SchubertSum, w::Vector{Int} )

  wi = findfirst( isequal(w), ss.schubs )

  wi==nothing ? cc=0 : cc=ss.coeffs[wi]

  return cc
end


# overload show for SchubertSum
#import Base.show

function Base.show(io::IO, p::SchubertSum{T}) where T
    if isempty(p.coeffs)
        print(io, "0")
        return
    end
    
    first = true
    for (coeff, schub) in zip(p.coeffs, p.schubs)
        # Handle sign and separation

        if first && string(coeff)[1]=='-'
            print(io, "-")
            coeff=-coeff
        end

        if !first
            if string(coeff)[1]=='-'
                print(io, " - ")
                coeff = -coeff
            else
                print(io, " + ")
            end
        end
        first = false
        
        # Print coefficient, handling multiple terms and suppressing '1' unless it's the only term
        str_coeff = string(coeff)
        if coeff != 1 && coeff != -1 || occursin(" + ", str_coeff) || occursin(" - ", str_coeff[2:end])
            if coeff == -1 && !isempty(p.coeffs)
                print(io, "-")
            elseif occursin(" + ", str_coeff) || occursin(" - ", str_coeff[2:end])  # Check for multiple terms
                print(io, "(", str_coeff, ")*")
            else
                print(io, str_coeff, "*")
            end
        end
        
        # Print the Schubert class, like S[2,3,1]
        print(io, "S", schub)
    end
end


# overload '+' for SchubertSum
import Base: +

function +(ss1::SchubertSum, ss2::SchubertSum) # add two Schubert sums

  length(ss1)>length(ss2) && return +( ss2, ss1 )

  sch = copy( ss2.schubs )
  cfs = copy( ss2.coeffs )

  for i in 1:length(ss1)

    vvi = findfirst(isequal(ss1.schubs[i]),sch)
    if vvi==nothing
      push!(sch,ss1.schubs[i])
      push!(cfs,ss1.coeffs[i])
    else
      cfs[vvi]=cfs[vvi]+ss1.coeffs[i]
    end
  end

  return condense( SchubertSum( cfs, sch ) )

end


# overload '-' for SchubertSum
import Base: -

function -(ss1::SchubertSum, ss2::SchubertSum) # difference of two Schubert sums

  return( +(ss1,-1*ss2) )

end


# overload '*' for SchubertSum
import Base: *

function *(cc::Int, ss::SchubertSum) # multiply Schubert sum by Int

  cc==0 && return condense( SchubertSum([0],[[1]]) )

  sch = ss.schubs
  cfs = cc*ss.coeffs

  return SchubertSum( cfs, sch )

end

function *(cc::ZZMPolyRingElem, ss::SchubertSum) # multiply Schubert sum by ZZMPolyRingElem

  cc==0 && return condense( SchubertSum([0],[[1]]) )

  sch = ss.schubs
  cfs = cc.*ss.coeffs

  return SchubertSum( cfs, sch )

end




function condense(ss::SchubertSum) # remove zeroes from Schubert sum, and trim indexing permutations

  ll = length(ss)
  sch = ss.schubs
  cfs = ss.coeffs


  i=1
  while i<=ll
     if cfs[i]==0
        cfs = vcat(cfs[1:i-1],cfs[i+1:ll])
        sch = vcat(sch[1:i-1],sch[i+1:ll])

        ll=ll-1
        continue
     else
        sch[i]=trimw(sch[i])

        if length(sch[i])==0
          push!(sch[i],1)
        end

        i+=1
     end
  end

  ssnew = SchubertSum( cfs, sch )

  return ssnew

end

