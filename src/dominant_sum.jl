# Type for writing Schubert classes in dominant basis
# David Anderson, June 2024.


export DominantSum, trimp

# TO DO: documentation


#########


function trimp( par::Vector{Int} )
# remove trailing zeros
  n=length(par)

  lam = deepcopy(par)
  while n>0 && lam[n]==0
    pop!(lam)
    n=length(lam)
  end

  return lam

end


# type for representing a sum of Schubert classes
mutable struct DominantSum{T<:Union{Int,ZZMPolyRingElem}}
    coeffs::Vector{T}
    pars::Vector{Vector{Int}}
end


# WARNING - so far no safety for possibly repeated elements in ss.pars
# always assume these are distinct


# a single Schubert class
function DominantSum( par::Vector{Int} )
  lam=trimp(par)
  length(lam)==0 && push!(lam,0)
  return DominantSum( [1], [lam] )
end


# FIXME meant to ensure length of coeffs and pars is the same
function DominantSum(coeffs::Vector{T}, pars::Vector{Vector{Int}}) where T <: Any
    cfs = Vector{Union{Int,ZZMPolyRingElem}}()
    for c in coeffs
      push!(cfs, c)
    end

    if length(cfs)==length(pars)
        return DominantSum(cfs, pars)

    elseif length(cfs) < length(pars)
        par = pars[1:length(cfs)]
        return DominantSum(cfs, par)

    elseif length(cfs) > length(pars)
        par = pars
        for i in length(pars)+1:length(cfs)
           push!(par,[0])
        end

        return DominantSum(cfs, par)
    end
  
end


function sort_dominant_sum( ss::DominantSum )

    p = sortperm( ss.pars )
    sorted_coeffs = ss.coeffs[p]
    sorted_pars = ss.pars[p]
    return DominantSum( sorted_coeffs, sorted_pars )

end


# overload identity for DominantSum
import Base: ==

function ==(a::DominantSum, b::DominantSum)

  a1 = condense(a)
  b1 = condense(b)

  a1 = sort_dominant_sum(a1)
  b1 = sort_dominant_sum(b1)

  return ( a1.pars == b1.pars && a1.coeffs == b1.coeffs )

end


# overload length for DominantSum
function Base.length(ss::DominantSum)
    length(ss.coeffs) 
end

# overload coeff for DominantSum
import Nemo.coeff
function coeff( ss::DominantSum, w::Vector{Int} )

  wi = findfirst( isequal(w), ss.pars )

  wi==nothing ? cc=0 : cc=ss.coeffs[wi]

  return cc
end


# overload show for DominantSum
#import Base.show

function Base.show(io::IO, p::DominantSum{T}) where T
    if isempty(p.coeffs)
        print(io, "0")
        return
    end
    
    first = true
    for (coeff, par) in zip(p.coeffs, p.pars)
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
        
        # Print the Schubert class, like ss[2,3,1]
        print(io, "ss", par)
    end
end


# overload '+' for DominantSum
import Base: +

function +(ss1::DominantSum, ss2::DominantSum) # add two Dominant sums

  length(ss1)>length(ss2) && return +( ss2, ss1 )

  par = copy( ss2.pars )
  cfs = copy( ss2.coeffs )

  for i in 1:length(ss1)

    vvi = findfirst(isequal(ss1.pars[i]),par)
    if vvi==nothing
      push!(par,ss1.pars[i])
      push!(cfs,ss1.coeffs[i])
    else
      cfs[vvi]=cfs[vvi]+ss1.coeffs[i]
    end
  end

  return condense( DominantSum( cfs, par ) )

end


# overload '-' for DominantSum
import Base: -

function -(ss1::DominantSum, ss2::DominantSum) # difference of two Dominant sums

  return( +(ss1,-1*ss2) )

end


# overload '*' for DominantSum
import Base: *

function *(cc::Int, ss::DominantSum) # multiply Dominant sum by Int

  cc==0 && return condense( DominantSum([0],[[0]]) )

  par = ss.pars
  cfs = cc*ss.coeffs

  return DominantSum( cfs, par )

end

function *(cc::ZZMPolyRingElem, ss::DominantSum) # multiply Dominant sum by ZZMPolyRingElem

  cc==0 && return condense( DominantSum([0],[[0]]) )

  par = ss.pars
  cfs = cc.*ss.coeffs

  return DominantSum( cfs, par )

end




function condense(ss::DominantSum) # remove zeroes from dominant sum, and trim indexing partitions

  ll = length(ss)
  par = ss.pars
  cfs = ss.coeffs


  i=1
  while i<=ll
     if cfs[i]==0
        cfs = vcat(cfs[1:i-1],cfs[i+1:ll])
        par = vcat(par[1:i-1],par[i+1:ll])

        ll=ll-1
        continue
     else
        par[i]=trimp(par[i])

        if length(par[i])==0
          push!(par[i],1)
        end

        i+=1
     end
  end

  ssnew = DominantSum( cfs, par )

  return ssnew

end

