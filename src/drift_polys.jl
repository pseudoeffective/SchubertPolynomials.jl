# Tools for computing drift polynomials
# David Anderson, May 2025.


export drift_poly, markconfig, dc2sd

# TO DO: review function drift2bin(d::Drift, R::DoublePolyRing) which acts like bpd2bin, records a product of binomials x[i]+y[j]
# TO DO: review function drift_poly(d::Drift, R::DoublePolyRing) which produces the drift polynomial from the iterator
# TO DO: clarify the logic in markbox
# TO DO: clarify the logic and clean up tableau_components

#####
# Drift polynomials
#####

function drift2bin( d::Drift, R::DoublePolyRing=xy_ring( size(d)[1], size(d)[2] )[1] )
# product of binomials for d
# requires DoublePolyRing
# can get single polyn by using no y_vars
  n,m=size(d)
  bin = R.ring(1)

  x = R.x_vars
  y = R.y_vars

  local aa=length(x)
  local bb=length(y)

  for i=1:n
    for j=1:m

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


function drift_poly( d::Drift, R::DoublePolyRing=xy_ring( size(d)[1] )[1]  )
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

  n,m = size(dc)

  dr=dc.m

  if dr[i,j]==7
    return(0,false,0)
  end

  if dr[i,j]!=0 && !isa(dr[i,j],Tuple)
    return dr[i,j]
  end

  if i==n || j==m
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
    while k2<k1 && i+k2+1<n && j+k2<m
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
    while k2<k1 && i+k2<n && j+k2+1<m
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

  local n,m=size(dc)

  local mm=Matrix{Any}(undef,n,m)

  for i=1:n
    for j=1:m
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

  local n,m=size(dc)

  local mm=Matrix{Any}(undef,n,m)

  for i=1:n
    for j=1:m
      mm[i,j]=unmarkbox(dc,i,j)
    end
  end

  return Drift(mm)

end



function schub_drifts( w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by drift class formula

  fbpds = flat_bpds(w)

  pol = R.ring(0)

  for b in fbpds
    b=markconfig(b)
    pol = pol+dc2sd(b,R)
  end

  return pol

end


function dc2sd( dc::Drift, R::DoublePolyRing=xy_ring( size(dc)[1], size(dc)[2] )[1]  )
# drift configuration to s-polynomial
# must take marked configuration as input

  local n,m=size(dc)

  for k=(n+m):-1:2
    for i=maximum([1,k-m]):minimum([n,k-1])
      if isa( dc.m[i,k-i], Tuple ) && dc.m[i,k-i][2]
        (dc1,dc2)=drift_split( dc, i, k-i )
        return ( dc2sd( dc1, R ) + dc2sd( dc2, R ) )
      end
    end
  end

  sd = R.ring(1)

  tcomps = tableau_components(dc)

  for tt in tcomps
    sd = sd*schur_poly( tt[1], tt[2], R; mu=tt[3], xoffset=tt[4][1], yoffset=tt[4][2], rowmin=true )
  end

  return sd

end





function tableau_components(dc::Drift)
# return labelled tableaux for a drift configuration dc
# must take marked configuration as input

  local n,m=size(dc)

  if !isflat(dc)
    return( tableau_components( nw_reset(dc) ) )
  end

  local lyds=Vector{Vector}([])

  local corners=Vector{Tuple{Int,Int}}([])

  for i=1:n
    for j=1:m
        if !( (i,j) in corners) && isa(dc.m[i,j],Tuple) && ((i,j)==(1,1) || (i>1 && j>1 && !isa( dc.m[i-1,j], Tuple) && !isa( dc.m[i,j-1],Tuple )  )) #find a new NW corner
        push!(corners,(i,j))

        local la=Vector{Int}([])
        local mu=Vector{Int}([])
        local rr=Vector{Int}([])

        local s=0
        while i+s<=n && isa( dc.m[i+s,j], Tuple )

          local k=0
          while j+k<=m && isa( dc.m[i+s,j+k], Tuple )  # find SE boxes
            k +=1
          end          
          push!(la,k)


          local kk=0
          while j-kk-1>0 && isa( dc.m[i+s,j-kk-1], Tuple )  # find skew boxes
            kk +=1
          end

          mu=mu+fill(kk,length(mu))
          la=la+fill(kk,length(la))
          push!(mu,0)

          if s>0 && i+s>1 && j-kk>1 && ( dc.m[i+s-1,j-kk-1]==1 || dc.m[i+s-1,j-kk-1]==9 || dc.m[i+s-1,j-kk-1]==6  )
            push!(corners,(i+s,j-kk) ) # record new corner
          end
          j=j-kk
          s +=1
        end

        rr=Vector{Vector{Int}}([])
        for el=1:length(la)
          push!(rr, fill(0,mu[el]) )
          for mm=mu[el]+1:la[el]
            push!(rr[end], dc.m[i-1+el,j-1+mm][1]+el )
          end
        end

        push!(lyds,[la,rr,mu,[i-1,j-1]])

      end
    end
  end

  return lyds
end



#=
# not used
# TO DO: bring all drift poly computations from schub_polys and put them here

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