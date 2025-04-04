# Tools for drawing BPDs in Julia
# David Anderson, October 2023.


export draw_bpd, print_all_bpds, print_all_Kbpds


# updated to use integer matrices
# "O" => 0
# "+" => 1
# "/" => 2
# "%" => 3
# "|" => 4
# "-" => 5
# "." => 6
# "*" => 7
# "" => 8
# "o" => 9




"""
    draw_bpd( b::Union{BPD,Drift}; saveto::String="none" img_size=(300,300), visible::Bool=true )

Display the bumpless pipedream `b`, and optionally save it to an image file `saveto`.

## Arguments
- `b::BPD`: a BPD
- `saveto::String`: the filename, with suffix specifying format.  (E.g., .png, .pdf)  Default is "none" for no file saved.
- `img_size`: an ordered pair specifying the image size.
- `visible::Bool` toggle whether the plot is displayed.  Default to `true`.

## Returns
`plot`: a plot object

## Example
```julia-repl
# Generate a BPD plot
julia> w = [1,4,5,3,2];

julia> b = Rothe(w)

/ - - - - 
| O O / - 
| O O | / 
| O / + + 
| / + + + 


julia> draw_bpd( b, saveto="bpd1.png" );

# Generate a Drift plot
julia> b = makeflat(b)

O O / - - 
O O | / - 
O / % | / 
/ % / + + 
| / + + + 


julia> dc = bpd2drift(b);

julia> draw_bpd( dc, saveto="dc1.png" )

# Also displays labeled drift diagrams
julia> dc = markconfig(dc)

julia> draw_bpd( dc, saveto="dc2.png" )
```
"""
function draw_bpd( B::Union{BPD,Drift} ; saveto::String="none", img_size=(300,300), visible::Bool=true )

    n, m = size(B.m)

  # set up plot
    p = plot(; xlim=(0, n), ylim=(0, m), aspect_ratio=:equal, legend=false, grid=true, framestyle=:none, tick_direction=:none, size=img_size)

  # light grid
    for i=1:n-1
      plot!([0,n],[i,i], linecolor=:black, linewidth=.25 )
    end
    for j=1:n-1
      plot!([j,j],[0,n], linecolor=:black, linewidth=.25 )
    end

  # place tiles    
    for i = 1:n
        for j = 1:n
            y, x = n-i, j-1  # Transpose and invert the y-coordinate
            tile!(p, B.m[i, j], x, y )
        end
    end

  # frame it
    plot!(; framestyle=:box, linecolor=:black, linewidth=3, ticks=nothing)

  # save to file
    if  saveto!="none"
      savefig(saveto)
    end

  # display
    if visible
      return p
    end
end



function print_all_bpds(w,fn,fmt)

  i=0

  for b in all_bpds(w)
    i+=1
    draw_bpd(b,saveto=string(fn,i,fmt),visible=false)
  end

end


function print_flat_bpds(w,fn,fmt)

  i=0

  for b in flat_bpds(w)
    i+=1
    draw_bpd(ds[i],saveto=string(fn,i,fmt),visible=false)
  end

end


function print_all_Kbpds(w,fn,fmt)

  i=0

  for b in all_Kbpds(w)
    i+=1
    draw_bpd(b,saveto=string(fn,i,fmt),visible=false)
  end

end



function draw_se_elbow_curve(x1, y1, x3, y3)
    x2, y2 = x3, y1
    t = range(0, stop=1, length=100)
    x_vals = @. (1-t)^2 * x1 + 2*(1-t)*t * x2 + t^2 * x3
    y_vals = @. (1-t)^2 * y1 + 2*(1-t)*t * y2 + t^2 * y3
    return x_vals, y_vals
end

function draw_nw_elbow_curve(x1, y1, x3, y3)
    x2, y2 = x1, y3
    t = range(0, stop=1, length=100)
    x_vals = @. (1-t)^2 * x1 + 2*(1-t)*t * x2 + t^2 * x3
    y_vals = @. (1-t)^2 * y1 + 2*(1-t)*t * y2 + t^2 * y3
    return x_vals, y_vals
end


function tile!(p, aa, xx, yy )
# insert the tile corresponding to entry aa at position (xx,yy) in plot p

  if aa==0  # "O"
    plot!(p,[xx, xx+1, xx+1, xx, xx], [yy, yy, yy+1, yy+1, yy], linecolor=:orange, linewidth=2, seriestype=:shape, fillalpha=0)

  elseif aa==1  # "+"
    plot!(p,[xx+0.5, xx+0.5], [yy, yy+1], linecolor=:blue, linewidth=2)
    plot!(p,[xx, xx+1], [yy+0.5, yy+0.5], linecolor=:blue, linewidth=2)

  elseif aa==4  # "|"
    plot!([xx+0.5, xx+0.5], [yy, yy+1], linecolor=:blue, linewidth=2)

  elseif aa==5  # "-"
    plot!([xx, xx+1], [yy+0.5, yy+0.5], linecolor=:blue, linewidth=2)

  elseif aa==2  # "/"
    x_vals, y_vals = draw_se_elbow_curve(xx+1, yy+0.5, xx+0.5, yy)
    plot!(x_vals, y_vals, linecolor=:blue, linewidth=2)

  elseif aa==3  # "%"
    x_vals, y_vals = draw_nw_elbow_curve(xx+0.5, yy+1, xx, yy+0.5)
    plot!(x_vals, y_vals, linecolor=:blue, linewidth=2)

  elseif aa==6 || aa==7  # "." or "*"
    scatter!([(xx + 0.5)], [(yy + 0.5)], markercolor=:blue, markersize=2, markerstrokewidth=0)

  elseif isa(aa, Tuple)
    plot!([xx, xx+1, xx+1, xx, xx], [yy, yy, yy+1, yy+1, yy], linecolor=:orange, linewidth=2, seriestype=:shape, fillalpha=0)
    if aa[2]
      annotate!([(xx + 0.5, yy + 0.5, text(string(aa[1]), :center, 10, :red))])
    else
      annotate!([(xx + 0.5, yy + 0.5, text(string(aa[1]), :center, 10))])
    end
  end

  if aa==7
    plot!([x, x+1, x+1, x, x], [y, y, y+1, y+1, y], linecolor=:orange, linewidth=2, seriestype=:shape, fillalpha=0)
  end

end



