# SchubertPolynomials

*A Julia package for Schubert calculus and related combinatorics*

[![Build Status](https://github.com/pseudoeffective/SchubertPolynomials.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pseudoeffective/SchubertPolynomials.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package provides functions for computing [Schubert polynomials](https://en.wikipedia.org/wiki/Schubert_polynomial) via *bumpless pipe dreams* and *drift polynomials*.  The former ("BPDs") are gadgets which enumerate terms in a Schubert polynomial, developed by [Lam-Lee-Shimozono](https://arxiv.org/abs/1806.11233).  The latter are generalizations of Schur polynomials that William Fulton and I have found useful.  The package also uses "standard" ways of computing (divided difference operators, transition trees).

To use the package, first install Julia.  If you are new to this language, I recommend Ulrich Thiel's page on [JuLie](https://ulthiel.github.io/JuLie.jl/dev/), which in addition to providing software, has a nice introduction for the algebraically-minded mathematician.

To install the SchubertPolynomials package do this:
```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/pseudoeffective/SchubertPolynomials.jl")
```

To start using the package, type
```julia-repl
julia> using SchubertPolynomials
```

The package provides several methods for computing (double) Schubert and Grothendieck polynomials.  For example,
```julia-repl
julia> w = [1,3,4,2];

# set the ambient ring, this one just in three x variables
julia> Rx = xy_ring(3)[1];

julia> sp = schub_poly(w,Rx)
x1*x2 + x1*x3 + x2*x3

# for the double Schubert polynomial, change the ring to also include three y variables
julia> Rxy = xy_ring(3,3)[1];

julia> sp = schub_poly(w,Rxy)
x1*x2 + x1*x3 + x1*y1 + x1*y2 + x2*x3 + x2*y1 + x2*y2 + x3*y1 + x3*y2 + y1^2 + y1*y2 + y2^2
```

There is a function for quickly counting the terms of a Schubert polynomial:
```julia
julia> w=[1,4,3,2,10,9,8,7,6,5];

julia> nschub(w)
4424420
```

One can compute back stable Schubert polynomials, in the "dominant" basis `s[lambda]`, indexed by partitions:
```julia
julia> R = xy_ring(4,4)[1];

julia> back_schub_poly( [1,3,4,2], R )
(x2*x3 + x2*y2 + x3*y2 + y2^2)*ss[] + ss[1, 1] + (x3 + y2)*ss[1]
```
For more info, see the help files:
```julia
julia> ?schub_poly

julia> ?groth_poly

julia> ?back_schub_poly
```




There are also functions for generating Young tableaux and computing (double,flagged) Schur polynomials.
```julia
julia> la = [2,1]; ff=2;

julia> ssyt( la, ff )
2-element Vector{Tableau}:
 
1 1 
2 

 
1 2 
2 

# define a polynomial ring in five x variables and no y variables
julia> R = xy_ring(5)[1];

julia> sp = schur_poly( la, ff, R )
x1^2*x2+x1*x2^2

# the argument ff can be a vector of integers, for a flagged Schur polynomial:
julia> sp = schur_poly( la, [2,3], R)
x1^2*x2 + x1^2*x3 + x1*x2^2 + x1*x2*x3 + x2^2*x3

# use the keyword argument `mu` for a skew Schur polynomial:
julia> ssp = schur_poly( la, [2,3], R, mu=[1] )
x1^2 + 2*x1*x2 + x1*x3 + x2^2 + x2*x3

# changing the ambient ring to one with y variables produces a double Schur polynomial
julia> R = xy_ring(5,5)[1];

julia> sp2 = schur_poly( la, 2, R )
x1^2*x2 + x1^2*y1 + x1*x2^2 + 2*x1*x2*y1 + x1*x2*y2 + x1*x2*y3 + x1*y1^2 + x1*y1*y2 + x1*y1*y3 + x2^2*y1 + x2*y1^2 + x2*y1*y2 + x2*y1*y3 + y1^2*y2 + y1^2*y3
```

There are some functions for multiplying Schubert classes.  For example, to multiply inside $Fl(3)$,
```julia-repl
julia> mult_2schub( [2,1],[1,3,2],3 )
S[3, 1, 2] + S[2, 3, 1]
```

To multiply equivariant classes, specify a ring with enough y variables.
```julia-repl
julia> R,x,y = xy_ring(3,3);

julia> mult_2schub( [2,3,1], [3,1,2], 4, R )
(y1 - y3)*S[3, 2, 1] + S[4, 2, 1, 3]
```

You can also expand a polynomial in the Schubert basis.
```julia-repl
julia> f = evaluate( schub_poly([2,3,1],R), [y[1],y[2],y[3]],[-x[3],-x[2],-x[1]] )
x1*x2 - x1*x3 - x2*x3 + x3^2

julia> ss = expand_schub(f,3,R)
(y1*y2 - y1*y3 - y2*y3 + y3^2)*S[1] - (2*y1 - y2 - y3)*S[1, 3, 2] + 3*S[2, 3, 1]
```

And you can get the Schubert structure constant $c_{u,v}^w$ using `lrc`:
```julia-repl
julia> u = [2,3,1];
julia> v = [3,1,2];
julia> w = [4,2,1,3];
julia> lrc( u,v,w )
1
```
The default is non-equivariant.  To get equivariant coefficients, specify the ring and the ambient rank.
```julia-repl
julia> w = [3,2,1];
julia> lrc( u,v,w )
0
julia> R = xy_ring(3,3)[1]
julia> lrc( u,v,w, 3, R )
y1 - y3
```
The SchubertPolynomials package requires [BumplessPipeDreams](https://github.com/pseudoeffective/BumplessPipeDreams.jl).


*Thanks to [Anders Buch](https://sites.math.rutgers.edu/~asbuch/) for convincing me to try out Julia.  Some of the functions in this package are based on his Maple code for multiplying Schubert polynomials.*

*Thanks to [Hugh Dennin](https://www.asc.ohio-state.edu/dennin.3/index.html) for his contributions!* 

*Also check out Avery St. Dizier's [Mathematica package](https://github.com/avstdi/Schubert-Polynomial-Package) for Schubert polynomials.*

