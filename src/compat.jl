# Compatibility / deprecation shims.
#
# Drift methods have moved to DriftPolynomials.jl.  SchubertPolynomials no longer
# depends on drifts in either direction.  The names below used to be exported
# here; they are now unexported migration stubs that point users to the new home.
# They are intentionally *not* exported, so that `using SchubertPolynomials,
# DriftPolynomials` does not produce export-collision warnings.

"""
    schur_poly(la, ff, R::DoublePolyRing; double=length(R.y_vars)>0, kwargs...)

Legacy positional-`DoublePolyRing` adapter onto the keyword `schur_poly` in
SemistandardTableaux.jl.  Prefer `schur_poly(la, ff; ring=schub_ring(...), double=)`.

The default `double=length(R.y_vars)>0` preserves the old behavior of this
positional form (a ring with y-variables yielded the double/factorial Schur
polynomial).  Single vs double is otherwise controlled by `double`.
"""
function schur_poly(la, ff, R::DoublePolyRing; double::Bool=length(R.y_vars)>0, kwargs...)
    Base.depwarn(_RING_POSITIONAL, :schur_poly)
    SemistandardTableaux.schur_poly(la, ff; ring=R.ring, double=double, kwargs...)
end


################################################################################
# Deprecated positional / DoublePolyRing call forms.
#
# The ring is now a keyword argument and is a bare MPolyRing.  These shims accept
# the old positional `DoublePolyRing` and forward to the keyword API.
################################################################################

const _RING_POSITIONAL =
    "passing the ring positionally as a `DoublePolyRing` is deprecated; build the " *
    "ring with `schub_ring` and pass it as the `ring` keyword"

function schub_poly(w::Vector{Int}, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :schub_poly)
    schub_poly(w; ring=R.ring, kwargs...)
end

function groth_poly(w::Vector{Int}, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :groth_poly)
    groth_poly(w; ring=R.ring, kwargs...)
end

function back_schub_poly(w, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :back_schub_poly)
    back_schub_poly(w; ring=R.ring, kwargs...)
end

function acoeff(w::Vector{Int}, lambda::Vector{Int}, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :acoeff)
    acoeff(w, lambda; ring=R.ring, kwargs...)
end

function mult_2schub(uu::Vector{Int}, vv::Vector{Int}, rnk::Int, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :mult_2schub)
    mult_2schub(uu, vv, rnk; ring=R.ring, kwargs...)
end

function lrc(uu::Vector{Int}, vv::Vector{Int}, ww::Vector{Int}, rnk::Int, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :lrc)
    lrc(uu, vv, ww, rnk; ring=R.ring, kwargs...)
end

function mult_schub(f::Union{Int,ZZMPolyRingElem}, ss::SchubertSum, rnk::Int, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :mult_schub)
    mult_schub(f, ss, rnk; ring=R.ring, kwargs...)
end

function mult_schub(f::Union{Int,ZZMPolyRingElem}, w::Vector{Int}, rnk::Int, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :mult_schub)
    mult_schub(f, w, rnk; ring=R.ring, kwargs...)
end

function expand_schub(f::Union{Int,ZZMPolyRingElem}, rnk::Int, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :expand_schub)
    expand_schub(f, rnk; ring=R.ring, kwargs...)
end

function localization(u::Vector{Int}, v::Vector{Int}, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :localization)
    localization(u, v; ring=R.ring, kwargs...)
end

function principal_specialization(pol::ZZMPolyRingElem, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :principal_specialization)
    principal_specialization(pol; ring=R.ring, kwargs...)
end

function principal_specialization(w::Vector{Int}, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :principal_specialization)
    principal_specialization(w; ring=R.ring, kwargs...)
end

function schub2dom(SS::SchubertSum, R::DoublePolyRing; kwargs...)
    Base.depwarn(_RING_POSITIONAL, :schub2dom)
    schub2dom(SS; ring=R.ring, kwargs...)
end


const _DRIFT_MOVED =
    "drift methods have moved to DriftPolynomials.jl; add it to your project and " *
    "`using DriftPolynomials`"

for f in (:Drift, :drift_class, :drift_poly, :markconfig, :dc2sd,
          :partition2drift, :bpd2drift, :nw_reset, :se_reset, :random_drift,
          :isintegrable, :iscancelable, :drift2rkmtx, :rkmtx2asm, :empty_drift)
    @eval function $f(args...; kwargs...)
        error("`", $(string(f)), "` has been removed from SchubertPolynomials.jl; ",
              _DRIFT_MOVED)
    end
end
