# Compatibility / deprecation shims.
#
# Drift methods have moved to DriftPolynomials.jl.  SchubertPolynomials no longer
# depends on drifts in either direction.  The names below used to be exported
# here; they are now unexported migration stubs that point users to the new home.
# They are intentionally *not* exported, so that `using SchubertPolynomials,
# DriftPolynomials` does not produce export-collision warnings.

"""
    schur_poly(la, ff, R::DoublePolyRing; double=length(R.y_vars)>0, kwargs...)

Legacy positional-ring adapter onto the keyword `schur_poly` in
SemistandardTableaux.jl.  Prefer `schur_poly(la, ff; ring=, double=)`.

The default `double=length(R.y_vars)>0` preserves the old behavior of this
positional form (a ring with y-variables yields the double/factorial Schur
polynomial).  When `double=false` the single Schur polynomial is returned in
`R.ring`, using only its x-variables even if `R` carries y-variables.
"""
function schur_poly(la, ff, R::DoublePolyRing; double::Bool=length(R.y_vars)>0, kwargs...)
    if double || isempty(R.y_vars)
        return SemistandardTableaux.schur_poly(la, ff; ring=R.ring, double=double, kwargs...)
    end
    # single polynomial in a ring that also has y-variables: build it in the
    # x-only subring, then embed into R.ring by mapping x-generators across
    xr, _ = polynomial_ring(base_ring(R.ring), ["x$i" for i in 1:length(R.x_vars)])
    p = SemistandardTableaux.schur_poly(la, ff; ring=xr, kwargs...)
    return evaluate(p, R.x_vars)   # map xr's generators onto R.x_vars, landing in R.ring
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
