# Compatibility / deprecation shims.
#
# Drift methods have moved to DriftPolynomials.jl.  SchubertPolynomials no longer
# depends on drifts in either direction.  The names below used to be exported
# here; they are now unexported migration stubs that point users to the new home.
# They are intentionally *not* exported, so that `using SchubertPolynomials,
# DriftPolynomials` does not produce export-collision warnings.

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
