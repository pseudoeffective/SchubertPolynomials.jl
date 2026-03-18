# Bit-packed permutation infrastructure for fast Schubert polynomial computation
# Supports UInt64 (n <= 16, 4 bits/entry) and UInt128 (16 < n <= 25, 5 bits/entry)
# Values stored 0-indexed internally (subtract 1 before packing, add 1 after unpacking)

const PackedPerm = Union{UInt64, UInt128}

# ---- UInt64 packing (n <= 16, 4 bits per entry) ----

@inline function pack_perm64(w::Vector{Int})
    code = UInt64(0)
    for i in 1:length(w)
        code |= UInt64(w[i] - 1) << (4 * (i - 1))
    end
    return code
end

@inline function unpack_perm64(code::UInt64, n::Int)
    w = Vector{Int}(undef, n)
    for i in 1:n
        w[i] = Int((code >> (4 * (i - 1))) & 0xF) + 1
    end
    return w
end

@inline function get_val(code::UInt64, i::Int)
    return Int((code >> (4 * (i - 1))) & 0xF) + 1
end

@inline function swap_positions(code::UInt64, a::Int, b::Int)
    shift_a = 4 * (a - 1)
    shift_b = 4 * (b - 1)
    va = (code >> shift_a) & UInt64(0xF)
    vb = (code >> shift_b) & UInt64(0xF)
    code &= ~(UInt64(0xF) << shift_a)
    code &= ~(UInt64(0xF) << shift_b)
    code |= vb << shift_a
    code |= va << shift_b
    return code
end

# ---- UInt128 packing (16 < n <= 25, 5 bits per entry) ----

@inline function pack_perm128(w::Vector{Int})
    code = UInt128(0)
    for i in 1:length(w)
        code |= UInt128(w[i] - 1) << (5 * (i - 1))
    end
    return code
end

@inline function unpack_perm128(code::UInt128, n::Int)
    w = Vector{Int}(undef, n)
    for i in 1:n
        w[i] = Int((code >> (5 * (i - 1))) & 0x1F) + 1
    end
    return w
end

@inline function get_val(code::UInt128, i::Int)
    return Int((code >> (5 * (i - 1))) & 0x1F) + 1
end

@inline function swap_positions(code::UInt128, a::Int, b::Int)
    shift_a = 5 * (a - 1)
    shift_b = 5 * (b - 1)
    va = (code >> shift_a) & UInt128(0x1F)
    vb = (code >> shift_b) & UInt128(0x1F)
    code &= ~(UInt128(0x1F) << shift_a)
    code &= ~(UInt128(0x1F) << shift_b)
    code |= vb << shift_a
    code |= va << shift_b
    return code
end

# ---- Helper functions (generic over PackedPerm) ----

"""
Check if (a,b) is a Bruhat cover going UP: w(a) < w(b), and no value
strictly between w(a) and w(b) appears at a position strictly between a and b.
"""
function is_bruhat_cover_up(code::PackedPerm, a::Int, b::Int)
    wa = get_val(code, a)
    wb = get_val(code, b)
    wa >= wb && return false
    wb - wa == 1 && return true
    b - a == 1 && return true
    for m in a+1:b-1
        wm = get_val(code, m)
        if wa < wm && wm < wb
            return false
        end
    end
    return true
end

"""
Find the transition pair (r, s). Scans r from right to left.
For each r, finds the minimum value to the left, then finds the largest
s > r with min_left < w(s) < w(r). Returns (0,0) if w is dominant.
"""
function find_transition_rs(code::PackedPerm, n::Int)
    for r in n-1:-1:2
        wr = get_val(code, r)
        min_left = wr + 1
        for i in 1:r-1
            wi = get_val(code, i)
            if wi < min_left
                min_left = wi
            end
        end
        best_s = 0
        for s in r+1:n
            ws = get_val(code, s)
            if ws < wr && min_left < ws
                best_s = s
            end
        end
        if best_s > 0
            return (r, best_s)
        end
    end
    return (0, 0)
end

"""
Find the first position i (1-indexed) where i + w(i) <= n.
Returns 0 if none (permutation is w0).
"""
function find_cotrans_index(code::PackedPerm, n::Int)
    for i in 1:n
        wi = get_val(code, i)
        if i + wi <= n
            return i
        end
    end
    return 0
end

"""
Compute the Coxeter length (number of inversions) of a packed permutation.
"""
function perm_length(code::PackedPerm, n::Int)
    ell = 0
    for i in 1:n
        wi = get_val(code, i)
        for j in i+1:n
            if wi > get_val(code, j)
                ell += 1
            end
        end
    end
    return ell
end
