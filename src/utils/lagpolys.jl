function _lagpoly(nodes, i, s, i0=1, i1=length(nodes))
    T = eltype(s)
    r = one(T)
    si = nodes[i]
    for j in i0:i1
        j == i && continue
        sj = nodes[j]
        r *= (s - sj) / (si - sj)
    end
    return r
end

function _lagpoly_diff(nodes, i, s, i0=1, i1=length(nodes))
    T = typeof(s)
    r = zero(T)
    si = nodes[i]
    for p in i0:i1
        p == i && continue
        sp = nodes[p]
        rp = one(T)
        for j in i0:i1
            j == i && continue
            j == p && continue
            sj = nodes[j]
            rp *= (s - sj) / (si - sj)
        end
        rp *= 1 / (si - sp)
        r += rp
    end
    return r
end

function _sylpoly(nodes, i, s)
    _lagpoly(nodes, i, s, 1, i)
end

function _sylpoly_diff(nodes, i, s)
    _lagpoly_diff(nodes, i, s, 1, i)
end

function _sylpoly_shift(nodes, i, s)
    _lagpoly(nodes, i, s, 2, i)
end

function _sylpoly_shift_diff(nodes, i, s)
    _lagpoly_diff(nodes, i, s, 2, i)
end

@inline function _gwp_nodes(::Type{T}, ::Val{N}) where {T,N}
    SVector{N,T}(ntuple(i -> T(i-1) / T(N-1), Val(N)))
end

# Tables of all Sylvester polynomials (and derivatives) with index m = 1..M at x
@inline function _sylpolys(nodes, x, ::Val{M}) where {M}
    R   = SVector{M}(ntuple(m -> _sylpoly(nodes, m, x),            Val(M)))
    dR  = SVector{M}(ntuple(m -> _sylpoly_diff(nodes, m, x),       Val(M)))
    Rs  = SVector{M}(ntuple(m -> _sylpoly_shift(nodes, m, x),      Val(M)))
    dRs = SVector{M}(ntuple(m -> _sylpoly_shift_diff(nodes, m, x), Val(M)))
    return R, dR, Rs, dRs
end