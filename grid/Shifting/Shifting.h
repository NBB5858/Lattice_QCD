#ifndef SHIFTING_H
#define SHIFTING_H

#include "../simd/simd_types.h"
#include "../simd/SimdFieldTypes.h"

template<typename FieldType, int d> class Lattice;
template<int d> class GridBase;

inline int mod_wrap(int x, int m) {
    int r = x % m;
    return (r < 0) ? (r + m) : r;
}


inline void copy_lane(vReal& dst, int lane_dst, const vReal& src, int lane_src) {
    float val = simd::extract_lane(src, lane_src);
    simd::insert_lane(dst, lane_dst, val);
}

inline void copy_lane(vComplex& dst, int lane_dst, const vComplex& src, int lane_src) {
    std::complex<float> val = simd::extract_lane(src, lane_src);
    simd::insert_lane(dst, lane_dst, val);
}

template<typename T, std::size_t N>
inline void copy_lane(std::array<T, N>& dst, int lane_dst, const std::array<T, N>& src, int lane_src) {
    for (std::size_t i = 0; i < N; ++i) {
        copy_lane(dst[i], lane_dst, src[i], lane_src);
    }
}

template<int Nc>
inline void copy_lane(CMat<Nc>& dst, int lane_dst, const CMat<Nc>& src, int lane_src) {
    for (int r = 0; r < Nc; ++r) {
        for (int c = 0; c < Nc; ++c) {
            copy_lane(dst.e[r][c], lane_dst, src.e[r][c], lane_src);
        }
    }
}

template<int Nc, int dd>
inline void copy_lane(Link<Nc, dd>& dst, int lane_dst, const Link<Nc, dd>& src, int lane_src) {
    for (int mu = 0; mu < dd; ++mu) {
        copy_lane(dst.U[mu], lane_dst, src.U[mu], lane_src);
    }
}

template<typename FieldType, int d>
auto Cshift(const Lattice<FieldType, d>& lat,
            std::uint64_t ss,
            const std::array<int, static_cast<std::size_t>(d)>& disp_sites)
{
    using Ret    = std::decay_t<decltype(lat(ss))>;
    using Coords = std::array<int, static_cast<std::size_t>(d)>;

    bool all_zero = true;
    for (int a = 0; a < d; ++a) all_zero &= (disp_sites[a] == 0);
    if (all_zero) return lat(ss);

    const auto* g  = lat.grid();
    const auto& st = g->stencil();

    constexpr int pack_dim = 0;
    constexpr int W = simd::W;

    const auto dims_sites = st.dims_sites();
    const Coords cb = st.ravel_block(static_cast<int>(ss));

    Ret out{};

    for (int lane = 0; lane < W; ++lane) {
        Coords cs = cb;
        cs[pack_dim] = cb[pack_dim] * W + lane;

        for (int a = 0; a < d; ++a) {
            cs[a] = mod_wrap(cs[a] + disp_sites[a], dims_sites[a]);
        }

        Coords tb = cs;
        const int lane_src = tb[pack_dim] % W;
        tb[pack_dim] /= W;

        const int ss_src = st.unravel_block(tb);
        const Ret src_block = lat(static_cast<std::uint64_t>(ss_src));

        copy_lane(out, lane, src_block, lane_src);
    }

    return out;
}
template<int d>
static std::array<int, static_cast<std::size_t>(d)>
unit_disp(int mu, int sgn) {
    std::array<int, static_cast<std::size_t>(d)> disp{};
    disp.fill(0);
    disp[mu] = sgn;
    return disp;
}

template<int d>
static std::array<int, static_cast<std::size_t>(d)>
disp2(int a, int sa, int b, int sb) {
    auto disp = unit_disp<d>(a, sa);
    disp[b] += sb;
    return disp;
}

#endif // SHIFTING_H




