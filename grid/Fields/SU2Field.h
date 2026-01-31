#ifndef SU2FIELD_H
#define SU2FIELD_H

#include <iostream>

#include "../simd/SimdFieldTypes.h"
#include "FieldBase.h"
#include "../RNG/RNG.h"
#include "../Lattice/Expression.h"


template<int d>
struct SU2Field {};

template<int d>
struct SU2Momentum {};

template<int d>
struct SU2FieldView {
    Link<2, d>* blk = nullptr;
    int lane = 0;

    CMat<2> get(int mu) const { return simd::extract_lane((*blk)[mu], lane); }
    void set(CMat<2> v, int mu) const { simd::insert_lane((*blk)[mu], lane, v); }
};

template<int d>
struct SU2FieldConstView {
    const Link<2, d>* blk = nullptr;
    int lane = 0;

    std::complex<float> get(int mu, int r, int c) const {
        return simd::extract_lane((*blk)[mu][r][c], lane);
    }

    void print_link(int mu) const {
        std::cout << "[";
        for (int r = 0; r < 2; ++r) {
            std::cout << "[";
            for (int c = 0; c < 2; ++c) {
                std::cout << get(mu, r, c);
                if (c + 1 < 2) std::cout << ",";
            }
            std::cout << "]";
            if (r + 1 < 2) std::cout << ",";
        }
        std::cout << "]";
    }
};


template<int d>
struct FieldTraits<SU2Field<d>> {
    using Storage   = Link<2, d>;
    using View      = SU2FieldView<d>;
    using ConstView = SU2FieldConstView<d>;

    static constexpr int groupdim = 2;

    using MomentumField = SU2Momentum<d>;

    static View view(Storage& x, int lane) { return View{ &x, lane }; }
    static ConstView view(const Storage& x, int lane) { return ConstView{ &x, lane }; }

    static void Randomize(Link<2, d>& U, RandNormal& rng) {
        for (int mu = 0; mu < d; ++mu) {

            simd::vReal a0 = simd::rand_vReal(rng);
            simd::vReal a1 = simd::rand_vReal(rng);
            simd::vReal a2 = simd::rand_vReal(rng);
            simd::vReal a3 = simd::rand_vReal(rng);

            simd::vReal r = simd::sqrt(a0*a0 + a1*a1 + a2*a2 + a3*a3);

            a0 = a0 / r; a1 = a1 / r; a2 = a2 / r; a3 = a3 / r;

            // Matrix:
            // [ a0 + i a3     a2 + i a1 ]
            // [ -a2 + i a1    a0 - i a3 ]
            vComplex u00{ a0,  a3 };
            vComplex u01{ a2,  a1 };
            vComplex u10{ -a2, a1 };
            vComplex u11{ a0, -a3 };

            U[mu] = CMat<2>{{
                {
                    { u00, u01 },
                   { u10, u11 }
                }
            }};
        }
    }
};

template<int d>
struct FieldTraits<SU2Momentum<d>> {
    using Storage   = Link<2, d>;
    using View      = SU2FieldView<d>;
    using ConstView = SU2FieldConstView<d>;

    static View view(Storage& x, int lane) { return View{ &x, lane }; }
    static ConstView view(const Storage& x, int lane) { return ConstView{ &x, lane }; }


    // P = i (p1 sigma1 + p2 sigma2 + p3 sigma3)
    //   = [ i p3        p2 + i p1 ]
    //     [ -p2 + i p1  -i p3     ]
    static void Randomize(Link<2, d>& P, RandNormal& rng) {
        for (int mu = 0; mu < d; ++mu) {

            simd::vReal p1 = simd::rand_vReal(rng);
            simd::vReal p2 = simd::rand_vReal(rng);
            simd::vReal p3 = simd::rand_vReal(rng);

            vComplex p00{ simd::vReal(0.0f),  p3 };
            vComplex p01{ p2,                p1 };
            vComplex p10{ -p2,               p1 };
            vComplex p11{ simd::vReal(0.0f), -p3 };

            P[mu] = CMat<2>{{
                {
                    { p00, p01 },
                   { p10, p11 }
                }
            }};
        }
    }

    template<class arg,
             std::enable_if_t<is_lattice<arg>::value || is_lattice_expr<arg>::value, int> = 0>
    static auto exp(const arg& A) -> decltype(SU2exp(A)) {
        return SU2exp(A);
    }
};



#endif //SU2FIELD_H
