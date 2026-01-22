#ifndef U1FIELD_H
#define U1FIELD_H

#include <iostream>

#include "../simd/SimdFieldTypes.h"
#include "FieldBase.h"
#include "../RNG/RNG.h"
#include "../Lattice/Expression.h"

template<int d>
struct U1Field {};

template<int d>
struct U1Momentum {};

template<int d>
struct U1FieldView {
    Link<1, d>* blk = nullptr;
    int lane = 0;

    CMat<1> get(int mu) const { return simd::extract_lane((*blk)[mu], lane); }
    void set( CMat<1> v, int mu ) const {simd::insert_lane((*blk)[mu], lane, v); }
};

template<int d>
struct U1FieldConstView {
    const Link<1, d>* blk = nullptr;
    int lane = 0;

    std::complex<float> get(int mu) const {
        return simd::extract_lane((*blk)[mu][0][0], lane);
    }
    void print() const {
        std::cout << "(";
        for(int i=0; i<d; ++i) std::cout << get(i) << ",";
        std::cout << ")";
    }
};


template<int d>
struct FieldTraits<U1Field<d>> {
    using Storage = Link<1, d>;
    using View = U1FieldView<d>;
    using ConstView = U1FieldConstView<d>;

    static constexpr int groupdim = 1;

    using MomentumField = U1Momentum<d>;

    static View view(Storage& x, int lane) { return View{ &x, lane }; }
    static ConstView view(const Storage& x, int lane) { return ConstView{ &x, lane }; }

    static void Randomize(Link<1, d>& U, RandNormal& rng) {
        for(int i=0; i<d; ++i) {
            simd::vReal x = simd::rand_vReal(rng);
            simd::vReal y = simd::rand_vReal(rng);
            simd::vReal r = simd::sqrt(x*x + y*y);

            simd::vComplex u1el{x/r, y/r};
            U[i] = CMat<1>{u1el};
        }
    }
};

template<int d>
struct FieldTraits<U1Momentum<d>> {
    using Storage = Link<1, d>;
    using View = U1FieldView<d>;
    using ConstView = U1FieldConstView<d>;

    static View view(Storage& x, int lane) { return View{ &x, lane }; }
    static ConstView view(const Storage& x, int lane) { return ConstView{ &x, lane }; }

    static void Randomize(Link<1, d>& P, RandNormal& rng) {
        for (int i = 0; i < d; ++i) {
            vComplex itheta{simd::vReal(0.0f), simd::rand_vReal(rng)};
            P[i] = CMat<1>{itheta};
        }
    }

    template<class arg,
             std::enable_if_t<is_lattice<arg>::value || is_lattice_expr<arg>::value, int> = 0>
    static auto exp(const arg& A) -> decltype(U1exp(A)) {
        return U1exp(A);
    }
};

#endif //U1FIELD_H