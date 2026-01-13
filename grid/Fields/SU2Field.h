#ifndef SU2FIELD_H
#define SU2FIELD_H

#include <iostream>
#include <complex>
#include <array>

#include "FieldBase.h"
#include "../ColorMatrix.h"
#include "../Expression.h"

constexpr std::complex<double> I (0.0,1.0);

template<int d>
struct SU2Field {};

template<int d>
struct SU2Momentum {};

template<int d>
struct SU2FieldView {
    using Storage = std::array<CMat<2>, d>;

    Storage* v;
    Storage& val() { return *v; }

    template<typename RNG>
    void Randomize(RNG& rng) {
        for(int i=0; i<d; ++i) {
            (*v)[i] = ////'CMat<2>{std::exp(I * rng.draw())};
        }
    }

};

template<int d>
struct SU2FieldConstView {
    using Storage = std::array<CMat<2>, d>;

    const Storage* v;
    const Storage& val() const { return *v; }
    void print() ////const { std::cout << "("; for(int i=0; i<d; ++i){std::cout << (*v)[i] << ",";} std::cout << ")"; }
};

template<int d>
struct FieldTraits<SU2Field<d>> {
    using Storage   = std::array<CMat<2>, d>;
    using View      = SU2FieldView<d>;
    using ConstView = SU2FieldConstView<d>;
    static constexpr int groupdim = 2;

    using MomentumField = SU2Momentum<d>;

    static View view(Storage& x) { return View{ &x }; }
    static ConstView view(const Storage& x) { return ConstView{ &x }; }

    static void set_zero(Storage& x) { for(int i=0; i<d; ++i) x[i] = CMat<2>::Zero(); };
};


template<int d>
struct SU2MomentumView {
    using Storage = std::array<CMat<2>, d>;

    Storage* v;

    Storage& val() { return *v; }

    template<typename RNG>
    void Randomize(RNG& rng) {
        for (int i = 0; i < d; ++i) {
            (*v)[i] = CMat<2>{ std::complex<double>{ 0.0, rng.draw()} };
        }
    }
};

template<int d>
struct SU2MomentumConstView {
    using Storage = std::array<CMat<2>, d>;

    const Storage* v;
    const Storage& val() const { return *v; }
    void print() const ///{ std::cout << "("; for(int i=0; i<d; ++i){std::cout << (*v)[i] << ",";} std::cout << ")"; }
};

template<int d>
struct FieldTraits<SU2Momentum<d>> {
    using Storage   = std::array<CMat<2>, d>;
    using View      = SU2MomentumView<d>;
    using ConstView = SU2MomentumConstView<d>;

    static View view(Storage& x) { return View{ &x }; }
    static ConstView view(const Storage& x) { return ConstView{ &x }; }

    static void set_zero(Storage& x) { for(int i=0; i<d; ++i) x[i] = CMat<2>::Zero(); };

    template<class arg,
            std::enable_if_t<is_lattice<arg>::value || is_lattice_expr<arg>::value, int> = 0>
    static auto exp(const arg& A) -> decltype(SU2exp(A)) { return SU2exp(A); }
};

#endif //SU2FIELD_H
