#ifndef U1FIELD_H
#define U1FIELD_H

#include <iostream>
#include <complex>
#include <array>

#include "FieldBase.h"
#include "../ColorMatrix.h"
#include "../Expression.h"

constexpr std::complex<double> I (0.0,1.0);

template<int d>
struct U1Field {};

template<int d>
struct U1Momentum {};

template<int d>
struct U1FieldView {
    using Storage = std::array<CMat<1>, d>;

    Storage* v;
    Storage& val() { return *v; }

    template<typename RNG>
    void Randomize(RNG& rng) {
        for(int i=0; i<d; ++i) {
            (*v)[i] = CMat<1>{std::exp(I * rng.draw())};
        }
    }

};

template<int d>
struct U1FieldConstView {
    using Storage = std::array<CMat<1>, d>;

    const Storage* v;
    const Storage& val() const { return *v; }
    void print() const { std::cout << "("; for(int i=0; i<d; ++i){std::cout << (*v)[i] << ",";} std::cout << ")"; }
};

template<int d>
struct FieldTraits<U1Field<d>> {
    using Storage   = std::array<CMat<1>, d>;
    using View      = U1FieldView<d>;
    using ConstView = U1FieldConstView<d>;
    static constexpr int groupdim = 1;

    using MomentumField = U1Momentum<d>;

    static View view(Storage& x) { return View{ &x }; }
    static ConstView view(const Storage& x) { return ConstView{ &x }; }

    static void set_zero(Storage& x) { for(int i=0; i<d; ++i) x[i] = CMat<1>::Zero(); };
};


template<int d>
struct U1MomentumView {
    using Storage = std::array<CMat<1>, d>;

    Storage* v;

    Storage& val() { return *v; }

    template<typename RNG>
    void Randomize(RNG& rng) {
        for (int i = 0; i < d; ++i) {
            (*v)[i] = CMat<1>{ std::complex<double>{ 0.0, rng.draw()} };
        }
    }
};

template<int d>
struct U1MomentumConstView {
    using Storage = std::array<CMat<1>, d>;

    const Storage* v;
    const Storage& val() const { return *v; }
    void print() const { std::cout << "("; for(int i=0; i<d; ++i){std::cout << (*v)[i] << ",";} std::cout << ")"; }
};

template<int d>
struct FieldTraits<U1Momentum<d>> {
    using Storage   = std::array<CMat<1>, d>;
    using View      = U1MomentumView<d>;
    using ConstView = U1MomentumConstView<d>;

    static View view(Storage& x) { return View{ &x }; }
    static ConstView view(const Storage& x) { return ConstView{ &x }; }

    static void set_zero(Storage& x) { for(int i=0; i<d; ++i) x[i] = CMat<1>::Zero(); };

    template<class arg,
            std::enable_if_t<is_lattice<arg>::value || is_lattice_expr<arg>::value, int> = 0>
    static auto exp(const arg& A) -> decltype(U1exp(A)) { return U1exp(A); }
};

#endif //U1FIELD_H
