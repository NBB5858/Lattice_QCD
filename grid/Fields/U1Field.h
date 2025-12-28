#ifndef U1FIELD_H
#define U1FIELD_H

#include <iostream>
#include <complex>
#include <array>

#include "FieldBase.h"
#include "../ColorMatrix.h"

using std::complex;
using std::array;
using std::exp;

constexpr complex<double> I (0.0,1.0);

template<int d>
struct U1FieldView {
    std::array<CMat<1>, d>* v;
    std::array<CMat<1>, d>& val() { return *v; }
    const std::array<CMat<1>, d>& val() const { return *v; }

    template<typename RNG>
    void Randomize(RNG& rng) {
        for(int i=0; i<d; ++i) {
            CMat<1> tmp{exp(I * rng.draw())};
            (*v)[i] = tmp;
        }
    }

};

template<int d>
struct U1FieldConstView {
    const std::array<CMat<1>, d>* v;
    const std::array<CMat<1>, d>& val() const { return *v; }
    void print() const { std::cout << "("; for(int i=0; i<d; ++i){std::cout << (*v)[i] << ",";} std::cout << ")"; }
};


template<int d>
struct U1Field {};

template<int d>
struct FieldTraits<U1Field<d>> {
    using Storage   = std::array<CMat<1>, d>;
    using View      = U1FieldView<d>;
    using ConstView = U1FieldConstView<d>;

    static View view(Storage& x) { return View{ &x }; }
    static ConstView view(const Storage& x) { return ConstView{ &x }; }
};



#endif //U1FIELD_H
