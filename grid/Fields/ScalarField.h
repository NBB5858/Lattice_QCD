
#ifndef SCALARFIELD_H
#define SCALARFIELD_H

#include <iostream>

#include "FieldBase.h"

struct ScalarFieldView {
    double* v;
    double& val() { return *v; }
    const double& val() const { return *v; }

    template<typename RNG>
    void Randomize(RNG& rng) { *v = rng.draw(); }
};

struct ScalarFieldConstView {
    const double* v;
    const double& val() const { return *v; }
    void print() const { std::cout << *v; }
};


struct ScalarField {};
template<>
struct FieldTraits<ScalarField> {
    using Storage   = double;
    using View      = ScalarFieldView;
    using ConstView = ScalarFieldConstView;

    static View view(Storage& x) { return View{ &x }; }
    static ConstView view(const Storage& x) { return ConstView{ &x }; }
};


// class ScalarField {
// public:
//     ScalarField() : _val(0.0) {}
//     explicit ScalarField(double val) : _val(val) {}
//     void print() const { std::cout << _val; }
//
//     double val() const { return _val; }
//
//     template<typename RNG>
//     void Randomize(RNG& rng) {
//         _val = rng.draw();
//     }
//
//     ScalarField& operator +=(const ScalarField& other) {
//         this->_val += other._val;
//         return *this;
//     }
//
//     ScalarField& operator -=(const ScalarField& other) {
//         this->_val -= other._val;
//         return *this;
//     }
//
//     friend ScalarField operator+(const ScalarField&  phi1, const ScalarField& phi2);
//     friend ScalarField operator-(const ScalarField& phi1, const ScalarField& phi2);
//     friend ScalarField operator*(const ScalarField& phi1, const ScalarField& phi2);
//     friend ScalarField operator/(const ScalarField& phi1, const ScalarField& phi2);
//
//     friend ScalarField operator*(double c, const ScalarField& phi);
//     friend ScalarField operator*(const ScalarField& phi, double c);
//     friend ScalarField operator/(const ScalarField& phi, double c);
//
//     friend double index_sum(const ScalarField& phi);
//
// private:
//     double _val;
// };
//
// ScalarField inline operator+(const ScalarField& phi1, const ScalarField& phi2){return ScalarField(phi1._val + phi2._val);}
// ScalarField inline operator-(const ScalarField& phi1, const ScalarField& phi2){return ScalarField(phi1._val - phi2._val);}
// ScalarField inline operator*(const ScalarField& phi1, const ScalarField& phi2){return ScalarField(phi1._val * phi2._val);}
// ScalarField inline operator/(const ScalarField& phi1, const ScalarField& phi2){return ScalarField(phi1._val / phi2._val);}
//
// ScalarField inline operator*(double c, const ScalarField& phi){return ScalarField(c*phi._val);}
// ScalarField inline operator*(const ScalarField& phi, double c){return ScalarField(c*phi._val);}
// ScalarField inline operator/(const ScalarField& phi, double c){return ScalarField(phi._val/c);}
//
// double inline index_sum(const ScalarField& phi) { return phi._val; }
//


#endif //SCALARFIELD_H
