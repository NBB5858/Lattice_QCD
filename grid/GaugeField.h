#ifndef GAUGEFIELD_H
#define GAUGEFIELD_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <iostream>
#include <complex>
#include <cmath>


#include "RNG.h"
#include "Fields/ScalarField.h"

// using namespace std::literals::complex_literals;
//
// class GroupElement {};
//
// template <typename T> using is_group = std::is_base_of<GroupElement, T>;
//
// class Z2GroupElement : public GroupElement {
// public:
//     Z2GroupElement() : _val(1){}
//     explicit Z2GroupElement(int val) : _val(val) {};
//
//     double val() const { return _val; }
//     static constexpr int groupdim = 1;
//
//     void Randomize(RandSign& rng) {_val = rng.draw();}
//
//     friend Z2GroupElement operator*(const Z2GroupElement& g1, const Z2GroupElement& g2);
//     friend Z2GroupElement adj(const Z2GroupElement& g);
//     friend Z2GroupElement conj(const Z2GroupElement& g);
//     friend ScalarField trace(const Z2GroupElement& g);
//
//     void print() const { std::cout << _val; }
//
// private:
//     int _val;
// };
//
// Z2GroupElement inline operator*(const Z2GroupElement& g1, const Z2GroupElement& g2) {return Z2GroupElement(g1._val * g2._val);}
// Z2GroupElement inline adj(const Z2GroupElement& g) {return g;}
// Z2GroupElement inline conj(const Z2GroupElement& g) {return g;}
// ScalarField inline trace(const Z2GroupElement& g) {return ScalarField(g._val);}
//
// // this is actually a ring element, which I think is
// class U1GroupElement : public GroupElement {
// public:
//     U1GroupElement() : _val(1.0+0.0i){}
//     explicit U1GroupElement(std::complex<double> val) : _val(val) {}
//
//     std::complex<double> val() const { return _val; }
//     static constexpr int groupdim = 1;
//
//     void Randomize(RandSign& rng) {_val = rng.draw();}
//
//     friend U1GroupElement operator*(const U1GroupElement& g1, const U1GroupElement& g2);
//     friend std::complex<double> operator+(const U1GroupElement& g1, const U1GroupElement& g2);
//     friend U1GroupElement adj(const U1GroupElement& g);
//     friend U1GroupElement conj(const U1GroupElement& g);
//     friend std::complex<double> trace(const U1GroupElement& g);
//
//     void print() const { std::cout << _val; }
//
// private:
//     std::complex<double> _val;
// };
//
// U1GroupElement inline operator*(const U1GroupElement& g1, const U1GroupElement& g2) {return U1GroupElement(g1._val * g2._val);}
// std::complex<double> inline operator+(const U1GroupElement& g1, const U1GroupElement& g2) {return g1._val + g2._val;}
// U1GroupElement inline adj(const U1GroupElement& g) {return g;}
// U1GroupElement inline conj(const U1GroupElement& g) {return g;}
// std::complex<double> inline trace(const U1GroupElement& g) {return g._val;}

#endif //GAUGEFIELD_H
