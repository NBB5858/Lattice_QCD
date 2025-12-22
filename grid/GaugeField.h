#ifndef GAUGEFIELD_H
#define GAUGEFIELD_H

#include <iostream>

#include "RNG.h"
#include "ScalarField.h"

class GroupElement {};

template <typename T> using is_group = std::is_base_of<GroupElement, T>;

class Z2GroupElement : public GroupElement {
public:
    Z2GroupElement() : _val(1){}
    explicit Z2GroupElement(int val) : _val(val) {};

    double val() const { return _val; }
    int groupdim() const { return _groupdim; }

    void Randomize(RandSign& rng) {_val = rng.draw();}

    friend Z2GroupElement operator*(const Z2GroupElement& g1, const Z2GroupElement& g2);
    friend Z2GroupElement adj(const Z2GroupElement& g);
    friend ScalarField trace(const Z2GroupElement& g);

    void print() const { std::cout << _val; }

private:
    int _val;
    int _groupdim = 1;
};

Z2GroupElement inline operator*(const Z2GroupElement& g1, const Z2GroupElement& g2) {return Z2GroupElement(g1._val + g2._val);}
Z2GroupElement inline adj(const Z2GroupElement& g) {return g;}
ScalarField inline trace(const Z2GroupElement& g) {return ScalarField(g._val);}



#endif //GAUGEFIELD_H
