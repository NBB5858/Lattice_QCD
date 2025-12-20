#ifndef GAUGEFIELD_H
#define GAUGEFIELD_H

#include <vector>
#include <iostream>

#include "RNG.h"

class Z2GroupElement {
public:
    Z2GroupElement() : _val(1){}
    explicit Z2GroupElement(int val) : _val(val) {};
    void print() const { std::cout << _val; }
    double val() const { return _val; }


    void Randomize(RandSign& rng) {
        _val = rng.draw();
    }

    friend Z2GroupElement operator*(const Z2GroupElement& g1, const Z2GroupElement& g2);
    friend Z2GroupElement adj(const Z2GroupElement& g);

private:
    int _val;
};

Z2GroupElement inline operator*(const Z2GroupElement& g1, const Z2GroupElement& g2) {return Z2GroupElement(g1._val + g2._val);}
Z2GroupElement inline adj(const Z2GroupElement& g) {return g;}


template<typename GroupType>
class GaugeField {

public:

    GaugeField() = default;
    explicit GaugeField(int d) : _val(d, GroupType()) {}
    GaugeField(int d, GroupType el) : _val(d, el) {}
    GaugeField(int d, double el) : _val(d) {
        for(int i=0; i<d; ++i) _val[i] = GroupType(el);
    }

    GroupType& operator[](int i) { return _val[i]; }

    const std::vector<GroupType>& val() const { return _val; }

    void print() { std::cout << "("; for( GroupType el : _val ){el.print(); std::cout << ",";} std::cout << ")"; }

    template<typename RNG>
        void Randomize(RNG& rng) {
        int d = this->_val.size();
        for(int i=0; i<d; ++i) _val[i].Randomize(rng);
    }

    friend GaugeField operator+(const GaugeField&  U1, const GaugeField& U2) {
        int d = U1._val.size();
        GroupType ret(d); for(int i=0; i < d; ++i) ret[i] = U1._val[i] * U2._val[i];
        return ret;
    }

private:
    std::vector<GroupType> _val;

};








#endif //GAUGEFIELD_H
