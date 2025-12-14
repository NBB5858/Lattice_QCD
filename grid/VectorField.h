#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include <iostream>


template<typename FieldType>
class VectorField {
public:
    VectorField() = default;
    explicit VectorField(int d) : _val(d, FieldType()) {}
    VectorField(int d, FieldType el) : _val(d, el) {}
    VectorField(int d, double el) : _val(d) {
        for(int i=0; i<d; ++i) _val[i] = FieldType(el);
    }

    FieldType& operator[](int i) { return _val[i]; }

    const std::vector<FieldType>& val() const { return _val; }

    void print() { std::cout << "("; for( FieldType el : _val ){el.print(); std::cout << ",";} std::cout << ")"; }

    // static inline void generate_momenta(ScalarField& P, pRNG) {
    //     gaussian(pRNG, P);
    //     P = std::sqrt(2) * P;
    // }

    template<typename RNG>
    void Randomize(RNG& rng) {
        int d = this->_val.size();
        for(int i=0; i<d; ++i) _val[i].Randomize(rng);
    }

    VectorField& operator +=(const VectorField& other) {
        int d = this->_val.size();
        for(int i=0; i<d; ++i) this->_val[i] += other._val[i];
        return *this;
    }

    VectorField& operator -=(const VectorField& other) {
        int d = this->_val.size();
        for(int i=0; i<d; ++i) this->_val[i] -= other._val[i];
        return *this;
    }

    friend VectorField operator+(const VectorField&  V1, const VectorField& V2) {
        int d = V1._val.size();
        VectorField ret(d); for(int i=0; i < d; ++i) ret[i] = V1._val[i] + V2._val[i];
        return ret;
    }

    friend VectorField operator-(const VectorField&  V1, const VectorField& V2) {
        int d = V1._val.size();
        VectorField ret(d); for(int i=0; i < d; ++i) ret[i] = V1._val[i] - V2._val[i];
        return ret;
    }

    friend VectorField operator*(double c, const VectorField& V) {
        int d = V._val.size();
        VectorField ret(d); for(int i=0; i < d; ++i) ret[i] = c*V._val[i];
        return ret;
    }

    friend VectorField operator*(const VectorField& V, double c) {
        int d = V._val.size();
        VectorField ret(d); for(int i=0; i < d; ++i) ret[i] = c*V._val[i];
        return ret;
    }

    friend VectorField operator/(const VectorField& V, double c) {
        int d = V._val.size();
        VectorField ret(d); for(int i=0; i < d; ++i) ret[i] = c*V._val[i];
        return ret;
    }

    friend FieldType dot(const VectorField& V1, const VectorField& V2) {
        int d = V1._val.size();
        FieldType ret(0.0); for(int i=0; i < d; ++i) ret = ret +  V1._val[i] * V2._val[i];
        return ret;
    }

private:
    std::vector<FieldType> _val;
};



#endif //VECTORFIELD_H
