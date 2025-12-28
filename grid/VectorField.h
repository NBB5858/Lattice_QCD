#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include <array>
#include <type_traits>

template<typename T, int d>
class VectorField {
public:
    using value_type = T;

    VectorField() { _val.fill(T{}); }
    explicit VectorField(const T& el) { _val.fill(el); }

    T& operator[](int i) { return _val[i]; }
    const T& operator[](int i) const { return _val[i]; }

    const std::array<T, d>& val() const { return _val; }

    VectorField& operator+=(const VectorField& other) {
        for (int i = 0; i < d; ++i) _val[i] += other._val[i];
        return *this;
    }
    VectorField& operator-=(const VectorField& other) {
        for (int i = 0; i < d; ++i) _val[i] -= other._val[i];
        return *this;
    }

    friend VectorField operator+(VectorField a, const VectorField& b) { return a += b; }
    friend VectorField operator-(VectorField a, const VectorField& b) { return a -= b; }

    friend VectorField operator*(double c, const VectorField& V) {
        VectorField ret;
        for (int i = 0; i < d; ++i) ret[i] = c * V._val[i];
        return ret;
    }
    friend VectorField operator*(const VectorField& V, double c) { return c * V; }

    friend VectorField operator/(const VectorField& V, double c) {
        VectorField ret;
        for (int i = 0; i < d; ++i) ret[i] = V._val[i] / c;
        return ret;
    }

    template<typename U = T, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
    friend U dot(const VectorField& V1, const VectorField& V2) {
        U ret = U{};
        for (int i = 0; i < d; ++i) ret += V1._val[i] * V2._val[i];
        return ret;
    }

private:
    std::array<T, d> _val;
};



// #include <iostream>
//
//
// template<typename FieldType, int d>
// class VectorField {
// public:
//
//     VectorField() {for(int i=0; i<d; ++i) _val[i] = FieldType();}
//     explicit VectorField(FieldType el) {for(int i=0; i<d; ++i) _val[i] = el;}
//
//     FieldType& operator[](int i) { return _val[i]; }
//     const FieldType& operator[](int i) const { return _val[i]; }
//
//     const std::array<FieldType, d>& val() const { return _val; }
//
//     // void print() const { std::cout << "("; for( FieldType el : _val ){el.print(); std::cout << ",";} std::cout << ")"; }
//
//     // template<typename RNG>
//     // void Randomize(RNG& rng) {
//     //     for(int i=0; i<d; ++i) _val[i].Randomize(rng);
//     // }
//
//     VectorField& operator +=(const VectorField& other) {
//         for(int i=0; i<d; ++i) this->_val[i] += other._val[i];
//         return *this;
//     }
//
//     VectorField& operator -=(const VectorField& other) {
//         for(int i=0; i<d; ++i) this->_val[i] -= other._val[i];
//         return *this;
//     }
//
//     friend VectorField operator+(const VectorField&  V1, const VectorField& V2) {
//         VectorField ret; for(int i=0; i < d; ++i) ret[i] = V1._val[i] + V2._val[i];
//         return ret;
//     }
//
//     friend VectorField operator-(const VectorField&  V1, const VectorField& V2) {
//         VectorField ret; for(int i=0; i < d; ++i) ret[i] = V1._val[i] - V2._val[i];
//         return ret;
//     }
//
//     friend VectorField operator*(double c, const VectorField& V) {
//         VectorField ret; for(int i=0; i < d; ++i) ret[i] = c*V._val[i];
//         return ret;
//     }
//
//     friend VectorField operator*(const VectorField& V, double c) {
//         VectorField ret; for(int i=0; i < d; ++i) ret[i] = c*V._val[i];
//         return ret;
//     }
//
//     friend VectorField operator/(const VectorField& V, double c) {
//         VectorField ret; for(int i=0; i < d; ++i) ret[i] = V._val[i]/c;
//         return ret;
//     }
//
//     friend FieldType dot(const VectorField& V1, const VectorField& V2) {
//         FieldType ret(0.0); for(int i=0; i < d; ++i) ret = ret +  V1._val[i] * V2._val[i];
//         return ret;
//     }
//
// private:
//     std::array<FieldType, d> _val;
// };

// template<typename Inner, int d>
// auto trace(const VectorField<Inner, d>& V) {
//     using InnerTrace = decltype(trace(V[0]));
//     VectorField<InnerTrace, d> ret;
//
//     for (int i = 0; i < d; ++i) ret[i] = trace(V[i]);
//     return ret;
// }

#endif //VECTORFIELD_H
