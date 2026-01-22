#ifndef LINKOVERLOADS_H
#define LINKOVERLOADS_H

#include "SimdFieldTypes.h"

template<int Nc, int d>
Link<Nc, d> adj(const Link<Nc, d>& A) {
    Link<Nc, d> result{};
    for (int mu = 0; mu < d; ++mu) {
        result[mu] = adj(A[mu]);
    }
    return result;
}

template<int Nc, int d>
Link<Nc, d> operator+(const Link<Nc, d>& A, const Link<Nc, d>& B) {
    Link<Nc, d> result{};
    for(int mu=0; mu < d; ++mu) {
        result[mu] = A[mu] + B[mu];
    }
    return result;
}



template<int Nc, int d>
Link<Nc, d> operator-(const Link<Nc, d>& A, const Link<Nc, d>& B) {
    Link<Nc, d> result{};
    for(int mu=0; mu < d; ++mu) {
        result[mu] = A[mu] - B[mu];
    }
    return result;
}


template<int Nc, int d>
Link<Nc, d> operator*(const Link<Nc, d>& A, const Link<Nc, d>& B) {
    Link<Nc, d> result{};
    for(int mu=0; mu < d; ++mu) {
        result[mu] = A[mu] * B[mu];
    }
    return result;
}

template<int Nc, int d>
Link<Nc, d> operator*(const vComplex& a, const Link<Nc, d>& B) {
    Link<Nc, d> result{};
    for(int mu=0; mu < d; ++mu) {
        result[mu] = a * B[mu];
    }
    return result;
}

template<int Nc, int d>
Link<Nc, d> operator*(const Link<Nc, d>& B, const vComplex& a) {
    return a * B;
}

template<int Nc, int d>
Link<Nc, d> operator*(const vReal& a, const Link<Nc, d>& B) {
    Link<Nc, d> result{};
    for(int mu=0; mu < d; ++mu) {
        result[mu] = a * B[mu];
    }
    return result;
}

template<int Nc, int d>
Link<Nc, d> operator*(const Link<Nc, d>& B, const vReal& a) {
    return a * B;
}

template<int Nc, int d>
Link<Nc, d> operator*(const float a, const Link<Nc, d>& B) {
    Link<Nc, d> result{};
    for(int mu=0; mu < d; ++mu) {
        result[mu] = vReal{a} * B[mu];
    }
    return result;
}

template<int Nc, int d>
Link<Nc, d> operator*(const Link<Nc, d>& B, const float a) {
    return a*B;
}


template<int Nc, int d>
CMat<Nc> Dot(const Link<Nc, d>& A, const Link<Nc, d>& B) {
    CMat<Nc> result{};
    for(int mu=0; mu < d; ++mu) {
        result += A[mu] * B[mu];
    }
    return result;
}

template<int Nc, int d>
Link<Nc, d> TA(const Link<Nc, d>& A) {
    Link<Nc, d> result{};
    for (int mu = 0; mu < d; ++mu) {
        result[mu] = TA(A[mu]);
    }
    return result;
}

template<int d>
Link<1, d> U1exp(const Link<1, d>& A) {
    Link<1, d> result{};
    for (int mu = 0; mu < d; ++mu) {
        result[mu] = U1exp(A[mu]);
    }
    return result;
}

template<int d>
Link<2, d> SU2exp(const Link<2, d>& A) {
    Link<2, d> result{};
    for (int mu = 0; mu < d; ++mu) {
        result[mu] = SU2exp(A[mu]);
    }
    return result;
}

#endif //LINKOVERLOADS_H

