#ifndef COLORMATRIX_H
#define COLORMATRIX_H


#include <Eigen/Dense>
#include <complex>

template<int Nc>
using CMat = Eigen::Matrix<std::complex<double>, Nc, Nc>;


// would like to put in a guard for type here; need to know
// the eigen expression base class?
template<typename T>
inline T adj(const T& A) {
    return A.adjoint();
}

template<typename T>
inline std::complex<double> trace(const T& A) {
    return A.trace();
}


template<int Nc>
inline CMat<Nc> Identity() {
    return CMat<Nc>::Identity();
}

template<int Nc>
inline CMat<Nc> Zero() {
    return CMat<Nc>::Zero();
}

template<typename T>
inline T TracelessAntiHermitian(const T& M) {
    int Nc = M.rows();
    T AH = 0.5 * (M - M.adjoint());              // anti-Hermitian
    std::complex<double> tr = AH.trace();
    return AH - (tr / static_cast<double>(Nc)) * T::Identity();
}

#endif //COLORMATRIX_H
