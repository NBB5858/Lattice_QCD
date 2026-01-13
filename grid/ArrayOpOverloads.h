#ifndef ARRAYOPOVERLOADS_H
#define ARRAYOPOVERLOADS_H

#include "ColorMatrix.h"

#include <array>

template<typename S, typename T, std::size_t N,
         std::enable_if_t<std::is_arithmetic_v<S>, int> = 0>
std::array<T, N> operator*(S c, const std::array<T, N>& V) {
    std::array<T, N> ret{};
    for (std::size_t i = 0; i < N; ++i) {ret[i] = c * V[i];}
    return ret;
}

template<typename S, typename T, std::size_t N,
         std::enable_if_t<std::is_arithmetic_v<S>, int> = 0>
std::array<T, N> operator*(const std::array<T, N>& V, S c) {
    return c * V;
}

template<typename S, typename T, std::size_t N,
         std::enable_if_t<std::is_floating_point_v<S>, int> = 0>
std::array<T, N> operator/(const std::array<T, N>& V, S c) {
    std::array<T, N> ret{};
    for (std::size_t i = 0; i < N; ++i) {ret[i] = V[i] / c;}
    return ret;
}

template<typename T, std::size_t N>
std::array<T, N>& operator+=(std::array<T, N>& a, const std::array<T, N>& b) {
    for (std::size_t i = 0; i < N; ++i) a[i] += b[i];
    return a;
}

template<typename T, std::size_t N>
std::array<T, N>& operator-=(std::array<T, N>& a, const std::array<T, N>& b) {
    for (std::size_t i = 0; i < N; ++i) a[i] -= b[i];
    return a;
}

template<typename T, std::size_t N>
std::array<T, N> operator+(std::array<T, N> a, const std::array<T, N>& b) { return a += b; }

template<typename T, std::size_t N>
std::array<T, N> operator-(std::array<T, N> a, const std::array<T, N>& b) { return a -= b; }

template<typename T, std::size_t N>
std::array<T, N> operator*(const std::array<T, N>& V1, const std::array<T, N>& V2) {
    std::array<T, N> ret{};
    for (std::size_t i = 0; i < N; ++i) ret[i] = V1[i] * V2[i];
    return ret;
}




template<typename T, std::size_t d>
T Dot(const std::array<T, d>& V1, const std::array<T, d>& V2) {
    T ret;
    if constexpr (std::is_arithmetic_v<T> || std::is_same_v<T, std::complex<double>>) {
        ret = T{0};
    } else {
        ret = T::Zero();
    }

    for (std::size_t i = 0; i < d; ++i) ret += V1[i] * V2[i];
    return ret;
}

template<class Mat, std::size_t N>
inline std::array<Mat, N> adj(const std::array<Mat, N>& V)
{
    std::array<Mat, N> ret{};
    for (std::size_t i = 0; i < N; ++i) ret[i] = adj(V[i]);
    return ret;
}

template<class Mat, std::size_t N>
inline std::array<std::complex<double>, N> trace(const std::array<Mat, N>& V)
{
    std::array<std::complex<double>, N> ret{};
    for (std::size_t i = 0; i < N; ++i) {
        ret[i] = trace(V[i]);
    }
    return ret;
}

template<class Mat, std::size_t N>
inline std::array<Mat, N> TA(const std::array<Mat, N>& V) {
    std::array<Mat, N> ret{};
    for (std::size_t i = 0; i < N; ++i) {
        ret[i] = TA(V[i]);
    }
    return ret;
}

template<class Mat, std::size_t N>
inline std::array<Mat, N> U1exp(const std::array<Mat, N>& V) {
    std::array<Mat, N> ret{};
    for (std::size_t i = 0; i < N; ++i) {
        ret[i] = U1exp(V[i]);
    }
    return ret;
}

#endif //ARRAYOPOVERLOADS_H
