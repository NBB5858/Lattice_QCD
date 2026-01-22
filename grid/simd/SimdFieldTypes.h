#ifndef SIMDFIELDTYPES_H
#define SIMDFIELDTYPES_H

#include "simd_types.h"

using vComplex = simd::vComplex;
using vReal = simd::vReal;

template<int Nc>
struct CMat {
  std::array<std::array<vComplex, Nc>, Nc> e{};

  const std::array<vComplex, Nc>& operator[](std::size_t row) const { return e[row]; }
  std::array<vComplex, Nc>& operator[](std::size_t row) { return e[row]; }

  static CMat Identity() {
      CMat result{};
      for (std::size_t i = 0; i < Nc; ++i) {
        result[i][i] = vComplex{vReal{1.0f}, vReal{0.0f}};
      }
      return result;
  }};

template<int Nc, int d>
struct Link {
  std::array<CMat<Nc>, d> U{};

  const CMat<Nc>& operator[](std::size_t mu) const { return U[mu]; }
  CMat<Nc>& operator[](std::size_t mu) { return U[mu]; }

  int static size() { return d; }
};


#endif //SIMDFIELDTYPES_H
