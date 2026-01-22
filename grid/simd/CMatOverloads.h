#ifndef CMATOVERLOADS_H
#define CMATOVERLOADS_H

#include "simd_types.h"
#include "SimdFieldTypes.h"

#include <array>

template<int Nc>
CMat<Nc> operator*(const vComplex& a, const CMat<Nc>& b) {
  CMat<Nc> result{};
  for (int i = 0; i < Nc; ++i) {
    for (int j = 0; j < Nc; ++j) {
      result[i][j] = a * b[i][j];
    }
  }
  return result;
}

template<int Nc>
CMat<Nc> operator*(const CMat<Nc>& b, const vComplex& a) {
  return a * b;
}

template<int Nc>
CMat<Nc> operator*(const vReal& a, const CMat<Nc>& b) {
  CMat<Nc> result{};
  for (int i = 0; i < Nc; ++i) {
    for (int j = 0; j < Nc; ++j) {
      result[i][j] = a * b[i][j];
    }
  }
  return result;
}

template<int Nc>
CMat<Nc> operator*(const CMat<Nc>& b, const vReal& a) {
  return a * b;
}

template<int Nc>
CMat<Nc> operator*(const float a, const CMat<Nc>& b) {
  return vReal{a} * b;
}

template<int Nc>
CMat<Nc> operator*(const CMat<Nc>& b, const float a) {
  return vReal{a} * b;
}

template<int Nc>
CMat<Nc> operator*(const CMat<Nc>& a, const CMat<Nc>& b) {
  CMat<Nc> result{};
  for (int i = 0; i < Nc; ++i) {
     for (int k = 0; k < Nc; ++k) {
        const vComplex aik = a[i][k];
        for (int j = 0; j < Nc; ++j) {
        result[i][j] += aik * b[k][j];
      }
    }
  }
  return result;
}


template<int Nc>
CMat<Nc> operator+(const CMat<Nc>& a, const CMat<Nc>& b) {
  CMat<Nc> result{};
  for (int i = 0; i < Nc; ++i) {
    for (int j = 0; j < Nc; ++j) {
        result[i][j] = a[i][j] + b[i][j];
    }
  }
  return result;
}

template<int Nc>
CMat<Nc>& operator+=(CMat<Nc>& a, const CMat<Nc>& b) {
  for (int i = 0; i < Nc; ++i) {
    for (int j = 0; j < Nc; ++j) {
      a[i][j] += b[i][j];
    }
  }
  return a;
}

template<int Nc>
CMat<Nc> operator-(const CMat<Nc>& a, const CMat<Nc>& b) {
  CMat<Nc> result{};
  for (int i = 0; i < Nc; ++i) {
    for (int j = 0; j < Nc; ++j) {
      result[i][j] = a[i][j] - b[i][j];
    }
  }
  return result;
}

template<int Nc>
CMat<Nc>& operator-=(CMat<Nc>& a, const CMat<Nc>& b) {
  for (int i = 0; i < Nc; ++i) {
    for (int j = 0; j < Nc; ++j) {
      a[i][j] -= b[i][j];
    }
  }
  return a;
}

template<int Nc>
vComplex trace(const CMat<Nc>& a) {
  vComplex result{};
  for (int i = 0; i < Nc; ++i) {
    result += a[i][i];
  }
  return result;
}

template<int Nc>
CMat<Nc> adj(const CMat<Nc>& a) {
  CMat<Nc> result{};
  for (int i = 0; i < Nc; ++i) {
    for (int j = 0; j < Nc; ++j) {
      result[i][j] = simd::conj(a[j][i]);
    }
  }
  return result;
}

template <int Nc>
CMat<Nc> TA(const CMat<Nc>& M) {

  CMat<Nc> AH = 0.5 * (M - adj(M));
  if (Nc == 1) return AH;

  const vComplex tr = trace(AH);
  return AH - (tr / vReal{static_cast<float>(Nc)}) * CMat<Nc>::Identity();
}


inline CMat<1> U1exp(const CMat<1>& A) {
  CMat<1> result{};
  vReal theta = (A[0][0]).im;
  auto cos_theta = simd::vReal{simd::cos(theta)};
  auto sin_theta = simd::vReal{simd::sin(theta)};

  result[0][0] = vComplex{cos_theta, sin_theta};
  return result;
}

inline CMat<2> SU2exp(const CMat<2>& A) {
  CMat<2> result{};


  // A = [ i a3        a2 + i a1 ]
  //     [ -a2 + i a1  -i a3     ]
  vReal a1 = A[0][1].im;
  vReal a2 = A[0][1].re;
  vReal a3 = A[0][0].im;

  vReal alpha = simd::sqrt(a1*a1 + a2*a2 + a3*a3);

  auto cos_a = simd::vReal{ simd::cos(alpha) };
  auto sin_a = simd::vReal{ simd::sin(alpha) };

  auto eps = vReal(1e-20f);
  vReal f = sin_a / (alpha + eps);

  result[0][0] = vComplex{cos_a, vReal(0.0f)} + A[0][0] * f;
  result[1][1] = vComplex{cos_a, vReal(0.0f)} + A[1][1] * f;
  result[0][1] = A[0][1] * f;
  result[1][0] = A[1][0] * f;

  return result;
}

#endif //CMATOVERLOADS_H
