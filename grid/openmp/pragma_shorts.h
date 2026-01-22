#ifndef PRAGMA_SHORTS_H
#define PRAGMA_SHORTS_H

#include <cstdint>
#include "../simd/simd_types.h"


#if defined(USE_OPENMP) && defined(_OPENMP)
  #define GRID_OMP 1
  #include <omp.h>
#else
  #undef GRID_OMP
#endif

#ifdef GRID_OMP
  #define DO_PRAGMA_(x) _Pragma(#x)
  #define DO_PRAGMA(x)  DO_PRAGMA_(x)

  #define thread_num(a) omp_get_thread_num()
  #define thread_max(a) omp_get_max_threads()

  #pragma omp declare reduction( \
  + : simd::vReal : omp_out += omp_in ) \
  initializer( omp_priv = simd::vReal(0.0f) )

  #pragma omp declare reduction( \
  + : simd::vComplex : omp_out += omp_in ) \
  initializer( omp_priv = simd::vComplex::zero() )

#else
  #define DO_PRAGMA_(x)
  #define DO_PRAGMA(x)

  // Same call style, serial fallbacks
  #define thread_num(a) (0)
  #define thread_max(a) (1)
#endif

#define thread_for(i, num, ...) \
DO_PRAGMA(omp parallel for schedule(static)) \
for (std::uint64_t i = 0; i < (std::uint64_t)(num); ++i) { __VA_ARGS__ }

#define thread_sum(i, num, out, ...) \
DO_PRAGMA(omp parallel for schedule(static) reduction(+:out)) \
for (std::uint64_t i = 0; i < (std::uint64_t)(num); ++i) { __VA_ARGS__ }

#endif // PRAGMA_SHORTS_H