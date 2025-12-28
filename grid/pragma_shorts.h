#ifndef PRAGMA_SHORTS_H
#define PRAGMA_SHORTS_H

// uncomment to turn OFF parallelization
//#define GRID_FORCE_SERIAL

#include "Fields/ScalarField.h"

#ifdef _OPENMP
  #ifndef GRID_FORCE_SERIAL
    #define GRID_OMP
    #include <omp.h>
  #endif
#endif

#ifdef GRID_OMP
#define DO_PRAGMA_(x) _Pragma (#x)
#define DO_PRAGMA(x) DO_PRAGMA_(x)
#define thread_num(a) omp_get_thread_num()
#define thread_max(a) omp_get_max_threads()
#else
#define DO_PRAGMA_(x)
#define DO_PRAGMA(x)
#define thread_num(a) (0)
#define thread_max(a) (1)
#endif

#define thread_for( i, num, ... ) \
DO_PRAGMA(omp parallel for schedule(static)) for(uint64_t i=0; i<num; i++ ) {__VA_ARGS__};

#define thread_sum( i, num, out, ... ) \
DO_PRAGMA(omp parallel for reduction(+:out)) for(uint64_t i=0; i<num; i++) {__VA_ARGS__};

// DO_PRAGMA(omp declare reduction (field_sum:ScalarField:omp_out=omp_out+omp_in) initializer(omp_priv=ScalarField(0.0)))

// #define thread_field_sum( i, num, out, ... ) \
// DO_PRAGMA(omp parallel for reduction(field_sum:out)) for(uint64_t i=0; i<num; i++) {__VA_ARGS__};

#define thread_for2d( i1, n1, i2, n2, ...) \
DO_PRAGMA(omp parallel for collapse(2)) \
for(uint64_t i1=0; i1<n1; i1++ ) { \
for(uint64_t i2=0; i2<n2; i2++ ) { \
__VA_ARGS__; \
}; \
};





#endif //PRAGMA_SHORTS_H
