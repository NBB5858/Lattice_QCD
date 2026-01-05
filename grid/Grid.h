#ifndef GRID_H
#define GRID_H

#include <array>

#include "pragma_shorts.h"
#include "Stencil.h"

class GridThread {
public:
    static int _threads;

    static void SetMaxThreads() {
#ifdef GRID_OMP
        _threads=omp_get_max_threads();
        omp_set_num_threads(_threads);
#else
        _threads = 1;
#endif
    };

    static int GetThreads() { return _threads; }
    static int ThreadBarrier() {
#ifdef GRID_OMP
#pragma omp barrier
        return omp_get_thread_num();
#else
        return 0;
#endif
    }
};


template<int d>
class GridBase : public GridThread {
public:
    template<class object, int dimension> friend class Lattice;
    explicit GridBase(const std::array<int, d>& dims) : _dims(dims), _stencil(dims) {
        _osites = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<>());
    }

    int osites() const { return _osites; }
    std::array<int, d> dims() const { return _dims; }

private:
    std::array<int, d> _dims;
    Stencil<d> _stencil;
    int _osites;
};


#endif //GRID_H