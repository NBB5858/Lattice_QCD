#ifndef GRID_H
#define GRID_H

#include "pragma_shorts.h"
#include "Stencil.h"

class GridThread {
public:
    static int _threads;

    static void SetMaxThreads() {
#ifdef GRID_OMP
        _threads = omp_get_max_threads();
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

class GridBase : public GridThread {
public:

    template<class object> friend class Lattice;
    explicit GridBase(const std::vector<int>& dims) : _dims(dims), _stencil(dims) {
        _osites = std::accumulate(begin(dims), end(dims), 1, std::multiplies<>());
    }

    int osites() const { return _osites; }

private:
    std::vector<int> _dims;
    Stencil _stencil;
    int _osites;
};


#endif //GRID_H