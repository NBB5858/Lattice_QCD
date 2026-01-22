#ifndef GRID_H
#define GRID_H

#include <array>

#include "openmp/pragma_shorts.h"
#include "Stencil/Stencil.h"

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

    explicit GridBase(const std::array<int, d>& dims)
        : _dims(dims), _stencil(dims)
    {
        _nsites = std::accumulate(_dims.begin(), _dims.end(), 1, std::multiplies<int>{});
        std::cout << "GridBase: total sites = " << _nsites << std::endl;

        constexpr int pack_dim = 0;

        if ((_dims[pack_dim] % W) != 0) {
            throw std::runtime_error("GridBase: dims[0] must be a multiple of simd::W when packing across x.");
        }

        const auto& db = _stencil.dims_blocks();
        _osites = std::accumulate(db.begin(), db.end(), 1, std::multiplies<int>{});
    }

    std::array<int, d> dims() const { return _dims; }
    const BlockStencil<d>& stencil() const { return _stencil; }
    int osites() const { return _osites; }
    int nsites() const { return _nsites; }

private:
    static constexpr int W = simd::W;
    std::array<int, d> _dims;
    BlockStencil<d> _stencil;
    int _osites = 0;
    int _nsites = 0;
};

#endif //GRID_H