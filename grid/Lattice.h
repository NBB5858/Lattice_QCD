
#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <Eigen/Dense>

#include "Grid.h"
#include "Expression.h"
#include "RNG.h"
#include "GaugeField.h"
#include "ColorMatrix.h"

class LatticeBase {};

template<typename FieldType, int d>
class VectorField;


template<class T>
using aligned_vector = std::vector<T, Eigen::aligned_allocator<T>>;


template<class FieldType, int d>
class Lattice : public LatticeBase {

    using Storage = typename FieldTraits<FieldType>::Storage;
    using View = typename FieldTraits<FieldType>::View;
    using ConstView = typename FieldTraits<FieldType>::ConstView;

    GridBase<d>* _grid;
    aligned_vector<Storage> _mem;

public:
    GridBase<d>* grid() const { return _grid; }

    // doesn't initialize correctly?
    explicit Lattice(GridBase<d>* grid) : _grid(grid) {
        _mem.resize(_grid->_osites);
    };


    const Storage& operator()(int ss) const { return _mem[ss]; }


    // Get a view to a site
    View site(int ss) { return FieldTraits<FieldType>::view(_mem[ss]); }
    ConstView site(int ss) const { return FieldTraits<FieldType>::view(_mem[ss]); }

    void FlatPrint() const {
        for (int ss = 0; ss < _grid->_osites; ++ss) {
            site(ss).print();
            std::cout << " ";
        }
        std::cout << "\n";
    }

    void Print2D() const {
        assert(_grid->_dims.size() == 2);

        const int nx = _grid->_dims[0];
        const int ny = _grid->_dims[1];

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                std::vector<int> coord = {i, j};
                int ss = _grid->_stencil.unravel(coord);
                site(ss).print();
                std::cout << " ";
            }
            std::cout << "\n";
        }
    }

    template<typename RNGType>
    void RandomizeLattice( std::vector<RNGType>& rngs ) {
        thread_for(ss, _grid->osites(), {
            int thread = thread_num(ss);
            Randomize(ss, rngs[thread]);
        });
    }

    template<typename RNG>
    void Randomize(int ss, RNG& rng) {
        site(ss).Randomize(rng);
    }

    void put(int i, const Storage& value) {
        _mem[i] = value;
    }

    template <typename T1,
              typename = std::enable_if_t<is_lattice_expr<T1>::value>>
    inline Lattice& operator=(const T1 &expr) {
        thread_for(ss, _grid->_osites, {
            auto tmp = eval(ss, expr);
            _mem[ss] = tmp;
        });
        return *this;
    };

    template <typename T1,
              typename std::enable_if<is_lattice<T1>::value || is_lattice_expr<T1>::value, int>::type = 0>
    inline Lattice& operator+=(const T1 &expr) {
        thread_for(ss, _grid->_osites, {
            auto tmp = eval(ss, expr);
            this->_mem[ss] += tmp;
        });
        return *this;
    };

    template <typename T1,
              typename std::enable_if<is_lattice<T1>::value || is_lattice_expr<T1>::value, int>::type = 0>
    inline Lattice& operator-=(const T1 &expr) {
        thread_for(ss, _grid->_osites, {
            auto tmp = eval(ss, expr);
            this->_mem[ss] -= tmp;
        });
        return *this;
    }

    friend std::array<Storage, d> Grad( int ss, const Lattice& lat ) {

        std::array<Storage, d> ret{};

        for(int dir=0; dir < d; ++dir) {
            int p = lat._grid->_stencil.neighbors[ss].p[dir];
            int m = lat._grid->_stencil.neighbors[ss].m[dir];

            ret[dir] = (lat._mem[p] - lat._mem[m]) / 2.0;
        }
        return ret;
    }

    friend Storage Box( int ss, const Lattice& lat ) {

        Storage ret{};

        for(int dir=0; dir < d; ++dir) {
            int p = lat._grid->_stencil.neighbors[ss].p[dir];
            int m = lat._grid->_stencil.neighbors[ss].m[dir];

            ret += (lat._mem[p] - 2.0*lat._mem[ss] + lat._mem[m]) / 1.0;
        }
        return ret;
    }

    friend double PlaqReTraceSum(int ss, const Lattice& U ) {
        double trace_sum = 0;
        for(int mu=0; mu<d; ++mu) {
            for(int nu=mu+1; nu<d; ++nu) {
                int right = U._grid->_stencil.neighbors[ss].p[mu];
                int up = U._grid->_stencil.neighbors[ss].p[nu];

                auto plaq = U._mem[ss][mu]
                      * U._mem[right][nu]
                      * adj(U._mem[up][mu])
                      * adj(U._mem[ss][nu]);

                trace_sum += std::real(trace(plaq));
            }
        }
        return trace_sum;
    }

    friend std::array<Storage, d> StapleSum(int ss, const Lattice& U ) {

        std::array<Storage, d> staple_sum;

        for(int mu=0; mu<d; ++mu) staple_sum[mu] = 0;

        for(int mu=0; mu<d; ++mu) {
            for(int nu=0; nu<d; ++nu) {
                if(nu==mu) continue;

                int x = ss;
                int x_plus_mu = U._grid->_stencil.neighbors[x].p[mu];
                int x_plus_nu = U._grid->_stencil.neighbors[x].p[nu];
                int x_minus_nu = U._grid->_stencil.neighbors[x].m[nu];

                staple_sum += U._mem[x][nu] * U._mem[x_plus_nu][mu] * adj(U._mem[x_plus_mu][nu]);

                int x_minus_nu_plus_mu = U._grid->_stencil.neighbors[x_minus_nu].p[mu];

                staple_sum += adj(U._mem[x_minus_nu][nu]) * U._mem[x_minus_nu][mu] * U._mem[x_minus_nu_plus_mu][nu];
            }
        }
        return staple_sum;
    }

    // Does not need to be extended unary. Fix.
    friend std::array<Storage, d> TA(int ss, const Lattice& U) {
        return TracelessAntiHermitian(ss);
    }

};


template<typename Storage, std::size_t d>
Storage Dot(const std::array<Storage, d>& V1, const std::array<Storage, d>& V2) {
    Storage ret{};
    for (int i = 0; i < d; ++i) ret += V1[i] * V2[i];
    return ret;
}

template<typename Storage, std::size_t d>
std::array<Storage, d> operator*(const std::array<Storage, d>& V1, const std::array<Storage, d>& V2) {
    std::array<Storage, d> ret{};
    for (int i = 0; i < d; ++i) ret[i] = V1[i] * V2[i];
    return ret;
}

#endif //LATTICE_H
