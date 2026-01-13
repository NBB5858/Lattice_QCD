
#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <vector>
#include <array>
#include <cassert>

#include "Grid.h"
#include "Expression.h"
#include "RNG.h"
#include "ColorMatrix.h"
#include "Fields/U1Field.h"

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
    aligned_vector<Storage> _mem{};

public:
    GridBase<d>* grid() const { return _grid; }

    explicit Lattice(GridBase<d>* grid) :
    _grid(grid),
    _mem(static_cast<size_t>(grid->_osites))
    {
        for (auto& x : _mem) FieldTraits<FieldType>::set_zero(x);
    }


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
              std::enable_if_t<is_lattice<T1>::value || is_lattice_expr<T1>::value, int> = 0>
    inline Lattice& operator+=(const T1& expr) {
        return *this = (*this + expr);
    }

    template <typename T1,
              std::enable_if_t<is_lattice<T1>::value || is_lattice_expr<T1>::value, int> = 0>
    inline Lattice& operator-=(const T1& expr) {
        return *this = (*this - expr);
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

    friend Storage StapleSum(int ss, const Lattice& U) {

        Storage staple_sum;
        FieldTraits<FieldType>::set_zero(staple_sum);

        const int x = ss;

        for (int mu = 0; mu < d; ++mu) {
            for (int nu = 0; nu < d; ++nu) {
                if (nu == mu) continue;

                const int x_plus_mu = U._grid->_stencil.neighbors[x].p[mu];
                const int x_plus_nu = U._grid->_stencil.neighbors[x].p[nu];
                const int x_minus_nu = U._grid->_stencil.neighbors[x].m[nu];

                staple_sum[mu] += U._mem[x][nu]
                                * U._mem[x_plus_nu][mu]
                                * adj(U._mem[x_plus_mu][nu]);


                const int x_minus_nu_plus_mu = U._grid->_stencil.neighbors[x_minus_nu].p[mu];
                staple_sum[mu] += adj(U._mem[x_minus_nu][nu])
                                * U._mem[x_minus_nu][mu]
                                * U._mem[x_minus_nu_plus_mu][nu];
            }
        }

        return staple_sum;
    }

    friend Storage exp(int ss, const Lattice& P) {
        return FieldTraits<FieldType>::exp(P._mem[ss]);
    }

};

#endif //LATTICE_H