
#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>

#include "Grid.h"
#include "Expression.h"
#include "RNG.h"

class LatticeBase {};

template<typename FieldType, int d>
class VectorField;

template<class FieldType, int d>
class Lattice : public LatticeBase {
public:
    GridBase<d>* grid() const { return _grid; }

    explicit Lattice(GridBase<d>* grid) : _grid(grid) {
        std::array<int, d> dims = grid->_dims;
        _mem.reserve(_grid->_osites);
        for(int i = 0; i < _grid->_osites; ++i) _mem.emplace_back();
    };

    template<typename RNG>
    void Randomize(int ss, RNG& rng) {
        _mem[ss].Randomize(rng);
    }

    void put(int i, FieldType field) {
        _mem[i] = field;
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

    const FieldType& operator()(int ss) const { return _mem[ss]; }

    void FlatPrint() const {
        for( FieldType SiteField : _mem) { SiteField.print(); std::cout << " ";};
        std::cout << std::endl;
    }

    void Print2D() const {
        // Make sure this really is a 2D lattice
        assert(_grid->_dims.size() == 2);

        const int nx = _grid->_dims[0];
        const int ny = _grid->_dims[1];

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                std::vector coord = {i, j};
                int ss = _grid->_stencil.unravel(coord);  // coord -> site index
                _mem[ss].print();
                std::cout << " ";
            }
            std::cout << "\n";
        }
    }

    // make new friends idiom; creates a new grad for each specialization of lattice
    // use SINFAE to ensure this doesn't operate on U_{mu}?
    friend VectorField<FieldType, d> Grad( int ss, const Lattice& lat ) {

        VectorField<FieldType, d> ret;

        for(int dir=0; dir < d; ++dir) {
            int p = lat._grid->_stencil.neighbors[ss].p[dir];
            int m = lat._grid->_stencil.neighbors[ss].m[dir];

            ret[dir] = (lat._mem[p] - lat._mem[m]) / 2.0;
        }
        return ret;
    }

    friend FieldType Box( int ss, const Lattice& lat ) {

        FieldType ret;

        for(int dir=0; dir < d; ++dir) {
            int p = lat._grid->_stencil.neighbors[ss].p[dir];
            int m = lat._grid->_stencil.neighbors[ss].m[dir];

            ret += (lat._mem[p] - 2.0*lat._mem[ss] + lat._mem[m]) / 1.0;
        }
        return ret;
    }

    //
    // if the lattice type is a vector of gauge fields, then this should be defined. use SINFAE?

    // template<typename GaugeField>
    // GaugeField plaquette(int ss, int mu, int nu, const Lattice<VectorField<GaugeField>>& U, GridBase* grid ) {
    //
    //     int ss;
    //     int right = grid->_stencil.neighbors[ss].p[mu];
    //     int up = grid->_stencil.neighbors[ss].p[nu];
    //     int right_up = grid->_stencil.neighbors[right].p[nu];
    //
    //     // each site stores just a field, not a field and it's inverse
    //     U(ss)[mu] * U(right)[nu] * inverse(U(up)[mu]) * inverse(U(ss)[nu]);
    //
    // }








private:
    GridBase<d>* _grid;
    std::vector<FieldType> _mem;
};

#endif //LATTICE_H
