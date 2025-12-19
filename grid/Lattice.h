
#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>

#include "Grid.h"
#include "Expression.h"
#include "RNG.h"

class LatticeBase {};

template<typename FieldType>
class VectorField;

template<class FieldType>
class Lattice : public LatticeBase {
public:
    GridBase* grid() const { return _grid; }
    //const GridBase* grid() const { return _grid; }

    explicit Lattice(GridBase* grid) : _grid(grid) {
        std::vector<int> dims = grid->_dims;
        _mem.reserve(_grid->_osites); // why can't this be reserve?
        for(int i = 0; i < _grid->_osites; ++i) _mem.emplace_back(FieldType(static_cast<int>(dims.size()))); // not used if FieldType is ScalarField
    };

    template<typename RNG>
    void Randomize( int ss, RNG& rng ) {
        _mem[ss].Randomize(rng);
    }

    // for debugging
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

//    FieldType operator()(int ss) const {return _mem[ss];}
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
                std::vector<int> coord = {i, j};
                int ss = _grid->_stencil.unravel(coord);  // coord -> site index
                _mem[ss].print();
                std::cout << " ";
            }
            std::cout << "\n";
        }
    }

    // make new friends idiom; creates a new grad for each specialization of lattice
    friend VectorField<FieldType> Grad( int ss, const Lattice& lat ) {

        int D = lat._grid->_dims.size();
        VectorField<FieldType> ret(D); // may rely on stupid constructor

        for(int dir=0; dir < D; ++dir) {
            int p = lat._grid->_stencil.neighbors[ss].p[dir];
            int m = lat._grid->_stencil.neighbors[ss].m[dir];

            ret[dir] = (lat._mem[p] - lat._mem[m]) / 2.0;
        }
        return ret;
    }

    friend FieldType Box( int ss, const Lattice& lat ) {

        int D = lat._grid->_dims.size();
        FieldType ret(lat._grid->_dims.size(), 0.0); // relies on stupid constructors

        for(int dir=0; dir < D; ++dir) {
            int p = lat._grid->_stencil.neighbors[ss].p[dir];
            int m = lat._grid->_stencil.neighbors[ss].m[dir];

            ret += (lat._mem[p] - 2.0*lat._mem[ss] + lat._mem[m]) / 1.0;
        }
        return ret;
    }

private:
    GridBase* _grid;
    std::vector<FieldType> _mem;
};

#endif //LATTICE_H
