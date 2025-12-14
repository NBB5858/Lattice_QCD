#ifndef ENERGY_H
#define ENERGY_H

#include <array>

template <typename EnergyModel, typename MagModel, int D>
class Lattice;


struct IsingE1D {
    IsingE1D( double J, double B ) : B_(B), J_(J), size_x(0) {};

    template<typename E, typename M, int D>
    void bind_geometry( const Lattice<E, M, D> &L ) { size_x = L.shape()[0]; }

    template<typename E, typename M, int D>
    double total_energy( const Lattice<E, M, D> &L ) const {
        double total_energy = 0;
        for(int i = 0; i < size_x; ++i) total_energy += -J_ * L(i) * L(wrap_l(i-1)) - B_ * L(i);
        return total_energy;
    };

    template<typename E, typename M, int D>
    double dE( const Lattice<E, M, D> &L, std::array<int,D> site) const {
        int x = site[0];
        return  2*J_*L(x) * ( L(wrap_l(x-1)) + L(wrap_r(x+1) ) + 2*B_*L(x) ) ;
    };

    int wrap_l(int p) const {return (p == -1 ? size_x - 1 : p);}
    int wrap_r(int p) const {return (p == size_x ? 0 : p);}

    double B_, J_;
    int size_x;
};


struct IsingE2D {
    IsingE2D( double J, double B ) : B_(B), J_(J), size_x(0), size_y(0), size(0){};

    template<typename E, typename M, int D>
    void bind_geometry( const Lattice<E, M, D> &L ) {
        auto [nx, ny] = L.shape(); size_x=nx, size_y=ny; size=size_x*size_y;
    }

    template<typename E, typename M, int D>
    double total_energy( const Lattice<E, M, D> &L ) const {

        double total_energy = 0;
        for(int i = 0; i < size_x; ++i) {
            for(int j = 0; j < size_y; ++j) {
                total_energy += -J_ * L(i, j) * ( L(wrap_l(i-1), j) + L(i, wrap_u(j-1)) ) - B_ * L(i, j);
            }
        }
        return total_energy;
    };

    template<typename E, typename M, int D>
    double dE( const Lattice<E, M, D> &L, std::array<int,D> xy) const {
        auto [i, j] = xy;

        double dlink_energy = 2*J_*L(i, j) * ( L(wrap_l(i-1), j) + L(wrap_r(i+1), j) + L(i, wrap_u(j-1)) + L(i, wrap_d(j+1)) );
        double dmag_energy = 2*B_*L(i, j);

        return  dlink_energy + dmag_energy;
    };

    int wrap_l(int p) const {return (p == -1 ? size_x - 1 : p);}
    int wrap_r(int p) const {return (p == size_x ? 0 : p);}

    int wrap_u(int p) const {return (p == -1 ? size_y - 1 : p);}
    int wrap_d(int p) const {return (p == size_y ? 0 : p);}

    double B_, J_;
    int size_x, size_y;
    double size;

};

#endif //ENERGY_H