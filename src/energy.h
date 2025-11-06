#ifndef ENERGY_H
#define ENERGY_H

#include <array>

template <typename EnergyModel, typename MagModel>
class Lattice;


class IsingE1D {
public:
    IsingE1D( double J, double B ) : B(B), J(J) {};

    template<typename E, typename M>
    double total_energy( const Lattice<E, M> &L ) const {
        int size_x = L.shape()[0];
        double total_energy = 0;
        for(int i = 0; i < size_x; ++i) {
            int left = (i == 0) ? size_x - 1 : i - 1;
            total_energy += -J * L(i) * L(left) - B * L(i);
        }
        return total_energy;
    };

    template<typename E, typename M>
    double dE( const Lattice<E, M> &L, std::array<int,3> xyz) const {
        int x = xyz[0];

        int left = (x == 0) ? L.shape()[0] - 1 : x - 1;
        int right = (x == L.shape()[0] - 1) ? 0 : x + 1;

        double dlink_energy = 2*J*L(x) * ( L(left) + L(right) );
        double dmag_energy = 2*B*L(x);

        return  dlink_energy + dmag_energy;
    };

private:
    double B;
    double J;
};


class IsingE2D {
public:
    IsingE2D( double J, double B ) : B(B), J(J) {};

    template<typename E, typename M>
    double total_energy( const Lattice<E, M> &L ) const {
        auto shape = L.shape();
        int x = shape[0];
        int y = shape[1];

        double total_energy = 0;
        for(int i = 0; i < x; ++i) {
            for(int j = 0; j < y; ++j) {
                int left = (i == 0) ? x - 1 : i - 1;
                int up = (j == 0) ? y - 1 : j - 1;

                total_energy += -J * L(i, j) * (L(left, j) + L(i, up)) - B * L(i, j);

            }
        }
        return total_energy;
    };

    template<typename E, typename M>
    double dE( const Lattice<E, M> &L, std::array<int,3> xyz) const {
        auto shape = L.shape();
        int size_x = shape[0];
        int size_y = shape[1];

        int i = xyz[0];
        int j = xyz[1];

        int left = (i == 0) ? size_x - 1 : i - 1;
        int right = (i == size_x - 1) ? 0 : i + 1;

        int up = (j == 0) ? size_y - 1 : j - 1;
        int down = (j == size_y - 1) ? 0 : j + 1;

        double dlink_energy = 2*J*L(i, j) * ( L(left, j) + L(right, j) + L(i, up) + L(i, down) );
        double dmag_energy = 2*B*L(i, j);

        return  dlink_energy + dmag_energy;
    };

private:
    double B;
    double J;
};

#endif //ENERGY_H