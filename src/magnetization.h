#ifndef MAGNETIZATION_H
#define MAGNETIZATION_H

template <typename EnergyModel, typename MagModel>
class Lattice;


class IsingM1D {
public:
    IsingM1D() = default;

    template<typename E, typename M>
    double total_mag( const Lattice<E, M> &L ) const {
        double total_mag = 0;
        for(int i = 0; i < L.shape()[0]; ++i) {
            total_mag += L(i);
        }
        return total_mag / static_cast<double>(L.size());
    };

    template<typename E, typename M>
    double dM( const Lattice<E, M> &L, std::array<int,3> xyz) const {
        int x = xyz[0];
        return -2*L(x) / static_cast<double>(L.size());
    };
};


class IsingM2D {
public:
    IsingM2D() = default;

    template<typename E, typename M>
    double total_mag( const Lattice<E, M> &L ) const {
        auto shape = L.shape();
        int x = shape[0];
        int y = shape[1];

        double total_mag = 0;
        for(int i = 0; i < x; ++i) {
            for(int j = 0; j < y; ++j) {
                total_mag += L(i, j);
            }
        }
        return total_mag / static_cast<double>(L.size());
    };

    template<typename E, typename M>
    double dM( const Lattice<E, M> &L, std::array<int,3> xyz) const {
        int x = xyz[0];
        int y = xyz[1];
        return -2*L(x, y) / static_cast<double>(L.size());
    };
};

#endif //MAGNETIZATION_H
