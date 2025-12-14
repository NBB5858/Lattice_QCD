#ifndef MAGNETIZATION_H
#define MAGNETIZATION_H

template <typename EnergyModel, typename MagModel, int D>
class Lattice;


struct IsingM1D {
    IsingM1D() = default;

    template<typename E, typename M, int D>
    void bind_geometry( const Lattice<E, M, D> &L ) { size_x = L.shape()[0]; }

    template<typename E, typename M, int D>
    double total_mag( const Lattice<E, M, D> &L ) const {
        double total_mag = 0; for(int i = 0; i < size_x; ++i) total_mag += L(i);
        return total_mag / static_cast<double>(size_x);
    };

    template<typename E, typename M, int D>
    double dM( const Lattice<E, M, D> &L, std::array<int,D> site) const {
        int x = site[0];
        return -2*L(x) / static_cast<double>(size_x);
    };

    int size_x;
};


struct IsingM2D {
    IsingM2D() = default;

    template<typename E, typename M, int D>
    void bind_geometry( const Lattice<E, M, D> &L ) {
        auto [nx, ny] = L.shape(); size_x=nx, size_y=ny; size=size_x*size_y;
    }

    template<typename E, typename M, int D>
    double total_mag( const Lattice<E, M, D> &L ) const {
        double total_mag = 0;
        for(int i = 0; i < size_x; ++i) {
            for(int j = 0; j < size_y; ++j) {
                total_mag += L(i, j);
            }
        }
        return total_mag / size;
    };

    template<typename E, typename M, int D>
    double dM( const Lattice<E, M, D> &L, std::array<int,D> xy) const {
        auto [x, y] = xy;
        return -2*L(x, y) / size;
    };

    int size_x, size_y;
    double size;

};

#endif //MAGNETIZATION_H
