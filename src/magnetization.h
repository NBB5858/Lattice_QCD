#ifndef MAGNETIZATION_H
#define MAGNETIZATION_H

template <typename EnergyModel, typename MagModel>
class Lattice;

class IsingM1D {
public:
    IsingM1D() = default;

    template<typename E, typename M>
    double total_mag( const Lattice<E, M> &L ) const;

    template<typename E, typename M>
    double dM( const Lattice<E, M> &L, int site, int spin) const;
};


template <typename E, typename M>
double IsingM1D::total_mag( const Lattice<E, M> &L )  const {

    double total_mag = 0;
    for(int i = 0; i < L.size(); ++i) {
        total_mag += L[i] / static_cast<double>(L.size());
    }

    return total_mag;
}

template <typename E, typename M>
double IsingM1D::dM( const Lattice<E, M> &L, int site, int spin ) const {
    return (spin - L[site]) / static_cast<double>(L.size());
}

#endif //MAGNETIZATION_H
