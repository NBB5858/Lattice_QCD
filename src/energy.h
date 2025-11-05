#ifndef ENERGY_H
#define ENERGY_H

template <typename EnergyModel, typename MagModel>
class Lattice;


class IsingE1D {
public:
    IsingE1D( double J, double B ) : B(B), J(J) {};

    template<typename E, typename M>
    double total_energy( const Lattice<E, M> &L ) const;

    template<typename E, typename M>
    double dE( const Lattice<E, M> &L, int site, int spin) const;

private:
    double B;
    double J;
};


template <typename E, typename M>
double IsingE1D::total_energy( const Lattice<E, M> &L )  const {

    double total_energy = 0;
    for(int i = 0; i < L.size(); ++i) {
        int left = (i == 0) ? L.size() - 1 : i - 1;
        total_energy += -L[i] * L[left] * J - L[i] * B;
    }

    return total_energy;
}

template <typename E, typename M>
double IsingE1D::dE( const Lattice<E, M> &L, int site, int spin ) const {

    int left = (site == 0) ? L.size() - 1 : site - 1;
    int right = (site == L.size() - 1) ? 0 : site + 1;

    double old_link_energy = -J * L[site] * ( L[left] + L[right] );
    double old_mag_energy = -B * L[site];

    double new_link_energy = -J * spin * ( L[left] + L[right] );
    double new_mag_energy = -B * spin;

    return  new_link_energy + new_mag_energy - old_link_energy - old_mag_energy;
}
#endif //ENERGY_H