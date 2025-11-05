#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <random>
#include <stdexcept>
#include <iostream>

template <typename EnergyModel, typename MagModel>
class Lattice {
public:
    Lattice(int lattice_size, double beta, EnergyModel energy_model, MagModel mag_model);

    int size() const {return lattice_size;};
    int operator[] (int i) const {return lattice_grid[i];};

    double beta() const {return lattice_beta;};
    double energy() const {return lattice_energy;};
    double mag() const {return lattice_mag;};

    void metro_iter();
    void print() const;

private:
    int lattice_size;
    double lattice_beta;
    double lattice_energy;
    double lattice_mag;

    EnergyModel energy_model;
    MagModel mag_model;

    std::vector<int> lattice_grid;
    std::mt19937 gen {std::random_device{}()};
    std::uniform_int_distribution<int> site_dist{0, size-1};
    std::uniform_int_distribution<int> spin_dist{0, 1};
    std::uniform_real_distribution<> uniform_dist{0, 1.0};

    int rand_site() {return site_dist(gen);};
    int rand_spin() {return spin_dist(gen) == 0 ? -1 : 1;};

    double compute_dE(int site, int spin) const { return energy_model.dE(*this, site, spin); };
    double compute_dM(int site, int spin) const { return mag_model.dM(*this, site, spin); };
};

template <typename E, typename M>
Lattice<E, M>::Lattice(int size, double beta, E energy_model, M mag_model)
: lattice_size(size),
  lattice_beta(beta),
  energy_model(energy_model),
  mag_model(mag_model),
  lattice_grid(size),
  site_dist(0, size-1)
{
    if (size <= 1)
        throw std::invalid_argument("Lattice size must be > 1");

    for(int i = 0; i < size; ++i)
        lattice_grid[i] = rand_spin();

    lattice_energy = energy_model.total_energy(*this);
    lattice_mag = mag_model.total_mag(*this);
}


template <typename E, typename M>
void Lattice<E, M>::metro_iter() {

    int site = rand_site();
    int spin = rand_spin();

    double dE = compute_dE(site, spin);
    double accep_prob = std::exp(-dE * lattice_beta);

    double u = uniform_dist(gen);
    if (u <= accep_prob) {
        double dM = compute_dM(site, spin);

        lattice_mag += dM;
        lattice_energy += dE;

        lattice_grid[site] = spin;
    }
}

template <typename E, typename M>
void Lattice<E, M>::print() const {
    for(int i = 0; i < lattice_size; ++i) {
        if( lattice_grid[i] == -1 )
            std::cout << "-";
        else if ( lattice_grid[i] == 1 ) {
            std::cout << "+";
        }
    }
    std::cout << std::endl;
}

#endif //LATTICE_H
