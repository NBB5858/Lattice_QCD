#include <iostream>
#include "lattice.h"

IsingLattice1D::IsingLattice1D(int size, double beta, double J, double B)
: lattice_size(size),
  lattice_beta(beta),
  lattice_J(J),
  lattice_B(B),
  lattice_energy(0.0),
  lattice_mag(0.0),
  lattice_grid(size),
  gen(std::random_device{}()),
  site_dist(0, size-1),
  spin_dist(0, 1),
  uniform_dist(0, 1)
{
    if (size <= 1)
        throw std::invalid_argument("Lattice size must be > 1");

    for(int i = 0; i < size; ++i)
        lattice_grid[i] = rand_spin();

    for(int i = 0; i < lattice_size; ++i) {
        int left = (i == 0) ? lattice_size - 1 : i - 1;
        lattice_energy += -lattice_grid[i] * lattice_grid[left] * lattice_J - lattice_grid[i] * B;
        lattice_mag += lattice_grid[i] / static_cast<double>(lattice_size);
    }
}

int IsingLattice1D::rand_site() {
    return site_dist(gen);
}

int IsingLattice1D::rand_spin() {
    return spin_dist(gen) == 0 ? -1 : 1;
}

double IsingLattice1D::compute_dE(int site, int spin) const {

    int left = (site == 0) ? lattice_size - 1 : site - 1;
    int right = (site == lattice_size-1) ? 0 : site + 1;

    double old_link_energy = -lattice_J * lattice_grid[site] * ( lattice_grid[left] + lattice_grid[right] );
    double old_mag_energy = -lattice_B * lattice_grid[site];

    double new_link_energy = -lattice_J * spin * ( lattice_grid[left] + lattice_grid[right] );
    double new_mag_energy = -lattice_B * spin;

    return  new_link_energy + new_mag_energy - old_link_energy - old_mag_energy;
}

void IsingLattice1D::update_mag(int site, int spin) {
    lattice_mag += (spin - lattice_grid[site]) / static_cast<double>(lattice_size);
}

void IsingLattice1D::metro_iter() {

    int site = rand_site();
    int spin = rand_spin();

    double dE = compute_dE(site, spin);
    double accep_prob = std::exp(-dE * lattice_beta);

    double u = uniform_dist(gen);
    if (u <= accep_prob) {
        update_mag(site, spin);
        lattice_energy += dE;

        lattice_grid[site] = spin;
    }
}

void IsingLattice1D::print() {
    for(int i = 0; i < lattice_size; ++i) {
        if( lattice_grid[i] == -1 )
            std::cout << "-";
        else if ( lattice_grid[i] == 1 ) {
            std::cout << "+";
        }
    }
    std::cout << std::endl;
}
