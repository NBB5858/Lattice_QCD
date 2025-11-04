#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <random>

class IsingLattice1D {
public:
    explicit IsingLattice1D(int lattice_size, double beta, double J, double B);
    int size() const {return lattice_size;};
    double beta() const {return lattice_beta;};
    double J() const {return lattice_J;};
    double B() const {return lattice_B;};
    double energy() const {return lattice_energy;};
    double mag() const {return lattice_mag;};
    void metro_iter();
    void print();

private:
    int lattice_size;
    double lattice_beta;
    double lattice_J;
    double lattice_B;
    double lattice_energy;
    double lattice_mag;
    std::vector<int> lattice_grid;
    std::mt19937 gen;
    std::uniform_int_distribution<int> site_dist;
    std::uniform_int_distribution<int> spin_dist;
    std::uniform_real_distribution<> uniform_dist;

    int rand_site();
    int rand_spin();
    double compute_dE(int site, int spin) const;
    void update_mag(int site, int spin);
};

#endif //LATTICE_H
