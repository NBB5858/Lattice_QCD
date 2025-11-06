#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <random>
#include <stdexcept>
#include <iostream>
#include <array>

template <typename EnergyModel, typename MagModel>
class Lattice {
public:
    Lattice( std::array<int,3> dims, double beta, EnergyModel energy_model, MagModel mag_model );

    int size() const {return size_;};
    std::array<int, 3> shape() const {return {dims_[0], dims_[1], dims_[2]};};

    int& operator()(int i, int j, int k) { return lattice_grid[i*strides_[0] + j*strides_[1] + k*strides_[2]];};
    int& operator()(int i, int j) { return lattice_grid[i*strides_[0] + j*strides_[1]];};
    int& operator()(int i) { return lattice_grid[i*strides_[0]];};

    int operator()(int i, int j, int k) const { return lattice_grid[i*strides_[0] + j*strides_[1] + k*strides_[2]];};
    int operator()(int i, int j) const { return lattice_grid[i*strides_[0] + j*strides_[1]];};
    int operator()(int i) const { return lattice_grid[i*strides_[0]];};

    double beta() const {return beta_;};
    double energy() const {return energy_;};
    double mag() const {return mag_;};

    void metro_iter();
    void print() const;

private:
    int size_;
    std::array<int,3> dims_;
    std::array<int,3> strides_;

    double beta_;
    double energy_;
    double mag_;

    EnergyModel energy_model_;
    MagModel mag_model_;

    std::vector<int> lattice_grid;
    std::mt19937 gen;
    std::uniform_int_distribution<int> i_dist, j_dist, k_dist;
    std::uniform_real_distribution<> uniform_dist;

    std::array<int,3> rand_site() {return {i_dist(gen), j_dist(gen), k_dist(gen)};};

    double compute_dE(std::array<int, 3> xyz) const { return energy_model_.dE(*this, xyz); };
    double compute_dM(std::array<int, 3> xyz) const { return mag_model_.dM(*this, xyz); };
};

template <typename E, typename M>
Lattice<E, M>::Lattice(std::array<int,3> dims, double beta, E energy_model, M mag_model)
: size_(dims[0]*dims[1]*dims[2]),
  dims_(dims),
  strides_{dims[1]*dims[2], dims[2], 1},
  beta_(beta),
  energy_model_(energy_model),
  mag_model_(mag_model),
  lattice_grid(size_),
  gen(std::random_device{}()),
  i_dist(0, dims_[0]-1),
  j_dist(0, dims_[1]-1),
  k_dist(0, dims_[2]-1),
  uniform_dist(0.0, 1.0)
{
    if (size_ <= 1)
        throw std::invalid_argument("Lattice size must be > 1");
    if (dims_[0] <= 0 || dims_[1] <= 0 || dims_[2] <= 0)
        throw std::invalid_argument("All lattice dimensions must be > 0");
    if (dims_[0] == 1 || (dims_[1] == 1 and dims_[2] != 1) )
        throw std::invalid_argument("Dims of size 1 must be trailing");


    for(int i = 0; i < size_; ++i)
        lattice_grid[i] = 1;

    energy_ = energy_model_.total_energy(*this);
    mag_ = mag_model_.total_mag(*this);
}


template <typename E, typename M>
void Lattice<E, M>::metro_iter() {

    std::array<int, 3> xyz = rand_site();
    auto [i,j,k] = xyz;

    double dE = compute_dE(xyz);
    double accep_prob = std::exp(-dE * beta_);

    double u = uniform_dist(gen);
    if (u <= accep_prob) {
        double dM = compute_dM(xyz);
        mag_ += dM;
        energy_ += dE;
        (*this)(i, j, k) = -(*this)(i, j, k);
    }
}

template <typename E, typename M>
void Lattice<E, M>::print() const {
    for (int i = 0; i < dims_[0]; ++i) {
        for (int j = 0; j < dims_[1]; ++j) {
            int s = (*this)(i,j);
            std::cout << (s == -1 ? '-' : '+');
        }
        std::cout << '\n';
    }
}

#endif //LATTICE_H
