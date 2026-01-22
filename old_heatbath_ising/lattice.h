#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <random>
#include <iostream>
#include <array>

template <typename EnergyModel, typename MagModel, int D>
class Lattice {
public:
    Lattice( std::array<int,D>& dims, double beta, EnergyModel energy_model, MagModel mag_model);

    int size() const {return size_;}
    std::array<int,D> shape() const { return dims_;}

    template <typename... Ints>
    const int& operator()(Ints... idxs) const {
        static_assert(sizeof...(idxs) == D);
        std::array<int, D> x{idxs...};
        return lattice_grid_[lin(x)];
    }

    double beta() const {return beta_;}
    double energy() const {return energy_;}
    double mag() const {return mag_;}

    void metro_iter();

private:
    int size_;
    std::array<int,D> dims_;
    std::array<int,D> strides_;
    std::vector<int> lattice_grid_;

    double beta_;
    double energy_;
    double mag_;

    EnergyModel energy_model_;
    MagModel mag_model_;

    double compute_dE(std::array<int, D> site) const { return energy_model_.template dE<EnergyModel, MagModel, D>(*this, site); }
    double total_E() const { return energy_model_.template total_energy<EnergyModel, MagModel, D>(*this); }

    double compute_dM(std::array<int, D> site) const { return mag_model_.template dM<EnergyModel, MagModel, D>(*this, site); }
    double total_M() const { return mag_model_.template total_mag<EnergyModel, MagModel, D>(*this); }

    std::mt19937 gen_;
    std::uniform_int_distribution<int> axis_dist_[D];
    std::uniform_real_distribution<> uniform_dist_;

    std::array<int,D> rand_site() {
        std::array<int,D> v; for (int a = 0; a < D; ++a) v[a] = axis_dist_[a](gen_);
        return v;
    }

    int lin(const std::array<int, D>& x) const {
        int idx=0; for (int a=0; a<D; ++a) idx += x[a]*strides_[a]; return idx;
    }
};

template <typename E, typename M, int D>
Lattice<E, M, D>::Lattice( std::array<int,D>& dims, double beta, E energy_model, M mag_model)
: dims_(dims),
  beta_(beta),
  energy_model_(energy_model),
  mag_model_(mag_model)
{
    strides_[D-1] = 1; for (int a = D - 2; a >= 0; --a) strides_[a] = strides_[a+1]*dims_[a+1];
    size_ = 1; for (int a = 0; a < D; ++a) size_ = size_*dims_[a];
    lattice_grid_.assign(size_, +1);

    for (int a = 0; a < D; ++a) axis_dist_[a] = std::uniform_int_distribution<int>(0, dims_[a]-1);
    uniform_dist_ = std::uniform_real_distribution<> (0.0, 1.0);
    gen_.seed(std::random_device{}());

    energy_model_.bind_geometry(*this);
    mag_model_.bind_geometry(*this);

    energy_ = total_E();
    mag_ = total_M();
}


template <typename E, typename M, int D>
void Lattice<E, M, D>::metro_iter() {

    std::array<int, D> site = rand_site();

    double dE = compute_dE(site);
    double accep_prob = std::exp(-beta_ * dE);

    double u = uniform_dist_(gen_);
    if (u <= accep_prob) {
        double dM = compute_dM(site);
        mag_ += dM;
        energy_ += dE;
        lattice_grid_[lin(site)] = -lattice_grid_[lin(site)];
    }
}


#endif //LATTICE_H
