
#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <vector>
#include <array>
#include <cassert>

#include "../Grid.h"
#include "Expression.h"
#include "../simd/SimdFieldTypes.h"
#include "../Shifting/Shifting.h"
#include "../openmp/pragma_shorts.h"
#include "../openmp/Reductions.h"

class LatticeBase {};

template<typename FieldType, int d>
class VectorField;


template<class FieldType, int d>
class Lattice : public LatticeBase {

    static constexpr int W = simd::W;

    using Traits = FieldTraits<FieldType>;
    using Storage = typename Traits::Storage;
    using View = typename Traits::View;
    using ConstView =  typename Traits::ConstView;

    int _nBlocks;   // number of blocks (of W sites)
    int _nSites;    // total number of sites

    GridBase<d>* _grid;
    std::vector<Storage> _mem{};

public:
    explicit Lattice(GridBase<d>* grid) :
    _grid(grid)
    {
        _nBlocks = grid->_osites;
        _nSites  = _nBlocks * W;
        _mem.resize(static_cast<size_t>(_nBlocks)); // remember: resize calls the default constructor
    }

    GridBase<d>* grid() const { return _grid; }
    int blocks() const { return _nBlocks; }
    int sites()  const { return _nSites; }

    // access underlying block storage
    Storage& operator()(uint64_t ss) { return _mem[ss]; }
    const Storage& operator()(uint64_t ss) const { return _mem[ss]; }

    View site(int s) {
        int b = s / W;
        int l = s % W;
        return Traits::view(_mem[b], l);
    }

    ConstView site(int s) const {
        int b = s / W;
        int l = s % W;
        return Traits::view(_mem[b], l);
    }

    void print2d() const {
        for(int i=0; i<_grid->dims()[0]; ++i) {
            for(int j=0; j<_grid->dims()[1]; ++j) {
                int ss = i + j * _grid->dims()[0];
                ConstView sv = site(ss);
                sv.print();
                std::cout << " ";
            }
            std::cout << std::endl;
        }
    }

    void set(int s, float v) {
        assert(s < _nSites && "Lattice::set: site index out of range");
        View sv = site(s);
        sv.set(v);
    }

    template<class RNGType>
    void RandomizeLattice( std::vector<RNGType>& rngs ) {
        thread_for (ss, _grid->osites(), {
            int thread = thread_num(ss);
            FieldTraits<FieldType>::Randomize(_mem[ss], rngs[thread]);
        });
    }

    template <typename T1,
              typename = std::enable_if_t<is_lattice_expr<T1>::value>>
    inline Lattice& operator=(const T1 &expr) {
        thread_for(ss, _grid->_osites, {
            _mem[ss] = eval(ss, expr);
        });
        return *this;
    };

    template <typename T1,
              std::enable_if_t<is_lattice<T1>::value || is_lattice_expr<T1>::value, int> = 0>
    inline Lattice& operator+=(const T1& expr) {
        return *this = (*this + expr);
    }

    template <typename T1,
              std::enable_if_t<is_lattice<T1>::value || is_lattice_expr<T1>::value, int> = 0>
    inline Lattice& operator-=(const T1& expr) {
        return *this = (*this - expr);
    }

    friend Storage grad_squared(int ss, const Lattice& lat) {
        Storage ret{};
        for (int mu = 0; mu < d; ++mu) {
            std::array<int, d> disp{};
            disp.fill(0);
            disp[mu] = +1;

            ret += square(Cshift(lat, ss, disp) - lat._mem[ss]);
        }
        return ret;
    }

    friend Storage box(int ss, const Lattice& lat) {
        Storage ret{};
        for (int mu = 0; mu < d; ++mu) {
            std::array<int, d> dp{};
            std::array<int, d> dm{};
            dp.fill(0);
            dm.fill(0);
            dp[mu] = +1;
            dm[mu] = -1;

            ret += (Cshift(lat, ss, dp) - 2.0 * lat._mem[ss] + Cshift(lat, ss, dm)) / 1.0;
        }
        return ret;
    }

    friend simd::vReal PlaqReTraceSum(int ss, const Lattice& U ) {
        simd::vReal trace_sum{};

        const auto ssu = static_cast<std::uint64_t>(ss);
        const auto Ux  = U._mem[ss];

        // Precompute x + e_mu for all mu (d gathers total)
        std::array<Storage, static_cast<std::size_t>(d)> Up{};
        for (int mu = 0; mu < d; ++mu) {
            Up[mu] = Cshift(U, ssu, unit_disp<d>(mu, +1));
        }

        for (int mu = 0; mu < d; ++mu) {
            for (int nu = mu + 1; nu < d; ++nu) {

                // U_mu(x)
                const auto& U_mu_x = Ux[mu];

                // U_nu(x+mu)
                const auto& U_nu_xpmu = Up[mu][nu];

                // U_mu(x+nu)
                const auto& U_mu_xpnu = Up[nu][mu];

                // U_nu(x)
                const auto& U_nu_x = Ux[nu];

                auto plaq = U_mu_x
                          * U_nu_xpmu
                          * adj(U_mu_xpnu)
                          * adj(U_nu_x);

                trace_sum += trace(plaq).re;
            }
        }
        return trace_sum;
    }

    friend Storage StapleSum(int ss, const Lattice& U) {
        Storage staple_sum{};

        const auto ssu = static_cast<std::uint64_t>(ss);

        // Cache the local block once (saves repeated loads / aliasing)
        const auto Ux = U._mem[ss];

        // 1) Precompute unit shifts: x + e_a and x - e_a for all directions a
        std::array<Storage, static_cast<std::size_t>(d)> Up{};
        std::array<Storage, static_cast<std::size_t>(d)> Um{};

        for (int a = 0; a < d; ++a) {
            Up[a] = Cshift(U, ssu, unit_disp<d>(a, +1));
            Um[a] = Cshift(U, ssu, unit_disp<d>(a, -1));
        }

        // 2) Precompute mixed shifts: x - e_nu + e_mu for all mu != nu
        //    (You only need [nu] out of this Storage later, but simplest is store full Storage.)
        std::array<std::array<Storage, static_cast<std::size_t>(d)>, static_cast<std::size_t>(d)> Umnu_pm{};

        for (int mu = 0; mu < d; ++mu) {
            for (int nu = 0; nu < d; ++nu) {
                if (nu == mu) continue;
                Umnu_pm[mu][nu] = Cshift(U, ssu, disp2<d>(nu, -1, mu, +1));
            }
        }

        // 3) Accumulate staples using cached shifts
        for (int mu = 0; mu < d; ++mu) {
            for (int nu = 0; nu < d; ++nu) {
                if (nu == mu) continue;

                // Forward staple: U_nu(x) U_mu(x+nu) U_nu^\dagger(x+mu)
                staple_sum[mu] += Ux[nu]
                                * Up[nu][mu]
                                * adj(Up[mu][nu]);

                // Backward staple: U_nu^\dagger(x-nu) U_mu(x-nu) U_nu(x-nu+mu)
                staple_sum[mu] += adj(Um[nu][nu])
                                * Um[nu][mu]
                                * Umnu_pm[mu][nu][nu];
            }
        }

        return staple_sum;
    }


    friend Storage exp(int ss, const Lattice& P) {
        return FieldTraits<FieldType>::exp(P._mem[ss]);
    }

};

#endif //LATTICE_H