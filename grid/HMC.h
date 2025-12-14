#ifndef HMC_H
#define HMC_H

#include<cmath>
#include<tuple>

#define MIN(a, b) a <= b ? a : b

#include "Action.h"
#include "RNG.h"
#include "Reductions.h"
#include "Observables.h"


template<typename FieldType, typename ActionType, typename RNGType, typename... Obs >
class HMC {
public:

    HMC(int nsteps, ActionType action, RNGType rng, GridBase* grid, Obs... obs) :
        _nsteps(nsteps),
        _action(action),
        _grid(grid),
        _phi(grid),
        _obs(std::move(obs)...)
    {
        int base_seed = 42;

        int nthreads = thread_max(0);
        _rngs.reserve(nthreads);

        for(int t=0; t<nthreads; ++t) {
            _rngs.emplace_back(rng);
            _rngs[t].reseed(base_seed + t);
        }
    }

    const Lattice<FieldType>& Lat() const {return _phi;}

    void RandomizeLattice( Lattice<ScalarField>& L ) {
        thread_for(ss, _grid->osites(), {
            int thread = thread_num(ss);
            L.Randomize(ss, _rngs[thread]);
        });
    }

    //debugging
    void RandomizeLattice() {
        thread_for(ss, _grid->osites(), {
            int thread = thread_num(ss);
            _phi.Randomize(ss, _rngs[thread]);
        });
    }

    void hmc_step() {

        Lattice<FieldType> phi_tmp(_grid);
        Lattice<FieldType> P(_grid);
        RandomizeLattice(P);

        double initial = sum(0.5*P*P, _grid) + _action.S(_phi);

        P   -= 0.5 * _action.Force(_phi);
        phi_tmp = _phi + P;
        P   -= 0.5 * _action.Force(phi_tmp);

        double final = sum(0.5*P*P, _grid) + _action.S(phi_tmp);

        std::cout << "initial: " << initial << "," << "final: " << final << std::endl;

        double prob_accept = MIN(1, std::exp(initial - final));
        double metro_draw = _urng.draw();
        if(prob_accept < metro_draw) _phi = phi_tmp;
    }

    void measure() {
        std::apply([&](auto& ob){ob.measure(_phi);}, _obs);
    }

    void Run() {
        RandomizeLattice();
        // _phi.Print2D();
        // std::cout << std::endl;

        for(int i = 0; i < _nsteps; ++i) {
            hmc_step();
            //measure();
        }
    }

private:
    int _nsteps;
    ActionType _action;
    GridBase* _grid;
    Lattice<FieldType> _phi;
    std::tuple<Obs...> _obs;
    std::vector<RNGType> _rngs;
    RandUniform _urng;
};

#endif //HMC_H
