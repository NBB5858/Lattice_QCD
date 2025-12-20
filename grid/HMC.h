#ifndef HMC_H
#define HMC_H

#include<cmath>
#include<tuple>
#include <chrono>


#include "Action.h"
#include "RNG.h"
#include "Reductions.h"
#include "Observables.h"

struct Logger {
    double acceptance_ratio() const {
        double output = static_cast<double>(accepts) / static_cast<double>(nsteps);
        return output;
    }

    int accepts = 0;
    int nsteps = 0;
};

// To Do: Break this up
template<typename FieldType, typename ActionType, typename RNGType, typename... Obs >
class HMC {
public:

    HMC(int nsteps, ActionType action, RNGType Prng, RNGType Phirng, double epsilon, int Lsteps, unsigned int base_seed, GridBase* grid, Obs&... obs) :
        _nsteps(nsteps),
        _action(action),
        _grid(grid),
        _phi(grid),
        _epsilon(epsilon),
        _Lsteps(Lsteps),
        _base_seed(base_seed),
        _urng(0, 1, base_seed),
        _obs(obs...)
    {
        int nthreads = thread_max(0);
        _Prngs.assign(nthreads, Prng);
        _Phirngs.assign(nthreads, Phirng);

        for (int t = 0; t < nthreads; ++t) {
            _Prngs[t].reseed(base_seed + 2*t);
            _Phirngs[t].reseed(base_seed + 2*t + 1);
        }
        _log.nsteps = _nsteps;
    }

    const Lattice<FieldType>& Lat() const {return _phi;}

    void RandomizeLattice( Lattice<ScalarField>& L, std::vector<RNGType>& rngs ) {
        thread_for(ss, _grid->osites(), {
            int thread = thread_num(ss);
            L.Randomize(ss, rngs[thread]);
        });
    }

    void hmc_step() {
        Lattice<FieldType> P(_grid);
        RandomizeLattice(P, _Prngs);

        // make a copy constructor
        Lattice<FieldType> phi_tmp(_grid);
        phi_tmp = _phi;

        double initial = sum(0.5*P*P, _grid) + _action.S(_phi);

        P   -= 0.5 * _action.Force(_phi) * _epsilon;
        for(int i=0; i<_Lsteps-1; ++i) {
            phi_tmp += P * _epsilon;
            P   -= _action.Force(phi_tmp) * _epsilon;
        }
        phi_tmp += P * _epsilon;
        P   -= 0.5 * _action.Force(phi_tmp) * _epsilon;

        double final = sum(0.5*P*P, _grid) + _action.S(phi_tmp);

        double prob_accept = std::min(1.0, std::exp((initial - final)));
        double metro_draw = _urng.draw();
        if(metro_draw < prob_accept) {
            _phi = phi_tmp;
            _log.accepts += 1;
        };
    }

    void measure() {
        std::apply([&](auto& ob){ob.measure(_phi);}, _obs);
    }

    void Run() {
        RandomizeLattice(_phi, _Phirngs);

        for(int i = 0; i < _nsteps; ++i) {
            hmc_step();
            measure();
        }
    }

    Logger& log() { return _log; };

private:
    int _nsteps;
    ActionType _action;
    GridBase* _grid;
    Lattice<FieldType> _phi;
    double _epsilon;
    int _Lsteps;
    unsigned int _base_seed;
    std::tuple<Obs&...> _obs;
    std::vector<RNGType> _Prngs;
    std::vector<RNGType> _Phirngs;
    RandUniform _urng;
    Logger _log;
};

#endif //HMC_H
