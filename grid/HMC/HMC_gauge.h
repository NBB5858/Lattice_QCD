#ifndef HMC_GAUGE_H
#define HMC_GAUGE_H


#include <vector>
#include <complex>
#include <algorithm>

#include "../openmp/pragma_shorts.h"
#include "../openmp/Reductions.h"
#include "../RNG/RNG.h"


struct Logger {
    double acceptance_ratio() const {
        double output = static_cast<double>(accepts) / static_cast<double>(nsteps);
        return output;
    }

    int accepts = 0;
    int nsteps = 0;
};




template<int d, typename FieldType, typename IntegratorType, typename PhiRNGType, typename PRNGType, typename... Obs >
class HMCGauge {

    using GaugeField = FieldType;
    using MomField   = typename FieldTraits<GaugeField>::MomentumField;

public:

    HMCGauge(int nsteps, IntegratorType integrator, PhiRNGType Phirng, PRNGType Prng, unsigned int base_seed, GridBase<d>* grid, Obs&... obs) :
        _nsteps(nsteps),
        _integrator(integrator),
        _grid(grid),
        _phi(grid),
        _base_seed(base_seed),
        _urng(0, 1, base_seed),
        _obs(obs...)
    {
        int nthreads = thread_max(0);
        _prngs.assign(nthreads, Prng);
        _phirngs.assign(nthreads, Phirng);

        for (int t = 0; t < nthreads; ++t) {
            _prngs[t].reseed(base_seed + 2*t);
            _phirngs[t].reseed(base_seed + 2*t + 1);
        }
        _log.nsteps = _nsteps;
    }

    const Lattice<FieldType, d>& Lat() const {return _phi;}

    void hmc_step() {
        Lattice<MomField, d> P(_grid);
        P.RandomizeLattice(_prngs);

        Lattice<FieldType, d> phi_tmp(_phi.grid());
        phi_tmp = _phi;

        std::complex<float> initial = _integrator.H(phi_tmp, P);
        _integrator.evolve(phi_tmp, P);
        std::complex<float> final = _integrator.H(phi_tmp, P);

        float delta = (final - initial).real();

        double prob_accept = std::min(1.0f, std::exp(-delta));
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
        _phi.RandomizeLattice(_phirngs);

        for(int i = 0; i < _nsteps; ++i) {
            hmc_step();
            measure();
        }
    }

    Logger& log() { return _log; };

private:
    int _nsteps;
    IntegratorType _integrator;
    GridBase<d>* _grid;
    Lattice<FieldType, d> _phi;
    unsigned int _base_seed;
    RandUniform _urng;
    std::tuple<Obs&...> _obs;
    std::vector<PhiRNGType> _phirngs;
    std::vector<PRNGType> _prngs;
    Logger _log;
};



#endif //HMC_GAUGE_H
