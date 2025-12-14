
#include <iostream>
#include <vector>

#include "Grid.h"
#include "Lattice.h"
#include "ScalarField.h"
#include "VectorField.h"
#include "Expression.h"
#include "Reductions.h"
#include "Action.h"
#include "Observables.h"
#include "HMC.H"


int GridThread::_threads = 1;

int main() {

    GridThread::SetMaxThreads();

    std::cout << "Threads: " << GridThread::GetThreads() << std::endl;

    const std::vector<int> dims({3, 3});
    GridBase Grid(dims);

    constexpr int nsteps = 100;
    const ScalarAction action(1.0, 1.0);
    const RandNormal rng(0,1);
    const Magnetization mag(nsteps);

    HMC<ScalarField, ScalarAction, RandNormal, Magnetization> hmc(nsteps, action, rng, &Grid, mag);

    hmc.Run();

    std::cout << mag.mean() << std::endl;
}

