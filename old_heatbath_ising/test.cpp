#include "lattice.h"
#include "energy.h"
#include "magnetization.h"
#include "io.h"
#include <execution>
#include <algorithm>

double run_simulation(const Params &P);

int main() {

    std::string filename = "../data/params.csv";
    std::vector<Params> params = read_params(filename);

    Params P = params[0];

    std::array<int,3> dims = {P.X, P.Y, P.Z};

    IsingE2D energy_model(P.J, P.B);
    IsingM2D mag_model;

    Lattice lattice(dims, P.beta, energy_model, mag_model);

    //std::cout << lattice.size();

    lattice.print();
}


