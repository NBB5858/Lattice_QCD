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

    std::vector<double> results(params.size());

    std::transform(
        std::execution::par,
        params.begin(),
        params.end(),
        results.begin(),
        run_simulation
        );

    write_results_csv("../data/results.csv", params, results);
}



double run_simulation(const Params &P) {

    IsingE2D energy_model(P.J, P.B);
    IsingM2D mag_model;

    std::array<int,2> dims = {P.X, P.Y};

    Lattice<IsingE2D, IsingM2D, 2> lattice(dims, P.beta, energy_model, mag_model);

    double mag = lattice.mag();
    double mag_sum = mag;
    for( long long _ = 1; _ < P.iter; _++) {
        lattice.metro_iter();

        mag = lattice.mag();
        mag_sum += mag;
    }

    return mag_sum / P.iter;
}

