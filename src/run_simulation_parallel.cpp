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

    IsingE1D energy_model(P.J, P.B);
    IsingM1D mag_model;

    Lattice lattice(P.N, P.beta, energy_model, mag_model);

    double mag = lattice.mag();

    double mag_sum = mag;
    for( long long _ = 1; _ < P.iter; _++) {
        lattice.metro_iter();

        mag = lattice.mag();
        mag_sum += mag;
    }

    return mag_sum / P.iter;
}

