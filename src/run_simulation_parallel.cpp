#include "lattice.h"
#include <fstream>
#include <sstream>
#include <string>
#include <execution>
#include <algorithm>

struct Params {
    int N;
    double beta;
    double B;
    double J;
    long long iter;
};

std::vector<Params> read_params(std::string& filename);
double run_simulation(Params params);
void write_results_csv(const std::string& filename,
                       const std::vector<Params>& params,
                       const std::vector<double>& results);


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

std::vector<Params> read_params(std::string& filename) {

    std::ifstream fin(filename);
    if (!fin)
        throw std::runtime_error("Error: could not open " + filename);

    std::vector<Params> params;

    std::string line;
    std::getline(fin, line);

    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        Params p;
        char comma;

        ss >> p.N >> comma >> p.beta >> comma >> p.B >> comma >> p.J >> comma >> p.iter;
        params.push_back(p);
    }
    return params;
}

double run_simulation(Params P) {

    IsingLattice1D lattice(P.N, P.beta, P.J, P.B);
    double mag = lattice.mag();

    double mag_sum = mag;
    for( long long _ = 1; _ < P.iter; _++) {
        lattice.metro_iter();

        mag = lattice.mag();
        mag_sum += mag;
    }

    return mag_sum / P.iter;
}

void write_results_csv(const std::string& filename,
                       const std::vector<Params>& params,
                       const std::vector<double>& results)
{
    std::ofstream fout(filename);
    if (!fout)
        throw std::runtime_error("Error: Could not open results file: " + filename);

    // header
    fout << "N,beta,B,J,iter,mag\n";

    // rows
    for (size_t i = 0; i < params.size(); ++i) {
        const Params& p = params[i];
        fout << p.N      << ","
             << p.beta   << ","
             << p.B      << ","
             << p.J      << ","
             << p.iter   << ","
             << results[i] << "\n";
    }
}