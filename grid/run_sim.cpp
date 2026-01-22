
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <array>

#include "Grid.h"
#include "Fields/ScalarField.h"
#include "Lattice.h"
#include "Expression.h"
#include "openmp/Reductions.h"
#include "Action.h"
#include "Observables/Magnetization.h"
#include "HMC.h"
#include "RNG.h"


int GridThread::_threads = 1;

struct Params {
    float mass2;
    float lambda;
    int nsteps;
    float epsilon;
    float Lsteps;
    float sigma;
};

std::vector<Params> read_params(const std::string& filename) {

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

        ss >> p.mass2 >> comma >> p.lambda >> comma >> p.nsteps >> comma >> p.epsilon >> comma >> p.Lsteps >> comma >> p.sigma;
        params.push_back(p);
    }
    return params;
}


int main() {

    /***************
     * File Setup
     ***************/
    std::string input_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/temp_data/phi4/input_params.csv";
    std::vector<Params> params = read_params(input_filename);

    std::string output_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/Results/phi4/phases.csv";
    std::ofstream fout(output_filename);

    if (!fout)
        throw std::runtime_error("Error: Could not open results file: " + output_filename);

    fout << "m^2,lambda,nsteps,mag,abs mag\n";

    /***************
     * Runner
     ***************/

    GridThread::SetMaxThreads();
    std::cout << "Threads: " << GridThread::GetThreads() << std::endl;

    constexpr int dimension = 4;
    constexpr std::array<int, dimension> dims({4, 4, 4, 4});


    GridBase<dimension> Grid(dims);

    for( Params P : params) {
        float mass2 = P.mass2;
        float lambda = P.lambda;
        int nsteps = P.nsteps;
        float epsilon = P.epsilon;
        float Lsteps = P.Lsteps;
        float sigma = P.sigma;


        std::cout << "Running for "
                  << "nsteps = " << nsteps << ", "
                  << "mass2 = " << mass2 << ", "
                  << "lambda = " << lambda <<", "
                  << "epsilon = " << epsilon <<", "
                  << "Lsteps = " << Lsteps << ", "
                  << "sigma = " << sigma <<", "
                  << std::endl;

        std::random_device rd;
        unsigned int base_seed = rd();

        const RandNormal Prng(0,sigma, 42); // copied in HMC, will be reseeded in each thread
        const RandNormal Phirng(0,1, 42); // copied in HMC, will be reseeded in each thread

        const ScalarAction action(mass2, lambda);
        Magnetization mag(nsteps);

        HMC<dimension, ScalarField, ScalarAction, RandNormal, Magnetization> hmc(nsteps, action, Prng, Phirng, epsilon, Lsteps,
                                                                      base_seed, &Grid, mag);

        hmc.Run();
        std::cout << "Acceptance Ratio: " << hmc.log().acceptance_ratio() << std::endl;

        std::cout << "mean " << mag.mean(0.1) << std::endl;

        fout << mass2 << ","
              << lambda << ","
              << nsteps << ","
              << mag.mean(0.1) << ","
              << mag.abs_mean(0.1) << ","
              << "\n";
    }

    std::cout << "Done." << std::endl;

}

