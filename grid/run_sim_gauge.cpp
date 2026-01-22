#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <array>
#include <chrono>


#include "Grid.h"
#include "Fields/U1Field.h"
#include "Fields/SU2Field.h"

#include "Actions/Action.h"
#include "Observables/AvgPlaquette.h"
#include "Integrators/LeapFrog.h"
#include "HMC/HMC_gauge.h"


int GridThread::_threads = 1;

struct Params {
    float beta;
    int nsteps;
    float epsilon;
    int Lsteps;
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

        ss >> p.beta >> comma >> p.nsteps >> comma >> p.epsilon >> comma >> p.Lsteps;
        params.push_back(p);
    }
    return params;
}


int main() {

    /***************
     * File Setup
     ***************/
    std::string input_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/temp_data/SU2/input_params.csv";
    std::vector<Params> params = read_params(input_filename);


    std::string output_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/Results/SU2/phases.csv";
    std::ofstream fout(output_filename);


    if (!fout)
        throw std::runtime_error("Error: Could not open results file: " + output_filename);

    fout << "beta,nsteps,avg_plaq,abs avg_plaq\n";

    /***************
     * Runner
     ***************/

    GridThread::SetMaxThreads();
    std::cout << "Threads: " << GridThread::GetThreads() << std::endl;

    constexpr int dimension = 4;
    constexpr std::array<int, dimension> dims({4, 4, 4, 4});

    GridBase<dimension> Grid(dims);

    for( Params P : params) {
        float beta = P.beta;
        int nsteps = P.nsteps;
        float epsilon = P.epsilon;
        int Lsteps = P.Lsteps;


        std::cout << "Running for "
                  << "nsteps = " << P.nsteps << ", "
                  << "beta = " << P.beta << ", "
                  << "epsilon = " << P.epsilon <<", "
                  << "Lsteps = " << P.Lsteps << ", "
                  << std::endl;

        std::random_device rd;
        unsigned int base_seed = rd();


        const RandNormal Phirng(0, 1, 43); // copied in HMC, will be reseeded in each thread
        const RandNormal Prng(0,1, 42); // copied in HMC, will be reseeded in each thread

        //using TheField = U1Field<dimension>;
        using TheField = SU2Field<dimension>;
        using TheAction = WilsonAction<TheField>;
        using TheIntegrator = LeapFrog<TheField, TheAction, dimension>;

        TheAction action(beta);
        TheIntegrator integrator(action, epsilon, Lsteps);
        AvgPlaquette avg_plaq(nsteps);

        HMCGauge<dimension, TheField, TheIntegrator, RandNormal, RandNormal, AvgPlaquette> hmc(nsteps, integrator, Phirng, Prng,
                                                                                    base_seed, &Grid, avg_plaq);

        auto t0 = std::chrono::steady_clock::now();
        hmc.Run();


        auto t1 = std::chrono::steady_clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "hmc.Run() wall time: " << dt.count() << " s\n";

        std::cout << "Acceptance Ratio: " << hmc.log().acceptance_ratio() << std::endl;

        std::cout << "mean " << avg_plaq.mean(0.1) << std::endl;

        fout << P.beta << ","
              << P.nsteps << ","
              << avg_plaq.mean(0.1) << ","
              << avg_plaq.abs_mean(0.1) << ","
              << "\n";
    }

    std::cout << "Done." << std::endl;

}
