
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <array>

#include "Grid.h"
#include "Action.h"
#include "Observables/AvgPlaquette.h"
#include "Integrators/LeapFrog.h"
#include "HMC_gauge.h"

int GridThread::_threads = 1;

struct Params {
    double beta;
    int nsteps;
    double epsilon;
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
    std::string input_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/temp_data/Gauge/input_params.csv";
    std::vector<Params> params = read_params(input_filename);

    // std::string output_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/stored_data/phases/5x5.csv";
    // std::ofstream fout(output_filename);

    // std::string output_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/temp_data/Gauge/burn_in.csv";
    // std::ofstream fout(output_filename);

    std::string output_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/temp_data/Gauge/phases_check.csv";
    std::ofstream fout(output_filename);


    if (!fout)
        throw std::runtime_error("Error: Could not open results file: " + output_filename);

    fout << "beta,nsteps,avg_plaq,abs avg_plaq\n";
    //fout << "avg_plaq\n";

    /***************
     * Runner
     ***************/

    GridThread::SetMaxThreads();
    std::cout << "Threads: " << GridThread::GetThreads() << std::endl;

    constexpr int dimension = 4;
    constexpr std::array<int, dimension> dims({5, 5, 5, 5});

    GridBase<dimension> Grid(dims);

    for( Params P : params) {
        double beta = P.beta;
        int nsteps = P.nsteps;
        double epsilon = P.epsilon;
        int Lsteps = P.Lsteps;


        std::cout << "Running for "
                  << "nsteps = " << nsteps << ", "
                  << "beta = " << beta << ", "
                  << "epsilon = " << epsilon <<", "
                  << "Lsteps = " << Lsteps << ", "
                  << std::endl;

        std::random_device rd;
        unsigned int base_seed = rd();


        const RandNormal Prng(0,1, 42); // copied in HMC, will be reseeded in each thread
        const RandUniform Phirng(0, 2*M_PI, 43); // copied in HMC, will be reseeded in each thread

        using TheField = U1Field<dimension>;
        using TheAction = WilsonAction<TheField>;
        using TheIntegrator = LeapFrog<TheField, TheAction, dimension>;

        TheAction action(beta);
        TheIntegrator integrator(action, epsilon, Lsteps);
        AvgPlaquette avg_plaq(nsteps);

        HMCGauge<dimension, TheField, TheIntegrator, RandUniform, RandNormal, AvgPlaquette> hmc(nsteps, integrator, Phirng, Prng,
                                                                                    base_seed, &Grid, avg_plaq);

        hmc.Run();
        std::cout << "Acceptance Ratio: " << hmc.log().acceptance_ratio() << std::endl;

        std::cout << "mean " << avg_plaq.mean(0.1) << std::endl;
        //for( double obs : avg_plaq.cache() ) fout << obs << "," << "\n";

        fout << beta << ","
              << nsteps << ","
              << avg_plaq.mean(0.1) << ","
              << avg_plaq.abs_mean(0.1) << ","
              << "\n";
    }

    std::cout << "Done." << std::endl;

}
