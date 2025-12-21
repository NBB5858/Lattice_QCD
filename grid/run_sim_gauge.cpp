#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

#include "Grid.h"
#include "Lattice.h"
#include "ScalarField.h"
#include "VectorField.h"
#include "Expression.h"
#include "Reductions.h"
#include "Action.h"
#include "GaugeField.h"
#include "Observables.h"
#include "HMC.H"


int GridThread::_threads = 1;

struct Params {
    double mass2;
    double lambda;
    int nsteps;
    double epsilon;
    double Lsteps;
    double sigma;
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
    // std::string input_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/temp_data/input_params.csv";
    // std::vector<Params> params = read_params(input_filename);

    // std::string output_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/stored_data/phases/5x5.csv";
    // std::ofstream fout(output_filename);

    // std::string output_filename = "/Users/noahbittermann/CLionProjects/Lattice_QCD/temp_data/burn_in.csv";
    // std::ofstream fout(output_filename);


    //
    // if (!fout)
    //     throw std::runtime_error("Error: Could not open results file: " + output_filename);
    //
    // fout << "m^2,lambda,nsteps,mag,abs mag\n";
    //fout << "mag\n";

    /***************
     * Runner
     ***************/

    GridThread::SetMaxThreads();
    std::cout << "Threads: " << GridThread::GetThreads() << std::endl;

    std::array<int, 2> dims({5, 5});
    GridBase<2> Grid(dims);


    Lattice<GaugeField<Z2GroupElement>, 2> U(&Grid);




        // fout << mass2 << ","
        //       << lambda << ","
        //       << nsteps << ","
        //       << mag.mean(0.1) << ","
        //       << mag.abs_mean(0.1) << ","
        //       << "\n";


    std::cout << "Done." << std::endl;

}

