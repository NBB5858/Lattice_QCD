#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <complex>
#include <cmath>

#include "Grid.h"
#include "Lattice.h"
#include "Fields/ScalarField.h"
#include "Fields/U1Field.h"
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

    constexpr int d = 2;
    constexpr std::array<int, d> dims({3, 3});

    GridBase<d> Grid(dims);

    // Lattice<ScalarField, d> phi(&Grid);

    const RandUniform Phirng(0,2*M_PI, 42); // copied in HMC, will be reseeded in each thread
    std::vector<RandUniform> Phirngs{Phirng};

    // phi.RandomizeLattice(Phirngs);
    //
    // phi.Print2D();

    //Lattice<ScalarField, d> output(&Grid);


    Lattice<U1Field<d>, d> U(&Grid);

    U.RandomizeLattice(Phirngs);


    U.Print2D();

    //double out = sum(plaqReTraceSum(U), &Grid);

    VectorField<std::array<CMat<1>, d>, d> staple_test;

    TA(U*stapleSum(U));


    //output.Print2D();






    // phi.put(0, 1);
    // phi.put(1, 2);
    // phi.put(2, 5);
    // phi.put(6, 7.8);
    //
    // phi.Print2D();
    //
    // Lattice<ScalarField, d> Gamma(&Grid);
    //
    // Gamma = dot(grad(phi), grad(phi));
    //
    // Gamma.Print2D();
    //
    // ScalarAction action(1, 1);
    //
    // Gamma = action.Force(phi);
    // Gamma.Print2D();
    //
    // std::cout << action.S(phi) << std::endl;


    //Gamma.Print2D();


    //
    // Lattice<VectorField<Z2GroupElement, d>, d> U(&Grid);
    //
    // RandSign rng(40);
    //
    // U.Randomize(0, rng);
    // U.Randomize(2, rng);
    // U.Randomize(4, rng);
    // U.Randomize(6, rng);
    // U.Randomize(8, rng);
    //
    // U.Print2D();
    //
    // GaugeAction<Z2GroupElement> action(1.0);
    //
    // double output = action.S(U);
    // std::cout << output << std::endl;
    //
    //
    // // Lattice<ScalarField, d> phi(&Grid);
    // //
    // // phi = plaqTraceSum(U);
    // //
    // // phi.Print2D();
    // //
    // // double output = sum(phi, &Grid).val();
    //
    // std::cout << output << std::endl;




    // mu=0,1 and nu=0,1. U{00} I guess should be zero. U_{01} should be

    // fout << mass2 << ","
    //       << lambda << ","
    //       << nsteps << ","
    //       << mag.mean(0.1) << ","
    //       << mag.abs_mean(0.1) << ","
    //       << "\n";


    std::cout << "Done." << std::endl;

}

