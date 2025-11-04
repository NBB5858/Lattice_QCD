#include "lattice.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

int main() {

    std::ifstream fin("../data/params.csv");
    if (!fin) {
        std::cerr << "Error: could not open params.csv\n";
        return 1;
    }


    std::string line;
    std::getline(fin, line);

    int counter = 0;
    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        int N;
        double beta, B, J;
        long long iter;
        char comma;

        ss >> N >> comma >> beta >> comma >> B >> comma >> J >> comma >> iter;

        // --- Now run simulation for this line --- //
        std::cout << "Running N=" << N
                  << ", beta=" << beta
                  << ", B=" << B
                  << ", J=" << J
                  << ", iter=" << iter << std::endl;

        // Output file per run
        std::string outname = "../data/magnetization_" + std::to_string(counter) + ".csv";

        std::ofstream file(outname);
        if (!file) {
            std::cerr << "Error: could not open " << outname << "\n";
            continue;
        }

        IsingLattice1D lattice(N, beta, J, B);
        double mag = lattice.mag();

        file << mag;

        double mag_sum = mag;
        for( long long _ = 1; _ < iter; _++) {
            lattice.metro_iter();

            mag = lattice.mag();
            mag_sum += mag;

            file << "," << mag;
        }

        std::cout << "Average magnetization = "
                      << (mag_sum / iter) << "\n\n";

        ++counter;
    }
    return 0;
}