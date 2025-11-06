#include "io.h"
#include <fstream>
#include <sstream>
#include <stdexcept>


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

        ss >> p.X >> comma >> p.Y >> comma >> p.Z >> comma >> p.beta >> comma >> p.B >> comma >> p.J >> comma >> p.iter;
        params.push_back(p);
    }
    return params;
}

void write_results_csv(const std::string& filename,
                       const std::vector<Params>& params,
                       const std::vector<double>& results)
{
    std::ofstream fout(filename);
    if (!fout)
        throw std::runtime_error("Error: Could not open results file: " + filename);

    // header
    fout << "X, Y, Z, beta,B,J,iter,mag\n";

    // rows
    for (size_t i = 0; i < params.size(); ++i) {
        const Params& p = params[i];
        fout << p.X      << ","
             << p.Y      << ","
             << p.Z      << ","
             << p.beta   << ","
             << p.B      << ","
             << p.J      << ","
             << p.iter   << ","
             << results[i] << "\n";
    }
}
