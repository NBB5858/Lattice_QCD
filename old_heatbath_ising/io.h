#ifndef IO_H
#define IO_H

#include <string>
#include <vector>

struct Params {
    int X;
    int Y;
    int Z;
    double beta;
    double B;
    double J;
    long long iter;
};

std::vector<Params> read_params(const std::string& filename);

void write_results_csv(const std::string& filename,
                       const std::vector<Params>& params,
                       const std::vector<double>& results);

#endif // IO_H
