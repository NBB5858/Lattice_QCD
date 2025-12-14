
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

void write_results_csv(const std::string& filename,
                       const std::vector<double>& results)
{
    std::ofstream fout(filename);
    if (!fout)
        throw std::runtime_error("Error: Could not open results file: " + filename);

    for (size_t i = 0; i < results.size(); ++i) {
        fout << results[i] << "\n";
    }
}

