#include "simd_types.h"
#include <vector>
#include <iostream>

// int main() {
//     constexpr int N = 1024;
//     std::vector<double> x(N), y(N);
//
//     for (int i = 0; i < N; ++i) { x[i] = 0.01 * i; y[i] = 1.0; }
//
//     // y = 2*x + 3*y
//     for (int i = 0; i < N; i += simd::W) {
//         auto xv = simd::loadu(&x[i]);
//         auto yv = simd::loadu(&y[i]);
//
//         auto out = simd::add(
//             simd::mul(simd::set1(2.0), xv),
//             simd::mul(simd::set1(3.0), yv)
//         );
//
//         simd::storeu(&y[i], out);
//     }
//
//     std::cout << "y[10]=" << y[10] << "\n";
// }

#include "simd_types.h"
#include <iostream>

#include "simd_types.h"
#include <iostream>
#include <complex>

int main() {
    using namespace simd;

    vReal x(2.0f);

    vReal result = sin(x);// * sin(x) + cos(x) * cos(x);

    alignas(16) float lanes[W];

    simd::storeu(lanes, result.v);

    std::cout << lanes[0] << std::endl;
    std::cout << lanes[1] << std::endl;
    std::cout << lanes[2] << std::endl;
    std::cout << lanes[3] << std::endl;

}
