//
// Created by Noah bittermann on 11/21/25.
//

// understand what this is actually doing
// cmake -S grid_test -B grid_test/build \
//   -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
//   -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++
// cmake --build grid_test/build


// understand what this is actually doing
// cmake -S grid -B grid/build \
//   -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
//   -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++
// cmake --build grid/build

// Goal: define a scalar field theory on the latt

// should define a scalar action which lives on a lattice.
// operations on the scalar lattice should be done using thread_for


#include <vector>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <cmath>
#include <cassert>






// template<typename FieldType>
// class Integrator {
// public:
//
//     Integrator(GridBase* grid, Action& action) : P(grid), _action(action) {};
//
//     void Refresh(FieldType& U, rng) {
//         FieldType::generate_momenta(FieldType& P, rng);
//     };
//
// private:
//     FieldType P;
//     Action& _action;
// };
//
// class Observeable {};
//
//
// // Build an HMC
// template<typename IntegratorType>
// class HMC {
// public:
//
//     HMC(GridBase* grid) : _Ugrid(grid) {};
//     IntegratorType::FieldType U(_Ugrid);
//
//     HMC(int nsteps) : _nsteps(nsteps) {};
//     // construct with obserables;
//
//     void hmc_step() {
//         TheIntegrator.Refresh(U, srand);
//
//
//
//     }
//
//     void measure() {
//
//     }
//
//     // start the HMC procedure
//     void Run() {
//         for(int i = 0; i < _nsteps; ++i) {
//             hmc_step();
//             measure();
//         }
//     }
//
// private:
//     int _nsteps;
//     IntegratorType TheIntegrator;
//     GridBase* _Ugrid;
// };
//
//



