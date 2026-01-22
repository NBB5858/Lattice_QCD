cmake_minimum_required(VERSION 3.16)
project(Grid CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(run_sim run_sim.cpp
                    Grid.h
                    Lattice/Lattice.h
                    Lattice/Expression.h
                    Fields/ScalarField.h
                    Stencil/Stencil.h
                    Actions/Action.h
                    Fields/FieldBase.h
                    Shifting/Shifting.h
                    openmp/pragma_shorts.h
                    openmp/Reductions.h
                    HMC/HMC.h
                    RNG/RNG.h
                    Observables/Magnetization.h
)

add_executable(run_sim_gauge run_sim_gauge.cpp
                            Grid.h
                            Lattice/Lattice.h
                            Lattice/Expression.h
                            Fields/ScalarField.h
                            Stencil/Stencil.h
                            openmp/pragma_shorts.h
                            openmp/Reductions.h
                            Actions/Action.h
                            Actions/Action.h
                            HMC/HMC.h
                            RNG/RNG.h
                            Observables/Magnetization.h
                            Fields/SU2Field.h
                            Fields/FieldBase.h
                            Shifting/Shifting.h
)

find_package(OpenMP REQUIRED)
target_link_libraries(run_sim PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(run_sim_gauge PRIVATE OpenMP::OpenMP_CXX)

find_package(Eigen3 REQUIRED)
target_link_libraries(run_sim PRIVATE Eigen3::Eigen)
target_link_libraries(run_sim_gauge PRIVATE Eigen3::Eigen)


find_package(sleef CONFIG REQUIRED)
target_link_libraries(run_sim PRIVATE sleef::sleef)
target_link_libraries(run_sim_gauge PRIVATE sleef::sleef)


target_compile_options(run_sim PRIVATE -march=native)
target_compile_options(run_sim_gauge PRIVATE -march=native)

target_compile_definitions(run_sim_gauge PRIVATE USE_SIMD_ACROSS_SITES)
target_compile_definitions(run_sim PRIVATE USE_SIMD_ACROSS_SITES)

target_compile_definitions(run_sim_gauge PRIVATE USE_OPENMP)
target_compile_definitions(run_sim PRIVATE USE_OPENMP)


