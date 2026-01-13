cmake_minimum_required(VERSION 3.16)
project(Grid CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(run_sim run_sim.cpp
                    Grid.h
                    lattice.h
        Fields/ScalarField.h
        ../scratch/VectorField.h
                    Stencil.h
                    Expression.h
                    pragma_shorts.h
                    Action.h
        Reductions.h
        Action.h
        HMC.h
        RNG.h
        Observables/Magnetization.h
        io.cpp
        ../scratch/GaugeField.h
        ColorMatrix.h
        Fields/U1Field.h
        Fields/FieldBase.h
        HMC_gauge.h
        ArrayOpOverloads.h
        Observables/AvgPlaquette.h
        Observables/ObservableBase.h
        Integrators/LeapFrog.h
        Fields/SU2Field.h
)

add_executable(run_sim_gauge run_sim_gauge.cpp
        Grid.h
        lattice.h
        Fields/ScalarField.h
        ../scratch/VectorField.h
        Stencil.h
        Expression.h
        pragma_shorts.h
        Action.h
        Reductions.h
        Action.h
        HMC.h
        RNG.h
        Observables/Magnetization.h
        io.cpp
        ../scratch/GaugeField.h
        ArrayOpOverloads.h
)

#add_executable(tests tests.cpp
#        Grid.h
#        lattice.h
#        ScalarField.h
#        VectorField.h
#        Stencil.h
#        Expression.h
#        pragma_shorts.h
#)

find_package(OpenMP REQUIRED)
target_link_libraries(run_sim PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(run_sim_gauge PRIVATE OpenMP::OpenMP_CXX)

find_package(Eigen3 REQUIRED)
target_link_libraries(run_sim PRIVATE Eigen3::Eigen)
target_link_libraries(run_sim_gauge PRIVATE Eigen3::Eigen)

target_compile_options(run_sim PRIVATE -march=native)
target_compile_options(run_sim_gauge PRIVATE -march=native)


