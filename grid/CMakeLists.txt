cmake_minimum_required(VERSION 3.16)
project(Grid CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(run_sim run_sim.cpp
                    Grid.h
                    lattice.h
        Fields/ScalarField.h
                    VectorField.h
                    Stencil.h
                    Expression.h
                    pragma_shorts.h
                    Action.h
        Reductions.h
        Action.h
        HMC.h
        RNG.h
        observables.h
        io.cpp
        GaugeField.h
        ColorMatrix.h
        Fields/U1Field.h
        Fields/FieldBase.h
)

add_executable(run_sim_gauge run_sim_gauge.cpp
        Grid.h
        lattice.h
        Fields/ScalarField.h
        VectorField.h
        Stencil.h
        Expression.h
        pragma_shorts.h
        Action.h
        Reductions.h
        Action.h
        HMC.h
        RNG.h
        observables.h
        io.cpp
        GaugeField.h
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
target_link_libraries(run_sim_gauge PRIVATE Eigen3::Eigen)
target_link_libraries(run_sim PRIVATE Eigen3::Eigen)

target_compile_options(run_sim_gauge PRIVATE -O3 -march=native)
target_compile_options(run_sim PRIVATE -O3 -march=native)



