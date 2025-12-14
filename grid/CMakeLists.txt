cmake_minimum_required(VERSION 3.16)
project(Grid CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(run_sim run_sim.cpp
                    Grid.h
                    lattice.h
                    ScalarField.h
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


