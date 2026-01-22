# Lattice Gauge Code

This respository implements a small lattice gauge simulation framework in C++.  The design is largely based off of [Grid](https://github.com/paboyle/Grid), though there are some differences in some design choices.
The core features of the framework are

1. __Focus on Hamiltonian Monte Carlo__. This is not necessary for quenched QCD, where a heatbath algorithm is simpler and faster.  However, HMC will become mandatory when fermions are implemented later.
2. __Heavily Templating__.  This allows for compiler optimizations and avoids runtime polymorphism.  The important constants of the simulation (the lattice dimension and gauge group size) are known at compile time, and templating on these parameters allows the compiler to optimize small matrix operations.  We keep the implementation generic by defining a FieldTraits<FieldType> struct, which defines the Storage for a field type, and functions to view, randomize, and exponentiate the field.  
2. __Lazy Evaluation & Expression Templates__.  Mathematical expressions involving lattices such as $2 * \phi$ are not evaluated until assignment into another lattice $\psi = 2 * \phi$.  This avoids the creation of temporary objects which are expensive to calculate and store.  Unary and Binary operations are implemendted as simple classes, which contain instructions on how to evaluate an operation. These instructions are only executed on assignment.  
3. __SIMD Vectorization__. Lattice sites are grouped into blocks.  Each site is packed into a lane of a  SIMD register, holding 4 floats altogether (for NEON float32x4_t on a Mac M2).  We build containers CMat and Link which hold SIMD leaf types vComplex, vReal, allowing for site-parallel matrix operations.    
4. __OpenMP Parallelization__.  Blocks of sites are process in parallel by separate threads using OpenMP.  These threads are dispatched whenever the execution enters a thread_for block. There is actually a high degree of overhead for constanly forking/joining  threads. This is actually the bulk of the runtime when OpenMP is active, but the cost is worth it when amortized over large lattices.  It may be possibe to avoid forking/joining by restructing for loops, but this would make the library less generic.  
5. __Cshift__. Many lattice operations require a shifted field (e.g. computing derivatives, plaquettes).  Naively, one could build a temporary lattice for each shift. We avoid by implementing a Cshift function which accepts a block of lattice sites, and locates their neighbors' data using a stencil for arbitrary displacements.  This is critical for computing plaquettes and staples efficiently.    
6. __Views__. Because data storage is packed into blocks of SIMD langes, it is non-trivial to peek at their values. We do this by implementing View objects, wich hold pointers to the underlying data, and provide accessors for viewing/printing.

### Requirements / Dependencies
 - A C++17 compatible compiler (e.g. g++ 7+, clang++ 7+)
 - AArch64 / Apple Silicon (M!/M2) when "USE_SIMD_ACROSS_SITES" is enabled (see CMake file)
 - OpenMP - a thread often "sleeps" after finishing its tasks while waiting for other threads. This creates overhead, and can be reduce by setting the OMP_WAIT_POLICY environment variables to ACTIVE and KMP_BLOCKTIME (the time a thread waits before sleeping) to some high number.
 - SLEEF - SIMD math library used for trig functions

Some sample results can be found in [notebooks/Sample_Results.ipynb](notebooks/Sample_Results.ipynb)



