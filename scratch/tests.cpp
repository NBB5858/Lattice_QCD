V = V1 + V2;


//Sum = -phi // test me; may need a unary operator

// phi1 = phi2; you never define an assignment operator tht takes in a lattice, jsut a lattice expression. test this


// Stencil test({2,2});
//
// // std::vector<int> coords;
// // coords = test.ravel(3);
// // std::cout << coords[0] << "," << coords[1] << std::endl;
//
// std::cout << test.unravel({0,0}) << std::endl;
// std::cout << test.unravel({0,1}) << std::endl;
// std::cout << test.unravel({1,0}) << std::endl;
// std::cout << test.unravel({1,1}) << std::endl;


// Stencil test({2,2});
//
// std::cout << test.neighbors[0].m[0] << std::endl;

// std::vector<int> coords;
// coords = test.ravel(3);
// std::cout << coords[0] << "," << coords[1] << std::endl;

// std::cout << test.unravel({0,0}) << std::endl;
// std::cout << test.unravel({0,1}) << std::endl;
// std::cout << test.unravel({1,0}) << std::endl;
// std::cout << test.unravel({1,1}) << std::endl;


// phi.Print2D();
//
// phi.put(0, ScalarField(2.5));
// phi.put(7, ScalarField(2.7));
// phi.put(8, ScalarField(3.0));
//
// std::cout << std::endl;
// phi.Print2D();

// phi.put(0, ScalarField(2.5));
// phi.put(7, ScalarField(2.7));
// phi.put(8, ScalarField(3.0));

// std::cout << std::endl;
// phi.Print2D();

// double integral = sum(phi + phi, &Grid);
// std::cout << integral << std::endl;



// Lattice<ScalarField> Gamma(&Grid);

// Gamma = dot(grad(phi), grad(phi));

// Gamma.Print2D();
// std::cout << std::endl;

// double integral = sum(dot(grad(phi), grad(phi)), &Grid);
// std::cout << integral << std::endl;

// check sum on different kinds of expressions
// do I want to be able to do sum on vectors? I don't think so, that wouldn't a double.


// std::cout << std::endl;

// Lattice<ScalarField> phi(&Grid);
//
// phi.put(0, ScalarField(2.5));
// phi.put(7, ScalarField(2.7));
// phi.put(8, ScalarField(3.0));
//
// phi.Print2D();
//
// std::cout << std::endl;
//
// Lattice<ScalarField> Gamma(&Grid);
//
// Gamma = box(phi);
//
// Gamma.Print2D();


//
// phi.put(0, ScalarField(2.5));
// phi.put(7, ScalarField(2.7));
// phi.put(8, ScalarField(3.0));
//
// phi.Print2D();
//
// std::cout << std::endl;
//
// Lattice<VectorField<ScalarField>> Gamma(&Grid);
// Lattice<VectorField<ScalarField>> Rho(&Grid);
//
// Gamma = grad(phi);
// Gamma.Print2D();
//
// std::cout << std::endl;
//
// Rho = box(Gamma);
// Rho.Print2D();


// //auto start = std::chrono::high_resolution_clock::now();
    // hmc.Run();
    //
    // //auto stop = std::chrono::high_resolution_clock::now();
    // //auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);