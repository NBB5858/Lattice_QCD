#include <random>
#include <fstream>
#include <iostream>

struct Params {
    int    N        = 20;
    double temp     = 1.0;
    double B   = 1.0;
    double J = 1.0;
    int    iter     = 1'000'000;
};

int rand_spin();
int rand_site( int lattice_size );
double compute_magnetization( int const *lattice);


int main() {

    Params P;

    int N = P.N;
    double temp = P.temp;
    double B = P.B;
    double J = P.J;
    double iter = P.iter;

    std::mt19937_64 rng;

    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);

    std::ofstream file("../data/magnetization.csv", std::ios::binary);

    rng.seed(42);

    int lattice[N];
    for(int i = 0; i < N; ++i)
        lattice[i] = rand_spin();

    double energy = 0;
    for(int i = 1; i < N; i++)
        energy += lattice[i] * lattice[i-1];

    double magnetization = compute_magnetization(lattice);
    file << magnetization;

    double mag_sum = magnetization / iter;
    for( int _ = 1; _ < iter; _++) {
        int i = rand_site(N);

        double old_link_energy;
        double new_link_energy;

        int new_spin = rand_spin();

        if(i == 0) {
            old_link_energy = J * lattice[i] * lattice[i+1] + B * lattice[i];
            new_link_energy = J * new_spin * lattice[i+1] + B * new_spin;
        }
        else if(i == N-1) {
            old_link_energy = J * lattice[i] * lattice[i-1] + B * lattice[i];
            new_link_energy = J * new_spin * lattice[i-1] + B * new_spin;
        }
        else {
            old_link_energy = J * lattice[i] * (lattice[i-1] + lattice[i+1]) + B * lattice[i];
            new_link_energy = J * new_spin * (lattice[i-1] + lattice[i+1]) + B * new_spin;
        }

        double new_energy = energy - old_link_energy + new_link_energy;
        double accep_prob = std::exp(-(new_energy - energy)/temp);

        double u = uniform_dist(rng);
        if (u <= accep_prob)
            lattice[i] = new_spin;

        magnetization = compute_magnetization(lattice);
        mag_sum += magnetization/iter;
        file << "," << magnetization;
    }

    std::cout << mag_sum / N << std::endl;
}

int rand_spin() {
    static std::mt19937 gen(std::random_device{}());
    static std::uniform_int_distribution<int> dist(0, 1);
    return dist(gen) == 0 ? -1 : 1;
}

int rand_site( int lattice_size) {
    static std::mt19937 gen(std::random_device{}());
    static std::uniform_int_distribution<int> dist(0, lattice_size-1);
    return dist(gen);
}

double compute_magnetization( int const *lattice ) {

    int N = sizeof(lattice);

    double magnetization = 0;
    for(int i = 0; i < N; i++)
        magnetization += lattice[i];

    return magnetization / N;
}
