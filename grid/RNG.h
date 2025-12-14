

#ifndef RNG_H
#define RNG_H

#include <random>


class RandNormal {
public:
    RandNormal(double mean=0.0, double stddev=1.0, unsigned int seed=42) :
    _gen(seed),
    _mean(mean),
    _stddev(stddev),
    _dist(mean, stddev) {};

    double draw() { return _dist(_gen);}
    void reseed(int seed) {_gen.seed(seed);}

private:
    std::mt19937 _gen;
    double _mean;
    double _stddev;
    std::normal_distribution<double> _dist;
};

class RandUniform {
public:
    RandUniform(double low=0.0, double high=1.0, unsigned int seed=42) :
    _low(low),
    _high(high),
    _gen(seed),
    _dist(low, high) {};

    double draw() { return _dist(_gen);}
    void reseed(int seed) {_gen.seed(seed);}

private:
    std::mt19937 _gen;
    double _low;
    double _high;
    std::uniform_real_distribution<double> _dist;
};



#endif //RNG_H
