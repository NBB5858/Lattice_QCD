

#ifndef RNG_H
#define RNG_H

#include <random>

class RNGBase {
public:
    explicit RNGBase(unsigned int seed) : _gen(seed) {}
    void reseed(int seed) { _gen.seed(seed); }

protected:
    std::mt19937 _gen;
};


class RandNormal : public RNGBase {
public:
    RandNormal(double mean, double stddev, unsigned int seed):
    RNGBase(seed),
    _mean(mean),
    _stddev(stddev),
    _dist(mean, stddev) {};

    double draw() { return _dist(_gen);}

private:
    double _mean;
    double _stddev;
    std::normal_distribution<double> _dist;
};

class RandUniform : public RNGBase {
public:
    RandUniform(double low, double high, unsigned int seed):
    RNGBase(seed),
    _low(low),
    _high(high),
    _dist(low, high) {}

    double draw() { return _dist(_gen);}

private:
    double _low;
    double _high;
    std::uniform_real_distribution<double> _dist;
};

class RandSign : public RNGBase {
public:
    explicit RandSign(unsigned int seed):
    RNGBase(seed),
    _dist(0.5) {}

    int draw() {return _dist(_gen) ? 1 : -1;}

private:
    std::bernoulli_distribution _dist;
};

#endif //RNG_H
