#ifndef ACTION_H
#define ACTION_H

#include "../Fields/ScalarField.h"
#include "../Fields/U1Field.h"
#include "../Lattice/Lattice.h"

class ScalarAction {
public:

    ScalarAction(float mass2, float lambda) : _mass2(mass2), _lambda(lambda) {}

    template<int d>
    float S(const Lattice<ScalarField, d>& phi) const {
        return sum(0.5f* grad_squared(phi)
                      + 0.5f * _mass2 * phi * phi
                      + (_lambda / 24.0f) * phi * phi * phi * phi,
                    phi.grid());
    }

    template<int d>
    auto Force(const Lattice<ScalarField, d>& phi) const  {
        return -(1.0f)*box(phi) + _mass2 * phi + (_lambda / 6.0f) * phi * phi * phi;
    }


private:
    float _mass2;
    float _lambda;
};


template<typename GaugeField>
class WilsonAction {
public:
    explicit WilsonAction(float beta) : _beta(beta) {};

    template<int d>
    float S(const Lattice<GaugeField, d>& U) const {
        const auto Nd = static_cast<float>(U.grid()->dims().size());
        const float nSites = static_cast<float>(U.grid()->nsites());
        const float Nplaq = nSites * Nd * (Nd - 1.0f) / 2.0f;
        return _beta * ( Nplaq - (1.0f/N) * sum(plaqReTraceSum(U), U.grid()) );
    }

    template<int d>
    auto Force(const Lattice<GaugeField, d>& U) const {
        return (_beta / N) * TA( U * adj(stapleSum(U)) );
    }

private:
    float _beta;
    float N = static_cast<float>(FieldTraits<GaugeField>::groupdim);
};








#endif //ACTION_H