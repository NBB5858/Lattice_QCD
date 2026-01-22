#ifndef LEAPFROG_H
#define LEAPFROG_H

#include "../Lattice/Lattice.h"

template<typename FieldType, typename ActionType, int d>
class LeapFrog {

    using GaugeField = FieldType;
    using MomField   = typename FieldTraits<GaugeField>::MomentumField;

public:
    LeapFrog(ActionType action, float epsilon, int Lsteps) :
        _action(action),
        _epsilon(epsilon),
        _lsteps(Lsteps) {}

    void evolve(Lattice<FieldType, d>& phi, Lattice<MomField, d>& P) {
        P   -= (0.5f) * _action.Force(phi) * _epsilon;
        for(int i=0; i<_lsteps-1; ++i) {
            phi = FieldTraits<MomField>::exp(_epsilon * P) * phi;
            P   -= _action.Force(phi) * _epsilon;
        }
        phi = FieldTraits<MomField>::exp(_epsilon * P) * phi;
        P   -= (0.5f) * _action.Force(phi) * _epsilon;
    }

    // std::complex<float> H(Lattice<FieldType, d>& phi, Lattice<MomField, d>& P) {
    //     return (0.5f) * sum(trace(dot(adj(P), P)), P.grid()) + _action.S(phi);
    // }
     std::complex<float> H(Lattice<FieldType, d>& phi, Lattice<MomField, d>& P) {
        constexpr float N = static_cast<float>(FieldTraits<FieldType>::groupdim);
        return (0.5f / N) * sum(trace(dot(adj(P), P)), P.grid()) + _action.S(phi);
    }

private:
    ActionType _action;
    float _epsilon;
    int _lsteps;
};


#endif //LEAPFROG_H
