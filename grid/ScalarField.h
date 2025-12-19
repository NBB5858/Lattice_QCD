
#ifndef SCALARFIELD_H
#define SCALARFIELD_H

#include <iostream>

class ScalarField {
public:
    ScalarField() : _val(0.0){}
    explicit ScalarField(int /*d*/) : _val(0.0){} // Get rid of this ASAP
    explicit ScalarField(int /*d*/, double val) : _val(val){} // Get rid of this ASAP
    explicit ScalarField(double val) : _val(val){}
    void print() const { std::cout << _val; }

    double val() const { return _val; }

    template<typename RNG>
    void Randomize(RNG& rng) {
        _val = rng.draw();
    }

    ScalarField& operator +=(const ScalarField& other) {
        this->_val += other._val;
        return *this;
    }

    ScalarField& operator -=(const ScalarField& other) {
        this->_val -= other._val;
        return *this;
    }

    friend ScalarField operator+(const ScalarField&  phi1, const ScalarField& phi2);
    friend ScalarField operator-(const ScalarField& phi1, const ScalarField& phi2);
    friend ScalarField operator*(const ScalarField& phi1, const ScalarField& phi2);
    friend ScalarField operator/(const ScalarField& phi1, const ScalarField& phi2);

    friend ScalarField operator*(double c, const ScalarField& phi);
    friend ScalarField operator*(const ScalarField& phi, double c);
    friend ScalarField operator/(const ScalarField& phi, double c);

private:
    double _val;
};

ScalarField inline operator+(const ScalarField& phi1, const ScalarField& phi2){return ScalarField(phi1._val + phi2._val);}
ScalarField inline operator-(const ScalarField& phi1, const ScalarField& phi2){return ScalarField(phi1._val - phi2._val);}
ScalarField inline operator*(const ScalarField& phi1, const ScalarField& phi2){return ScalarField(phi1._val * phi2._val);}
ScalarField inline operator/(const ScalarField& phi1, const ScalarField& phi2){return ScalarField(phi1._val / phi2._val);}

ScalarField inline operator*(double c, const ScalarField& phi){return ScalarField(c*phi._val);}
ScalarField inline operator*(const ScalarField& phi, double c){return ScalarField(c*phi._val);}
ScalarField inline operator/(const ScalarField& phi, double c){return ScalarField(phi._val/c);}




#endif //SCALARFIELD_H
