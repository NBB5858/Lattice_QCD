#ifndef SCALARFIELD_H
#define SCALARFIELD_H

#include<iostream>

#include "../simd/simd_types.h"
#include "FieldBase.h"
#include "../RNG/RNG.h"

struct ScalarField {};

struct ScalarFieldView {
    simd::vReal* blk = nullptr;
    int lane = 0;

    float get() const { return simd::extract_lane(*blk, lane); }
    void set(float v) const { simd::insert_lane(*blk, lane, v); }
};


struct ScalarFieldConstView {
    const simd::vReal* blk = nullptr;
    int lane = 0;

    float get() const { return simd::extract_lane(*blk, lane);}
    void print() const { std::cout << get() << ","; }
};

template<>
struct FieldTraits<ScalarField> {
    using Storage   = simd::vReal;
    using View = ScalarFieldView;
    using ConstView = ScalarFieldConstView;
    using MomentumField = ScalarField;


    static View view(Storage& x, int lane) { return View{ &x, lane }; }
    static ConstView view(const Storage& x, int lane) { return ConstView{ &x, lane }; }

    static void Randomize(simd::vReal& x, RandNormal& rng) {
        x = simd::rand_vReal(rng);
    }
};

#endif //SCALARFIELD_H
