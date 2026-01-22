#ifndef SIMD_TYPES_H
#define SIMD_TYPES_H

#include <sleef.h>
#include <complex>
#include <cmath>

namespace simd {
#ifdef USE_SIMD_ACROSS_SITES

#if !defined(__aarch64__)
#error "This NEON backend assumes aarch64 architecture"
#endif

#include <arm_neon.h>

    static constexpr int W = 4;
    using vfloat = float32x4_t;            // ARM Neon type; SIMD register holding 4 lanes of floats

    inline vfloat set1(float x) { return vdupq_n_f32(x); }
    inline vfloat loadu(const float* p) { return vld1q_f32(p); }        // load 4 floats into register
    inline void   storeu(float* p, vfloat v) { vst1q_f32(p, v); }       // store 4 floats in memory

    inline vfloat add(vfloat a, vfloat b) { return vaddq_f32(a, b); }
    inline vfloat sub(vfloat a, vfloat b) { return vsubq_f32(a, b); }
    inline vfloat mul(vfloat a, vfloat b) { return vmulq_f32(a, b); }
    inline vfloat div(vfloat a, vfloat b) { return vdivq_f32(a, b); }
    inline vfloat neg(vfloat a) { return vnegq_f32(a); }
    inline vfloat sqrt(vfloat a) { return vsqrtq_f32(a); }
    inline vfloat square(vfloat a) { return mul(a, a); }
    inline vfloat sin(vfloat x) { return Sleef_sinf4_u10(x); }
    inline vfloat cos(vfloat x) { return Sleef_cosf4_u10(x); }

    // horizontal sum of vector elements
    inline float hsum(vfloat a) { return vaddvq_f32(a); }

    // Lane manipulation
    inline vfloat rotl1(vfloat a) { return vextq_f32(a, a, 1); }                 // rotate left by 1 lane
    inline vfloat shift_plus1(vfloat cur, vfloat next) { return vextq_f32(cur, next, 1); } // [cur1,cur2,cur3,next0]
    inline vfloat shift_minus1(vfloat prev, vfloat cur) { return vextq_f32(prev, cur, 3); } // [prev3,cur0,cur1,cur2]

#else

    static constexpr int W = 1;
    using vfloat = float;

    inline vfloat set1(float x) { return x; }
    inline vfloat loadu(const float* p) { return *p; }
    inline void  storeu(float* p, vfloat v) { *p = v; }

    inline vfloat add(vfloat a, vfloat b) { return a+b; }
    inline vfloat sub(vfloat a, vfloat b) { return a-b; }
    inline vfloat mul(vfloat a, vfloat b) { return a*b; }
    inline vfloat div(vfloat a, vfloat b) { return a/b; }
    inline vfloat neg(vfloat a) {return -a; }
    inline vfloat sqrt(vfloat a) { return std::sqrt(a); }
    inline vfloat square(vfloat a) { return mul(a, a); }
    inline vfloat cos(vfloat a) { return std::cos(a); }
    inline vfloat sin(vfloat a) { return std::sin(a); }

    inline vfloat hsum(vfloat a) { return a; }

    //TODO: check if rotation is correct
    inline vfloat rotl1(vfloat a) { return a; }
    inline vfloat shift_plus1(vfloat cur, vfloat next) { return next; }
    inline vfloat shift_minus1(vfloat prev, vfloat cur) { return prev; }

#endif // USE_SIMD_ACROSS_SITES

    struct vReal {
        vfloat v;
        vReal() : v(set1(0.0f)) {};

        explicit vReal(float x) : v(set1(x)) {}

#ifdef USE_SIMD_ACROSS_SITES
        // when SIMD is off, this collides with prior constructor
        explicit vReal(vfloat r) : v(r) {}
#endif

        static vReal loadu(const float* p) { return vReal{ simd::loadu(p) }; }
        static void  storeu(float* p, vReal x) { simd::storeu(p, x.v); }
    };

    inline vReal operator+(vReal a, vReal b) { return vReal{ add(a.v, b.v) }; }
    inline vReal operator-(vReal a, vReal b) { return vReal{ sub(a.v, b.v) }; }
    inline vReal operator*(vReal a, vReal b) { return vReal{ mul(a.v, b.v) }; }
    inline vReal operator/(vReal a, vReal b) { return vReal{ div(a.v, b.v) }; }
    inline vReal operator-(vReal a) {return vReal{ neg(a.v) }; }

    inline vReal operator*(vReal a, float b) { return a * vReal(b); }
    inline vReal operator*(float b, vReal a) { return vReal(b) * a; }
    inline vReal operator/(vReal a, float b) { return a / vReal{b};}

    inline vReal square(vReal a){ return vReal{ square(a.v) }; }
    inline vReal sqrt(vReal a){ return vReal{ sqrt(a.v) }; }
    inline vReal cos(vReal a) { return vReal{ cos(a.v) }; }
    inline vReal sin(vReal a) { return vReal{ sin(a.v) }; }

    inline vReal& operator+=(vReal& a, vReal b) { a = a+b; return a; }
    inline vReal& operator-=(vReal& a, vReal b) { a = a-b; return a; }

    inline vReal conj(vReal a) { return a; }

    //TODO: Not truly vectorized; fix later
    template<class RNG>
    inline vReal rand_vReal(RNG& rng) {
        alignas(16) float tmp[W];
        for (int lane = 0; lane < W; ++lane) tmp[lane] = static_cast<float>(rng.draw());
        return vReal(simd::loadu(tmp));
    }

    inline float sum_lanes(vReal a) { return hsum(a.v); }

    inline vReal rotl1(vReal a) { return vReal{ simd::rotl1(a.v) }; }
    inline vReal shift_plus1(vReal cur, vReal next) { return vReal{ simd::shift_plus1(cur.v, next.v) }; }
    inline vReal shift_minus1(vReal prev, vReal cur) { return vReal{ simd::shift_minus1(prev.v, cur.v) }; }


#ifdef USE_SIMD_ACROSS_SITES
    inline float extract_lane(const vReal& x, int lane) {
        switch (lane) {
            case 0: return vgetq_lane_f32(x.v, 0);
            case 1: return vgetq_lane_f32(x.v, 1);
            case 2: return vgetq_lane_f32(x.v, 2);
            default: return vgetq_lane_f32(x.v, 3);
        }
    }

    inline void insert_lane(vReal& x, int lane, float val) {
        switch (lane) {
            case 0: x.v = vsetq_lane_f32(val, x.v, 0); break;
            case 1: x.v = vsetq_lane_f32(val, x.v, 1); break;
            case 2: x.v = vsetq_lane_f32(val, x.v, 2); break;
            default: x.v = vsetq_lane_f32(val, x.v, 3); break;
        }
    }
#else
    inline float extract_lane(const vReal& x, int /*lane*/) { return x.v; }
    inline void  insert_lane(vReal& x, int /*lane*/, float val) { x.v = val; }
#endif

    struct vComplex {
        vReal re;
        vReal im;

        vComplex() = default;
        vComplex(vReal r, vReal i) : re(r), im(i) {}
        vComplex(vfloat r, vfloat i) : re(vReal{r}), im(vReal{i}) {}
        explicit vComplex(std::complex<float> z) : re(z.real()), im(z.imag()) {}

        static vComplex zero() { return vComplex{ vReal(0.0f), vReal(0.0f) }; }
    };

    inline vComplex operator+(vComplex a, vComplex b) { return vComplex{ a.re + b.re, a.im + b.im }; }
    inline vComplex operator-(vComplex a, vComplex b) { return vComplex{ a.re - b.re, a.im - b.im }; }
    inline vComplex operator*(vComplex a, vComplex b) {
        vReal ac = a.re * b.re;
        vReal bd = a.im * b.im;
        vReal ad = a.re * b.im;
        vReal bc = a.im * b.re;
        return vComplex{ ac - bd, ad + bc };
    }

    inline vComplex operator*(vComplex a, float b) { return vComplex{ a.re * b, a.im * b }; }
    inline vComplex operator*(float b, vComplex a) { return a * b; }
    inline vComplex operator/(vComplex a, float b) { return vComplex{ a.re / b, a.im / b }; }

    inline vComplex operator*(vComplex a, vReal b) { return vComplex{ a.re * b, a.im * b }; }
    inline vComplex operator*(vReal b, vComplex a) { return a * b; }
    inline vComplex operator/(vComplex a, vReal b) { return vComplex{ a.re / b, a.im / b }; }

    inline vComplex& operator+=(vComplex& a, vComplex b) { a = a+b; return a; }
    inline vComplex& operator-=(vComplex& a, vComplex b) { a = a-b; return a; }

    inline vComplex conj(vComplex a) { return vComplex{ a.re, -a.im}; }

    //TODO: Replace with fma operation
    inline vReal norm2(vComplex a) { return a.re * a.re + a.im * a.im; }

    inline std::complex<float> sum_lanes(vComplex a) { return {hsum(a.re.v), hsum(a.im.v)}; }

    inline vComplex rotl1(const vComplex& a) { return vComplex{ rotl1(a.re), rotl1(a.im) }; }

    inline vComplex shift_plus1(const vComplex& cur, const vComplex& next) {
        return vComplex{ shift_plus1(cur.re, next.re),
                         shift_plus1(cur.im, next.im) };
    }

    inline vComplex shift_minus1(const vComplex& prev, const vComplex& cur) {
        return vComplex{ shift_minus1(prev.re, cur.re),
                         shift_minus1(prev.im, cur.im) };
    }

    inline vComplex loadu_soa(const float* re, const float* im) {
        return vComplex{ vReal::loadu(re), vReal::loadu(im) };
    }
    inline void storeu_soa(float* re, float* im, vComplex z) {
        vReal::storeu(re, z.re);
        vReal::storeu(im, z.im);
    }

    // ---------------------------
    // FAST vComplex lane ops (built from vReal lane ops)
    // ---------------------------
    inline std::complex<float> extract_lane(const vComplex& z, int lane) {
        float r = extract_lane(z.re, lane);
        float i = extract_lane(z.im, lane);
        return {r, i};
    }

    inline void insert_lane(vComplex& z, int lane, std::complex<float> val) {
        insert_lane(z.re, lane, val.real());
        insert_lane(z.im, lane, val.imag());
    }

} // namespace simd

#endif // SIMD_TYPES_H