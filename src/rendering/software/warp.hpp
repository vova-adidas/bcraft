#pragma once

#include <xsimd/xsimd.hpp>

struct Warp {

    using arch = xsimd::avx2;

    using Int = xsimd::batch<int, arch>;
    using IntMask = xsimd::batch_bool<int, arch>;
    using FloatMask = xsimd::batch_bool<float, arch>;
    using Float = xsimd::batch<float, arch>;

    constexpr static int SIZE = Int::size;
    constexpr static int ALIGNMENT = SIZE * 4;

    static Float scalar(float v) { return Float(v); }
    static Float load(float* src) { return Float::load_aligned(src); }
    static Float fma(Float a, Float b, Float c) { return xsimd::fma(a, b, c); }
    static Float mul(Float a, Float b) { return a * b; }
    static Float div(Float a, Float b) { return a / b; }
    static Float add(Float a, Float b) { return a + b; }
    static Float reciprocal(Float a) { return xsimd::reciprocal(a); }
    static Float floor(Float a) { return xsimd::floor(a); }
    static Float min(Float a, Float b) { return xsimd::min(a, b); }

    static Int scalar(int v) { return Int(v); }
    static Int load(int* v) { return Int::load_aligned(v); }
    static void store(Int a, int* v) { return a.store_aligned(v); }
    static Int mul(Int a, Int b) { return a * b; }
    static Int add(Int a, Int b) { return a + b; }
    static Int neg(Int a) { return -a; }
    static int reduce_max(Int a) { return xsimd::reduce_max(a); }
    static int reduce_min(Int a) { return xsimd::reduce_min(a); }
    static Int min(Int a, Int b) { return xsimd::min(a, b); }
    static Int max(Int a, Int b) { return xsimd::max(a, b); }

    static IntMask and_(IntMask a, IntMask b) { return a && b; }
    static IntMask andnot(IntMask a, IntMask b) { return xsimd::bitwise_andnot(a, b); }
    static FloatMask and_(FloatMask a, FloatMask b) { return a && b; }
    static IntMask greater(Int a, Int b) { return a > b; }
    static IntMask less(Int a, Int b) { return a < b; }
    static IntMask greater_equals(Int a, Int b) { return a >= b; }
    static IntMask less_equals(Int a, Int b) { return a <= b; }

    static FloatMask greater(Float a, Float b) { return a > b; }

    static Int select(IntMask mask, Int tr, Int fls) {
        return xsimd::select(mask, tr, fls);
    }

    static IntMask mask_from_bool(bool flag) {
        return IntMask(flag);
    }
};