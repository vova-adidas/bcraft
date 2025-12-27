#pragma once

/**
 * @brief Scalar (non-SIMD) implementation of a Warp.
 *
 * `Warp1` is a simple, scalar version of the Warp interface intended primarily
 * for testing, debugging, and validation of rasterization logic.
 *
 * @note
 * This implementation is **not suitable for real use**
 * due to lack of vectorization and low performance.
 */
struct Warp1 {

	using Int = int;
	using IntMask = int;
	using FloatMask = float;
	using Float = float;

	constexpr static inline int SIZE = 1;
	constexpr static inline int ALIGNMENT = 4;

	inline constexpr static Float scalar(float v) { return v; }
	inline constexpr static Float load(float* src) { return *src; }
	inline constexpr static Float fma(Float a, Float b, Float c) { return a * b + c; }
	inline constexpr static Float mul(Float a, Float b) { return a * b; }
	inline constexpr static Float div(Float a, Float b) { return a / b; }
	inline constexpr static Float reciprocal(Float a) { return 1.f / a; }
	inline constexpr static Float add(Float a, Float b) { return a + b; }

	inline constexpr static Int scalar(int v) { return v; }
	inline constexpr static Int load(int* v) { return *v; }
	inline constexpr static Int mul(Int a, Int b) { return a * b; }
	inline constexpr static Int add(Int a, Int b) { return a + b; }

	inline constexpr static IntMask and_(IntMask a, IntMask b) { return (a == 1) && (b == 1); }
	inline constexpr static IntMask greater(Int a, Int b) { return a > b ? 1 : 0; }
};